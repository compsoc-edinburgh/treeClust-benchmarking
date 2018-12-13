# Script to test the impact of data size, outliers, noise and missing values on the detection of linear correlations by treeClust 
# 2018 G. Kustatscher, Rappsilber lab
# Note: Execution of the complete script will take several hours on a standard office computer

# Load the required libraries
library(data.table); library(treeClust); library(ggplot2); library(PRROC); library(WGCNA); library(viridis); library(grid); library(gridExtra)

# Set seed to make the results of this script reproducible
set.seed(123)

#### Create a set of necessary functions ####

## Function to create a synthetic dataset of experiments and proteins by randomly drawing values from a normal distribution
## The dataset is created in such a way that a defined percentage of protein pairs show a positive linear relationship, while the
## remaining pairs have no (= a random) relationship. The function accepts the following input parameters:
##        n_samples: Number of synthetic samples / measurements
##        n_proteins: Number of synthetic proteins
##        pc_related_pairs: Percentage of protein pairs that should have a positive linear relationship
##        pc_outliers: Percentage of samples that should constitute outlier measurements
##        pc_NA: Percentage of measurements that should be missing
##        j_amount: Amoung of measurement noise added to each protein

create_synthetic_dataset <- function(n_samples, n_proteins, pc_related_pairs, pc_outliers, pc_NA, j_amount){
  
  # Set the number of seed proteins ("chunk size") in which the dataset will be created
  seed_proteins <- 100
  
  # Throw error if n_proteins is not a multiple of the selected seed proteins
  if( n_proteins %% seed_proteins != 0 ){
    stop( "The value chosen for n_proteins is not valid. It needs to be a multiple of 100" )
  }
  
  # Create the first chunk of the synthetic dataset
  start_data <- replicate( seed_proteins, rnorm( n_samples ))   
  colnames(start_data) <- paste("start", 1:seed_proteins, sep = "_")

  # How many chunks of the same size need to be added to achieve user-specific dataset size?
  n_iter <- ( n_proteins/seed_proteins ) - 1
  
  # Add required chunks in iterations, by transforming the starting proteins so that each new protein with the same ID
  # is linearly correlated with its starting protein (with a certain amount of noise).
  data <- start_data
  for(i in 1:n_iter){
    iter_data <- apply( start_data, 2, function(x){ jitter( x , amount = j_amount ) })
    iter_index <- paste("iter", i, sep = "")
    colnames(iter_data) <- paste(iter_index, 1:seed_proteins, sep = "_")
    data <- cbind(data, iter_data)
    }
  
  # Assess pairwise associations created by these proteins
  n_total_pairs <- choose( n_proteins, 2 )                               # All pairwise combinations in the dataset
  n_related_pairs <- ((n_iter^2 + n_iter) / 2 ) * seed_proteins          # Pairwise combinations with a defined relationship
  n_rel_pai_req <- ceiling( n_total_pairs * ( pc_related_pairs / 100 ))  # Pairs with defined relationship that are required to fulfill user specification
  if( n_rel_pai_req >= n_related_pairs ) stop( "Cannot produce so many defined relationships with the current data dimensions. Try to adjust data dimensions" )
  
  # Calculate number of proteins per iteration / chunk that need to have a linear relationship, so that the user-specified number of defined pairs is achieved
  n_prots_per_iter <- ceiling( n_rel_pai_req / ((n_iter^2 + n_iter) / 2 ))
  if( n_prots_per_iter >= seed_proteins ) stop( "Cannot produce so many defined relationships with the current data dimensions. Try to adjust data dimensions." )
  if( n_prots_per_iter == 0 ) stop( "Cannot produce dataset with so few defined relationships")
  
  # For each iteration / chunk, replace values of surplus defined proteins by sampling from a normal distribution, reverting them to "start" proteins
  idcs_to_randomise_per_iter <- (n_prots_per_iter + 1):seed_proteins     # Protein indices per added iteration that need to be randomised
  idcs_to_randomise_dataset <- integer()                                 # Get the actual indices for "data" that correpond to these proteins
  for(i in 1:n_iter){ idcs_to_randomise_dataset <- c(idcs_to_randomise_dataset, (i*seed_proteins) + idcs_to_randomise_per_iter)}
  data[, idcs_to_randomise_dataset ] <- replicate( length(idcs_to_randomise_dataset) , rnorm( n_samples ))  # Replace these proteins with the same type of random values as used for the start set
  colnames(data)[ idcs_to_randomise_dataset ] <- paste("start",                                             # Adjust column names accordingly
                                                      (seed_proteins+1):(seed_proteins+length(idcs_to_randomise_dataset)),
                                                       sep = "_")
  
  # Insert a defined percentage of outliers
  n_outliers <-  ceiling( n_samples * pc_outliers/100 )                               # User-specified n of samples that should be turned into outliers
  if( n_outliers == 1) stop("Cannot introduce such few outliers. Choose more samples or less outliers.")
  if( n_outliers >= 2){
    data[ 1:n_outliers, ] <- apply( data[ 1:n_outliers, ], 2,                         # Create outliers by adding extra noise drawn from a larger normal distribution
                                    function(x){ x + rnorm(n = length(x), sd = 2.5) })  
  }
 
  # Insert a defined percentage of missing values
  n_values <- nrow(data) * ncol(data)      # Total n of data points
  n_NA_values <- n_values * pc_NA/100      # User-specified n of data points that should be NA
  pos_NA <- sample(n_values, n_NA_values)  # Randomly select positions in the data matrix to be turned into NA
  data[ pos_NA ] <- NA                     # And turn them to NA
  
  ## Test the data created and print stats
  # print( "A dataset with the following characteristics has been produced:" )
  # print( paste( nrow(data), "rows (samples)" ))
  # print( paste( ncol(data), "columns (proteins)" ))
  
  def_prot_per_iter <- max( as.integer( gsub("iter.+_", "", colnames(data)[ grep("iter", colnames(data)) ])))
  n_def_combis <- ((n_iter^2 + n_iter) / 2 ) * def_prot_per_iter
  pc_combis <- round( n_def_combis / choose( ncol(data), 2) * 100 , 2)
  # print( paste( pc_combis, "% of protein pairs expected to have positive linear relationship", sep = ""))
  warning_message <- paste("% defined pairwise relationships achieved was approximate. Setting:", pc_related_pairs, "Achieved:", pc_combis, sep = " ")
  if( pc_combis > pc_related_pairs*1.03 | pc_combis < pc_related_pairs*0.97 ) warning( warning_message )
  
  # print( paste( round( n_outliers/(n_samples/100), 2), "% samples are outliers" , sep = ""))
  # 
  # pc_NA_tested <- sum(is.na(data)) / (sum(is.na(data)) + sum(!is.na(data))) * 100 
  # print( paste( pc_NA_tested, "% missing values", sep = ""))
  # 
  # print( paste( j_amount, "noise setting" ))

  # Return the dataset
  return(data)
  }


## Functions to learn dissimilarities: treeClust, Pearson correlation (PCC), Spearman correlation (RHO) and bicor robust correlation
f_tC_dist <- function(x){
  temp_data <- as.data.frame( t( x ))                                # Transform matrix for treeClust learning
  temp_dist <- treeClust.dist(temp_data, d.num = 2, verbose = FALSE) # Learn treeClust dissimilarities
  tc_dist <- as.data.table( melt( as.matrix( temp_dist )))           # Turn dist object into a long data.table
  tc_dist[, tC_sim := 1-value ]                 # Turn treeClust dissimilarities into similarities
  tc_dist[, Var1 := as.character(Var1) ]        # Turn protein names into character vector (from factors)
  tc_dist[, Var2 := as.character(Var2) ]        # Turn protein names into character vector (from factors)
  tc_dist <- tc_dist[, .(Var1, Var2, tC_sim) ]  # Keep only relevant columns
  tc_dist <- tc_dist[ Var1 > Var2 ]             # Remove duplicates (incl. self-references)
  tc_dist[ gsub(".+_", "", Var1) == gsub(".+_", "", Var2),   # Select pairs derived from the same protein 
           Class := "defined"]                               # And label them as "defined" 
  tc_dist[ is.na(Class), Class := "random" ]    # Label remaining (i.e. undefined) pairs as random
  tc_dist }

f_cor_dist <- function(x, corType){
  temp_data <- stats::cor( x , method = corType, use = "pairwise.complete.obs")  # Get correlation matrix (using R's default stats function)
  cor_dist <- as.data.table( melt( temp_data ))   # Turn into long format data table 
  cor_dist[, Var1 := as.character(Var1) ]         # Turn protein names into character vector (from factors)
  cor_dist[, Var2 := as.character(Var2) ]         # Turn protein names into character vector (from factors)
  cor_dist <- cor_dist[ Var1 > Var2 ]             # Remove duplicates (incl. self-references)
  cor_dist[ gsub(".+_", "", Var1) == gsub(".+_", "", Var2),   # Select pairs derived from the same protein 
            Class := "defined"]                               # And label them as "defined" 
  cor_dist[ is.na(Class), Class := "random" ]     # Label remaining (i.e. undefined) pairs as random
  cor_dist }

f_bicor_dist <- function(x){
  temp_data <- bicor( x , use = "pairwise.complete.obs")  # Get robust correlation matrix
  cor_dist <- as.data.table( melt( temp_data ))   # Turn into long format data table 
  cor_dist[, Var1 := as.character(Var1) ]         # Turn protein names into character vector (from factors)
  cor_dist[, Var2 := as.character(Var2) ]         # Turn protein names into character vector (from factors)
  cor_dist <- cor_dist[ Var1 > Var2 ]             # Remove duplicates (incl. self-references)
  cor_dist[ gsub(".+_", "", Var1) == gsub(".+_", "", Var2),   # Select pairs derived from the same protein 
            Class := "defined"]                               # And label them as "defined" 
  cor_dist[ is.na(Class), Class := "random" ]     # Label remaining (i.e. undefined) pairs as random
  cor_dist }


# Plot formatting settings
my_plot_theme <-  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                        legend.position = "none")
my_cols <- c( bicor = "mediumblue", PCC = "steelblue2", RHO = "springgreen", treeClust = "magenta2" )


#### Test 1a: Number of samples ####

# Define which sample sizes to test
test_values <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = i, 
                                      n_proteins = 500,
                                      pc_related_pairs = 0.3,
                                      pc_outliers = 0,
                                      pc_NA = 0,
                                      j_amount = 0.2)
    
    # Learn distances / correlations
    TCL <- f_tC_dist(    sdata  )
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
  
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
  
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
  
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pN_samples <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
                geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 1, size = 0.25)+              
                geom_point( size = 0.2 )+
                geom_path( size = 0.25 )+
                scale_colour_manual( values = my_cols)+
                scale_x_continuous( limits = c(0, 50), breaks = seq(0, 50, 10))+
                scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
                xlab("Number of samples")+
                ylab("AUPRC")+
                my_plot_theme


#### Test 1b: Percentage of defined relationships ####

# Define which percentage of defined relationships to test
test_values <- seq(0.1, 0.5, 0.05)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 50, 
                                      n_proteins = 500,
                                      pc_related_pairs = i,
                                      pc_outliers = 0,
                                      pc_NA = 0,
                                      j_amount = 0.2)
    
    # Learn distances / correlations
    TCL <- f_tC_dist(    sdata  )
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pPct_rel <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
              geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 0.01, size = 0.25)+              
              geom_point( size = 0.2 )+
              geom_path( size = 0.25 )+
              scale_colour_manual( values = my_cols)+
              scale_x_continuous( limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1))+
              scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
              xlab("Co-expressed protein pairs [%]")+
              ylab("AUPRC")+
              my_plot_theme


#### Test 1c: Inter-dependence of sample number and pct of defined relationships (2D benchmark) ####

# To show that treeClust "learns" better if there are either more defined relationships or more samples
# I perform a "2D" test, where the AUPRC is shown through colour and both n_samples and pc_related_pairs are varied

# Define the range and combinations of values to test
test_values_n_samples <- seq(20,300,20)
test_values_pc_related_pairs <- seq(0.1, 0.5, 0.05)
test_values <- expand.grid(test_values_n_samples, test_values_pc_related_pairs)

# Initialise result table
results <- data.table()

# Get AUPRC in for loop
for(i in 1:nrow(test_values)){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = test_values[i, "Var1"], 
                                      n_proteins = 500,
                                      pc_related_pairs = test_values[i, "Var2"],
                                      pc_outliers = 0,
                                      pc_NA = 0,
                                      j_amount = 0.2)
    
    # Learn treeClust distances
    TCL <- f_tC_dist(    sdata  )
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
   
    results <- rbind(results, data.table(Var1 = test_values[i, "Var1"], Var2 = test_values[i, "Var2"], treeClust = TCL_auprc))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "------------------------------"))
    
}


## Plot the results
pTile1 <- ggplot(results, aes(x = Var2, y = Var1, fill = treeClust))+
           geom_tile()+
           scale_fill_viridis(option = "E", guide_legend(title="AUPRC", keywidth = unit(0.7, "cm"), keyheight = unit(1.4, "cm")),
                              limits = c(0,1))+
           xlab("Co-expressed protein pairs [%]")+
           ylab("Number of samples")+
           scale_x_continuous( breaks = seq(0.1, 1, 0.1), expand = c(0,0))+
           scale_y_continuous( breaks = seq(0,300,20)  , expand = c(0,0))+
           theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                 axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                 legend.position = c(0.6, 0.6), legend.background = element_rect(color = "black", fill = "white", 
                 size = 0.25, linetype = "solid"), legend.text = element_text(size=5), legend.title = element_text(size=6))


#### Test 1d: Number of proteins ####

# Define which sample sizes to test
test_values <- seq(200,1000,100)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 20, 
                                      n_proteins = i,
                                      pc_related_pairs = 0.3,
                                      pc_outliers = 0,
                                      pc_NA = 0,
                                      j_amount = 0.2)
    
    # Learn distances / correlations
    TCL <- f_tC_dist(    sdata  )
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pN_proteins <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
                geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 20, size = 0.25)+              
                geom_point( size = 0.2 )+
                geom_path( size = 0.25 )+
                scale_colour_manual( values = my_cols)+
                scale_x_continuous( limits = c(0, 1001), breaks = seq(0, 2000, 200))+
                scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
                xlab("Number of proteins")+
                ylab("AUPRC")+
                my_plot_theme


## Combine and print these top 4 figure panels

# Print on screen
p1 <- arrangeGrob(pN_samples, pPct_rel, pTile1, pN_proteins, nrow = 1)

# Clear workspace
rm( list = ls()[! ls() %in% c("create_synthetic_dataset", "f_tC_dist", "f_cor_dist", "f_bicor_dist", "my_plot_theme", "my_cols",
                              "p1")] )
gc()


#### Test 2: Amount of noise ####

# Define which sample sizes to test
test_values <- seq(0,1.2,0.2)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 50, 
                                      n_proteins = 500,
                                      pc_related_pairs = 0.5,
                                      pc_outliers = 0,
                                      pc_NA = 0,
                                      j_amount = i)
    
    # Learn distances / correlations
    TCL <- f_tC_dist(    sdata  )
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pNoise <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
            geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 0.05, size = 0.25)+              
            geom_point( size = 0.2 )+
            geom_path( size = 0.25 )+
            scale_colour_manual( values = my_cols)+
            scale_x_continuous( limits = c(0, 1.25), breaks = seq(0, 2, 0.2))+
            scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
            xlab("Amount of noise (jitter)")+
            ylab("AUPRC")+
            my_plot_theme


## Plot some example scatterplots to illustrate the increasing amount of noise

# Get three examples pairs
sdata_0.2_noise <- create_synthetic_dataset(j_amount = 0.2, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 0, pc_NA = 0)
sdata_0.6_noise <- create_synthetic_dataset(j_amount = 0.6, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 0, pc_NA = 0)
sdata_1.0_noise <- create_synthetic_dataset(j_amount = 1.0, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 0, pc_NA = 0)

# Create empty base plot
pB <- ggplot( data.frame(), aes(x = start_1, y = iter1_1))+
        scale_x_continuous( limits = c(-3, 3), breaks = seq(-2, 2, 2))+
        scale_y_continuous( limits = c(-3, 3), breaks = seq(-2, 2, 2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Populate plot with actual data
pNoise_1 <- pB + geom_point( data = as.data.frame(sdata_0.2_noise) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 0.2, size = 2.1)
pNoise_2 <- pB + geom_point( data = as.data.frame(sdata_0.6_noise) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 0.6, size = 2.1)
pNoise_3 <- pB + geom_point( data = as.data.frame(sdata_1.0_noise) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 1.0, size = 2.1)

# Combine scatterplots with AUPRC plots  
p2 <- arrangeGrob(pNoise_1, pNoise_2, pNoise_3, pNoise, nrow = 1)

# Clear workspace
rm( list = ls()[! ls() %in% c("create_synthetic_dataset", "f_tC_dist", "f_cor_dist", "f_bicor_dist", "my_plot_theme", "my_cols",
                              "p1", "p2")] )
gc()


#### Test 3a: Percentage of missing values (using 50 samples and 500 proteins) ####

# Define which sample sizes to test
test_values <- seq(0, 40, 5)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 50, 
                                      n_proteins = 500,
                                      pc_related_pairs = 0.5,
                                      pc_outliers = 0,
                                      pc_NA = i,
                                      j_amount = 0.2)
    
    # Learn distances / correlations (As missing values in the input data sometimes cause treeClust to fail, 
    # it is executed here in a while-tryCatch loop which will repeat until execution was successful
    TCL <- NULL
    while( is.null(TCL)) { TCL <- tryCatch( f_tC_dist(    sdata  ), error = function(e) NULL ) 
       if( is.null(TCL)) print("treeClust execution failed and was repeated (ignore `max(me)` warning)") } 
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Missing values in the input data interfere with PR calculation, so remove them with a message
    if( TCL[ !complete.cases(TCL) , .N] > 0 ) print( paste( TCL[ !complete.cases(TCL) , .N], "NAs removed from TCL" ))
    if( PCC[ !complete.cases(PCC) , .N] > 0 ) print( paste( PCC[ !complete.cases(PCC) , .N], "NAs removed from PCC" ))
    if( RHO[ !complete.cases(RHO) , .N] > 0 ) print( paste( RHO[ !complete.cases(RHO) , .N], "NAs removed from RHO" ))
    if( BIC[ !complete.cases(BIC) , .N] > 0 ) print( paste( BIC[ !complete.cases(BIC) , .N], "NAs removed from BIC" ))
    TCL <- TCL[ complete.cases(TCL) ]
    PCC <- PCC[ complete.cases(PCC) ]
    RHO <- RHO[ complete.cases(RHO) ]
    BIC <- BIC[ complete.cases(BIC) ]
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pNA_1 <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
          geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 1, size = 0.25)+              
          geom_point( size = 0.2 )+
          geom_path( size = 0.25 )+
          scale_colour_manual( values = my_cols)+
          scale_x_continuous( limits = c(0, 41), breaks = seq(0, 100, 10))+
          scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
          xlab("Missing values [%]")+
          ylab("AUPRC")+
          annotate("text", x = 10, y = 0.2, label = "50 samples, 500 proteins", size = 2)+
          my_plot_theme


#### Test 3b: Percentage of missing values (using 200 samples and 2000 proteins) ####

# Define which sample sizes to test
test_values <- seq(0, 40, 5)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 100, 
                                      n_proteins = 1000,
                                      pc_related_pairs = 0.5,
                                      pc_outliers = 0,
                                      pc_NA = i,
                                      j_amount = 0.2)
    
    # Learn distances / correlations (As missing values in the input data sometimes cause treeClust to fail, 
    # it is executed here in a while-tryCatch loop which will repeat until execution was successful
    TCL <- NULL
    while( is.null(TCL)) { TCL <- tryCatch( f_tC_dist(    sdata  ), error = function(e) NULL ) 
       if( is.null(TCL)) print("treeClust execution failed and was repeated (ignore `max(me)` warning)") } 
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Missing values in the input data interfere with PR calculation, so remove them with a message
    if( TCL[ !complete.cases(TCL) , .N] > 0 ) print( paste( TCL[ !complete.cases(TCL) , .N], "NAs removed from TCL" ))
    if( PCC[ !complete.cases(PCC) , .N] > 0 ) print( paste( PCC[ !complete.cases(PCC) , .N], "NAs removed from PCC" ))
    if( RHO[ !complete.cases(RHO) , .N] > 0 ) print( paste( RHO[ !complete.cases(RHO) , .N], "NAs removed from RHO" ))
    if( BIC[ !complete.cases(BIC) , .N] > 0 ) print( paste( BIC[ !complete.cases(BIC) , .N], "NAs removed from BIC" ))
    TCL <- TCL[ complete.cases(TCL) ]
    PCC <- PCC[ complete.cases(PCC) ]
    RHO <- RHO[ complete.cases(RHO) ]
    BIC <- BIC[ complete.cases(BIC) ]
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pNA_2 <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
          geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 1, size = 0.25)+              
          geom_point( size = 0.2 )+
          geom_path( size = 0.25 )+
          scale_colour_manual( values = my_cols)+
          scale_x_continuous( limits = c(0, 41), breaks = seq(0, 100, 10))+
          scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
          xlab("Missing values [%]")+
          ylab("AUPRC")+
          annotate("text", x = 10, y = 0.2, label = "100 samples and 1000 proteins", size = 2)+
          my_plot_theme


#### Test 3c: Inter-dependence of sample number and pct of missing values (2D benchmark) ####

# Define the range and combinations of values to test
test_values_n_samples <- seq(50,400,50)
test_values_pc_NA <- seq(10,40,5)
test_values <- expand.grid(test_values_n_samples, test_values_pc_NA)

# Initialise result table
results <- data.table()

# Get AUPRC in for loop
for(i in 1:nrow(test_values)){
  
  # Create synthetic dataset
  sdata <- create_synthetic_dataset(n_samples = test_values[i, "Var1"], 
                                    n_proteins = 1000,
                                    pc_related_pairs = 0.5,
                                    pc_outliers = 0,
                                    pc_NA = test_values[i, "Var2"],
                                    j_amount = 0.2)
  
  # Learn distances / correlations (As missing values in the input data sometimes cause treeClust to fail, 
  # it is executed here in a while-tryCatch loop which will repeat until execution was successful
  TCL <- NULL
  while( is.null(TCL)) { TCL <- tryCatch( f_tC_dist(    sdata  ), error = function(e) NULL ) 
     if( is.null(TCL)) print("treeClust execution failed and was repeated (ignore `max(me)` warning)") } 
  
  # Get areas under the corresponding PR curves
  TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
  
  results <- rbind(results, data.table(Var1 = test_values[i, "Var1"], Var2 = test_values[i, "Var2"], treeClust = TCL_auprc))
  
  # Print a line to make it easier to read which messages belong to which iteration
  print( paste("--------- This was i =", i, "------------------------------"))
  
}


## Plot the results
pTileNA <- ggplot(results, aes(x = Var2, y = Var1, fill = treeClust))+
            geom_tile()+
            scale_fill_viridis(option = "E", guide_legend(title="AUPRC", keywidth = unit(0.7, "cm"), keyheight = unit(1.4, "cm")),
                               limits = c(0,1))+
            xlab("Missing values [%]")+
            ylab("Number of samples")+
            scale_x_continuous( breaks = seq(0,100,5), expand = c(0,0))+
            scale_y_continuous( breaks = seq(0,500,100), expand = c(0,0))+
            theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                  axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                  legend.position = c(0.6, 0.6), legend.background = element_rect(color = "black", fill = "white", size = 0.25, linetype = "solid"), 
                  legend.text = element_text(size=5), legend.title = element_text(size=6))


#### Test 3d: Inter-dependence of protein number and pct of missing values (2D benchmark) ####

# Set seed
set.seed(123)

# Define the range and combinations of values to test
test_values_n_proteins <- seq(300,2400,300)
test_values_pc_NA <- seq(10,40,5)
test_values <- expand.grid(test_values_n_proteins, test_values_pc_NA)

# Initialise result table
results <- data.table()

# Get AUPRC in for loop
for(i in 1:nrow(test_values)){
  
  # Create synthetic dataset
  sdata <- create_synthetic_dataset(n_samples = 150, 
                                    n_proteins = test_values[i, "Var1"],
                                    pc_related_pairs = 0.5,
                                    pc_outliers = 0,
                                    pc_NA = test_values[i, "Var2"],
                                    j_amount = 0.2)
  
  # Learn distances / correlations (As missing values in the input data sometimes cause treeClust to fail, 
  # it is executed here in a while-tryCatch loop which will repeat until execution was successful
  TCL <- NULL
  while( is.null(TCL)) { TCL <- tryCatch( f_tC_dist(    sdata  ), error = function(e) NULL ) 
     if( is.null(TCL)) print("treeClust execution failed and was repeated (ignore `max(me)` warning)") } 
  
  # Get areas under the corresponding PR curves
  TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
  
  results <- rbind(results, data.table(Var1 = test_values[i, "Var1"], Var2 = test_values[i, "Var2"], treeClust = TCL_auprc))
  
  # Print a line to make it easier to read which messages belong to which iteration
  print( paste("--------- This was i =", i, "------------------------------"))
  
}


## Plot the results
pTileNA2 <- ggplot(results, aes(x = Var2, y = Var1, fill = treeClust))+
            geom_tile()+
            scale_fill_viridis(option = "E", guide_legend(title="AUPRC", keywidth = unit(0.7, "cm"), keyheight = unit(1.4, "cm")),
                               limits = c(0,1))+
            xlab("Missing values [%]")+
            ylab("Number of proteins")+
            scale_x_continuous( breaks = seq(0,100,5), expand = c(0,0))+
            scale_y_continuous( breaks = seq(0,3000,300)  , expand = c(0,0))+
            theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                  axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                  legend.position = c(0.6, 0.6), legend.background = element_rect(color = "black", fill = "white", size = 0.25, linetype = "solid"), 
                  legend.text = element_text(size=5), legend.title = element_text(size=6))


## Combine the plots of this part
p3 <- arrangeGrob(pNA_1, pNA_2, pTileNA, pTileNA2, nrow = 1)

# Clear workspace
rm( list = ls()[! ls() %in% c("create_synthetic_dataset", "f_tC_dist", "f_cor_dist", "f_bicor_dist", "my_plot_theme", "my_cols",
                              "p1", "p2", "p3")] )
gc()


#### Test 4: Percentage of outliers ####

# Define which sample sizes to test
test_values <- seq(0, 90, 5)

# Initialise result table
results <- data.table()

# Collect three replicate measurements each
for(k in 1:3){
  
  # Get AUPRC in for loop
  for(i in test_values){
    
    # Create synthetic dataset
    sdata <- create_synthetic_dataset(n_samples = 100, 
                                      n_proteins = 500,
                                      pc_related_pairs = 0.5,
                                      pc_outliers = i,
                                      pc_NA = 0,
                                      j_amount = 0.2)
    
    # Learn distances / correlations
    TCL <- f_tC_dist(    sdata  )
    PCC <- f_cor_dist(   sdata  , "pearson")
    RHO <- f_cor_dist(   sdata  , "spearman")
    BIC <- f_bicor_dist( sdata  )
    
    # Get areas under the corresponding PR curves
    TCL_auprc <- pr.curve( TCL[ Class == "defined", tC_sim] , TCL[ Class == "random", tC_sim] )$auc.integral
    PCC_auprc <- pr.curve( PCC[ Class == "defined", value]  , PCC[ Class == "random", value]  )$auc.integral
    RHO_auprc <- pr.curve( RHO[ Class == "defined", value]  , RHO[ Class == "random", value]  )$auc.integral
    BIC_auprc <- pr.curve( BIC[ Class == "defined", value]  , BIC[ Class == "random", value]  )$auc.integral
    
    results <- rbind(results,
                     data.table(test_values = i, replicate = k, treeClust = TCL_auprc, PCC = PCC_auprc, RHO = RHO_auprc, bicor = BIC_auprc ))
    
    # Print a line to make it easier to read which messages belong to which iteration
    print( paste("--------- This was i =", i, "(replicate =", k, ") -----------------------"))
    
  }
}

# Calculate the means and SEMs
res1 <- melt( results[, lapply(.SD, mean),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
res2 <- melt( results[, lapply(.SD, function(x){ sd(x) / sqrt(length(x)) }),
                      by = test_values , .SDcols = c("treeClust", "PCC", "RHO", "bicor") ] , id.vars = "test_values")
results <- merge(res1, res2[, .(test_values, variable, sem = value)] )


## Plot the results
pOutlier <- ggplot(results, aes(x = test_values, y = value , colour = variable))+
            geom_errorbar( aes( ymin = value - sem, ymax = value + sem), width = 1, size = 0.25)+              
            geom_point( size = 0.2 )+
            geom_path( size = 0.25 )+
            scale_colour_manual( values = my_cols)+
            scale_x_continuous( breaks = seq(0, 90, 10))+
            scale_y_continuous( limits = c(-0.02, 1.02), breaks = seq(0,1,0.2))+
            xlab("Outlier data [%]")+
            ylab("AUPRC")+
            my_plot_theme


## Plot some example scatterplots to illustrate the increasing amount of noise

# Get three examples pairs
sdata_00_out <- create_synthetic_dataset(j_amount = 0.2, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 0, pc_NA = 0)
sdata_30_out <- create_synthetic_dataset(j_amount = 0.2, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 30, pc_NA = 0)
sdata_60_out <- create_synthetic_dataset(j_amount = 0.2, n_samples = 100, n_proteins = 500, pc_related_pairs = 0.5, pc_outliers = 60, pc_NA = 0)

# Create empty base plot
pB <- ggplot( data.frame(), aes(x = start_1, y = iter1_1))+
      scale_x_continuous( limits = c(-6, 6), breaks = seq(-6, 6, 2))+
      scale_y_continuous( limits = c(-6, 6), breaks = seq(-6, 6, 2))+
      xlab("Ratio protein 1")+
      ylab("Ratio protein 2")+
      my_plot_theme

# Populate plot with actual data
pOut_1 <- pB + geom_point( data = as.data.frame(sdata_00_out) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 00, size = 2.1)
pOut_2 <- pB + geom_point( data = as.data.frame(sdata_30_out) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 20, size = 2.1)
pOut_3 <- pB + geom_point( data = as.data.frame(sdata_60_out) , size = 0.1 ) + annotate("text", x = -2, y = 2, label = 40, size = 2.1)

# Combine scatterplots with AUPRC plots  
p4 <- arrangeGrob(pOut_1, pOut_2, pOut_3, pOutlier, nrow = 1)

# Clear workspace
rm( list = ls()[! ls() %in% c("create_synthetic_dataset", "f_tC_dist", "f_cor_dist", "f_bicor_dist", "my_plot_theme", "my_cols",
                              "p1", "p2", "p3", "p4")] )
gc()


#### Create combined output plot ####

# Create the plot
pAll <- arrangeGrob(p1, p2, p3, p4, nrow = 4)
grid.newpage()
grid.draw(pAll)

# Save the plot
ggsave("Benchmark_plot.pdf", pAll, width = 18, height = 20, units = "cm")

