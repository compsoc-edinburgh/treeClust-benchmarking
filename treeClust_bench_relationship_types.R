# Script to test which type of relationships treeClust can detect
# 2018 G. Kustatscher, Rappsilber lab

# Load the required libraries
library(data.table); library(treeClust); library(ggplot2); library(PRROC); library(WGCNA); library(grid); library(gridExtra)

# Set seed to make the results of this script reproducible
set.seed(123)

#### Illustrate Anscombe's quartet ####

# Plot formatting settings
my_plot_theme <-  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                        legend.position = "none")

# Get R - inbuilt Anscombe data in proper shape
dt_ans <- rbind( data.table( x = anscombe$x1, y = anscombe$y1 , n_facet = 1),
                 data.table( x = anscombe$x2, y = anscombe$y2 , n_facet = 2),
                 data.table( x = anscombe$x3, y = anscombe$y3 , n_facet = 3),
                 data.table( x = anscombe$x4, y = anscombe$y4 , n_facet = 4),
                 data.table( x = NA,          y = NA ,          n_facet = 5))  # The last is an empty facet just to ensure plotting dimensions are the same as other plots in this script

# Get correlation values
ans_cors_1 <- data.table( round( cor( anscombe$x1, anscombe$y1 , method = "pearson"  ), 2),
                          round( cor( anscombe$x1, anscombe$y1 , method = "spearman" ), 2),
                          round( bicor( anscombe$x1, anscombe$y1                     ), 2),
                          1)

ans_cors_2 <- data.table( round( cor( anscombe$x2, anscombe$y2 , method = "pearson"  ), 2),
                          round( cor( anscombe$x2, anscombe$y2 , method = "spearman" ), 2),
                          round( bicor( anscombe$x2, anscombe$y2                     ), 2),
                          2)

ans_cors_3 <- data.table( round( cor( anscombe$x3, anscombe$y3 , method = "pearson"  ), 2),
                          round( cor( anscombe$x3, anscombe$y3 , method = "spearman" ), 2),
                          round( bicor( anscombe$x3, anscombe$y3                     ), 2),
                          3)

ans_cors_4 <- data.table( round( cor( anscombe$x4, anscombe$y4 , method = "pearson"  ), 2),
                          round( cor( anscombe$x4, anscombe$y4 , method = "spearman" ), 2),
                          "n.d.",
                          4)

names(ans_cors_1) <- c("PCC", "rho", "bicor", "n_facet")
names(ans_cors_2) <- c("PCC", "rho", "bicor", "n_facet")
names(ans_cors_3) <- c("PCC", "rho", "bicor", "n_facet")
names(ans_cors_4) <- c("PCC", "rho", "bicor", "n_facet")

ans_cors <- rbind(ans_cors_1, ans_cors_2, ans_cors_3, ans_cors_4)

# Create the plot
pAns <- ggplot( dt_ans, aes(x, y))+
        geom_smooth( method = "lm", se = FALSE, size = 0.25, colour = "orange")+
        geom_point( size = 0.1 )+
        geom_text( data = ans_cors[, .(PCC, n_facet)], aes( label = paste("PCC =", PCC), x = 2, y = 2), size = 2, hjust = 0, colour = "steelblue2")+
        geom_text( data = ans_cors[, .(rho, n_facet)], aes( label = paste("rho =", rho), x = 2, y = 1.6), size = 2, hjust = 0, colour = "springgreen")+
        geom_text( data = ans_cors[, .(bicor, n_facet)], aes( label = paste("bicor =", bicor), x = 2, y = 1.2), size = 2, hjust = 0, colour = "mediumblue")+
        facet_wrap(~n_facet, nrow = 1, scales = "free_y")+
        scale_x_continuous( limits = c(0, 20), breaks = seq(0,20,5))+
        scale_y_continuous( limits = c(0, 13), breaks = seq(0,14, 2))+
        xlab("x - values")+
        ylab("y - values")+
        my_plot_theme + theme( strip.background = element_blank(), strip.text = element_blank())


#### Define functions to learn dissimilarities: treeClust, Pearson correlation (PCC), Spearman correlation (RHO) and bicor robust correlation ####
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
  print( paste( round( tc_dist[ Class == "defined" , .N] / tc_dist[, .N] * 100, 2), "% defined pairs" ))
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
  print( paste( round( cor_dist[ Class == "defined" , .N] / cor_dist[, .N] * 100, 2), "% defined pairs" ))
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
  print( paste( round( cor_dist[ Class == "defined" , .N] / cor_dist[, .N] * 100, 2), "% defined pairs" ))
  cor_dist }


#### Create a range of synthetic gene expression datasets spiked with linear or various non-linear pairwise relationships #### 

# Note that the treeClust_benchmarking.R script contains a function that - specifically for positive linear correlations -
# can design more complex synthetic datasets (incl. outliers, missing values etc)

# Set the parameters 
n_samples <- 100        # Number of synthetic experiments
n_proteins <- 100       # Number of synthetic proteins
n_test <- 100           # Number of test proteins with defined (non-random) relationships to spike into the random dataset
j_amount <- 0.2         # Amount of jitter (noise) to be added to the test proteins

# Create a synthetic dataset of random experiments and proteins (no defined relationships between these proteins)
data <- replicate( n_proteins, rnorm( n_samples ))       # Expression values are random samples of normal distributions (mean = 0, sd = 1)
colnames(data) <- paste("ran", 1:n_proteins, sep = "_")  # Name the random proteins

# Add proteins with a *positive linear* relationship to existing ones
data_pos_linear <- apply( data[, 1:n_test] , 2, function(x){ jitter( x , amount = j_amount ) })  # Create n_test synthetic proteins with a defined relationship to the first n_test proteins in the dataset
colnames(data_pos_linear) <- paste("a_pos_lin", 1:n_test, sep = "_")                             # Name the proteins with the defined relationship
data_pos_linear <- cbind(data, data_pos_linear)                                                  # Add the proteins with the defined relationship to the dataset
                                                                                                 # E.g. protein pos_lin_1 has a positive linear relationship with ran_1, pos_lin_2 with ran_2, etc.
                                                                                                 # In total, there are n_test defined pairwise relationships in the dataset

# Add proteins with a *negative linear* relationship to existing ones
data_neg_linear <- apply( data[, 1:n_test] , 2, function(x){ jitter( -x , amount = j_amount ) })        # See explanation above
colnames(data_neg_linear) <- paste("b_neg_lin", 1:n_test, sep = "_")
data_neg_linear <- cbind(data, data_neg_linear) 

# Add proteins with an *exponential* relationship to existing ones
data_exponential <- apply( data[, 1:n_test] , 2, function(x){ jitter(exp(1)^x , amount = j_amount ) })  # See explanation above
colnames(data_exponential) <- paste("c_exp", 1:n_test, sep = "_")
data_exponential <- cbind(data, data_exponential) 

# Add proteins with an *sigmoid* relationship to existing ones (using the logistic function)
data_sigmoid <- apply( data[, 1:n_test] , 2, function(x){ jitter( (4 / (1 + exp(1)^(-5*(x-0)))) , amount = j_amount)  })  # See explanation above
colnames(data_sigmoid) <- paste("d_sig" , 1:n_test, sep = "_")
data_sigmoid <- cbind(data, data_sigmoid) 

# Add proteins with a *quadratic* relationship to existing ones
data_quadratic <- apply( data[, 1:n_test] , 2, function(x){ jitter( x^2 , amount = j_amount ) })        # See explanation above
colnames(data_quadratic) <- paste("e_qua", 1:n_test, sep = "_")
data_quadratic <- cbind(data, data_quadratic) 


#### Apply dissimilarity functions to the synthetic datasets ####

# Apply treeClust learning to the datasets spiked with various types of defined relationships
TCL_pos_linear  <- f_tC_dist( data_pos_linear  )
TCL_neg_linear  <- f_tC_dist( data_neg_linear  )
TCL_exponential <- f_tC_dist( data_exponential )
TCL_sigmoid     <- f_tC_dist( data_sigmoid     )
TCL_quadratic   <- f_tC_dist( data_quadratic   )

# Apply Pearson correlation to the datasets spiked with various types of defined relationships
PCC_pos_linear  <- f_cor_dist( data_pos_linear  , "pearson")
PCC_neg_linear  <- f_cor_dist( data_neg_linear  , "pearson")
PCC_exponential <- f_cor_dist( data_exponential , "pearson")
PCC_sigmoid     <- f_cor_dist( data_sigmoid     , "pearson")
PCC_quadratic   <- f_cor_dist( data_quadratic   , "pearson")

# Apply Spearman correlation to the datasets spiked with various types of defined relationships
RHO_pos_linear  <- f_cor_dist( data_pos_linear  , "spearman")
RHO_neg_linear  <- f_cor_dist( data_neg_linear  , "spearman")
RHO_exponential <- f_cor_dist( data_exponential , "spearman")
RHO_sigmoid     <- f_cor_dist( data_sigmoid     , "spearman")
RHO_quadratic   <- f_cor_dist( data_quadratic   , "spearman")

# Apply bicor correlation to the datasets spiked with various types of defined relationships
BIC_pos_linear  <- f_bicor_dist( data_pos_linear  )
BIC_neg_linear  <- f_bicor_dist( data_neg_linear  )
BIC_exponential <- f_bicor_dist( data_exponential )
BIC_sigmoid     <- f_bicor_dist( data_sigmoid     )
BIC_quadratic   <- f_bicor_dist( data_quadratic   )


#### Get the areas under the corresponding PR curves (AUPRC) ####

# For treeClust
TCL_pos_linear_pr <- pr.curve( TCL_pos_linear[ Class == "defined", tC_sim] , TCL_pos_linear[ Class == "random", tC_sim], curve = TRUE)
TCL_neg_linear_pr <- pr.curve( TCL_neg_linear[ Class == "defined", 1-tC_sim] , TCL_neg_linear[ Class == "random", 1-tC_sim], curve = TRUE)    # Invert the scores for negative correlation (otherwise PR plot won't work)
TCL_exponential_pr <- pr.curve( TCL_exponential[ Class == "defined", tC_sim] , TCL_exponential[ Class == "random", tC_sim], curve = TRUE)
TCL_sigmoid_pr <- pr.curve( TCL_sigmoid[ Class == "defined", tC_sim] , TCL_sigmoid[ Class == "random", tC_sim], curve = TRUE)
TCL_quadratic_pr <- pr.curve( TCL_quadratic[ Class == "defined", tC_sim] , TCL_quadratic[ Class == "random", tC_sim], curve = TRUE)

# For Pearson correlation
PCC_pos_linear_pr <- pr.curve( PCC_pos_linear[ Class == "defined", value] , PCC_pos_linear[ Class == "random", value], curve = TRUE)
PCC_neg_linear_pr <- pr.curve( PCC_neg_linear[ Class == "defined", 1-value] , PCC_neg_linear[ Class == "random", 1-value], curve = TRUE)
PCC_exponential_pr <- pr.curve( PCC_exponential[ Class == "defined", value] , PCC_exponential[ Class == "random", value], curve = TRUE)
PCC_sigmoid_pr <- pr.curve( PCC_sigmoid[ Class == "defined", value] , PCC_sigmoid[ Class == "random", value], curve = TRUE)
PCC_quadratic_pr <- pr.curve( PCC_quadratic[ Class == "defined", value] , PCC_quadratic[ Class == "random", value], curve = TRUE)

# For Spearman correlation
RHO_pos_linear_pr <- pr.curve( RHO_pos_linear[ Class == "defined", value] , RHO_pos_linear[ Class == "random", value], curve = TRUE)
RHO_neg_linear_pr <- pr.curve( RHO_neg_linear[ Class == "defined", 1-value] , RHO_neg_linear[ Class == "random", 1-value], curve = TRUE)
RHO_exponential_pr <- pr.curve( RHO_exponential[ Class == "defined", value] , RHO_exponential[ Class == "random", value], curve = TRUE)
RHO_sigmoid_pr <- pr.curve( RHO_sigmoid[ Class == "defined", value] , RHO_sigmoid[ Class == "random", value], curve = TRUE)
RHO_quadratic_pr <- pr.curve( RHO_quadratic[ Class == "defined", value] , RHO_quadratic[ Class == "random", value], curve = TRUE)

# For bicor correlation
BIC_pos_linear_pr <- pr.curve( BIC_pos_linear[ Class == "defined", value] , BIC_pos_linear[ Class == "random", value], curve = TRUE)
BIC_neg_linear_pr <- pr.curve( BIC_neg_linear[ Class == "defined", 1-value] , BIC_neg_linear[ Class == "random", 1-value], curve = TRUE)
BIC_exponential_pr <- pr.curve( BIC_exponential[ Class == "defined", value] , BIC_exponential[ Class == "random", value], curve = TRUE)
BIC_sigmoid_pr <- pr.curve( BIC_sigmoid[ Class == "defined", value] , BIC_sigmoid[ Class == "random", value], curve = TRUE)
BIC_quadratic_pr <- pr.curve( BIC_quadratic[ Class == "defined", value] , BIC_quadratic[ Class == "random", value], curve = TRUE)


#### Create a histogram for illustration ####

# For treeClust
pH_TCL <- ggplot( TCL_pos_linear, aes( x = tC_sim, fill = Class))+
            geom_histogram( bins = 100 )+
            scale_fill_manual( values = c("orange", "grey70"))+
            scale_y_continuous( limits = c(0, 2500), expand = c(0,0))+
            xlim(0,1)+
            xlab("treeClust similarity")+
            ylab("Protein pair count")+
            my_plot_theme

# For treeClust (zoom)
pH_TCL_Zoom <- ggplot( TCL_pos_linear, aes( x = tC_sim, fill = Class))+
                geom_histogram( bins = 100 )+
                scale_fill_manual( values = c("orange", "grey70"))+
                scale_y_continuous( limits = c(0, 2500), expand = c(0,0))+
                xlim(0,1)+
                coord_cartesian( ylim = c(0,100) )+
                xlab("treeClust similarity")+
                ylab("Protein pair count")+
                my_plot_theme

# For PCC
pH_PCC <- ggplot( PCC_pos_linear, aes( x = value, fill = Class))+
            geom_histogram( bins = 100 )+
            scale_fill_manual( values = c("orange", "grey70"))+
            scale_y_continuous( limits = c(0, 2000), expand = c(0,0))+
            xlim(-1,1)+
            xlab("Correlation coefficient")+
            ylab("Protein pair count")+
            my_plot_theme


#### Create example scatterplots for all types of pairwise relationships ####

# Random
pRd <- ggplot( as.data.frame(data_pos_linear), aes( ran_1, ran_2 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Positive linear
pPL <- ggplot( as.data.frame(data_pos_linear), aes( ran_1, a_pos_lin_1 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Negative linear
pNL <- ggplot( as.data.frame(data_neg_linear), aes( ran_1, b_neg_lin_1 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Exponential
pEx <- ggplot( as.data.frame(data_exponential), aes( ran_1, c_exp_1 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-0.5,9), breaks = seq(0,10,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Sigmoid
pSi <- ggplot( as.data.frame(data_sigmoid), aes( ran_1, d_sig_1 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-0.5,4.5), breaks = seq(-2,4,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme

# Quadratic
pQu <- ggplot( as.data.frame(data_quadratic), aes( ran_1, e_qua_1 ))+
        geom_point( size = 0.1 )+
        scale_x_continuous( limits = c(-2.5,2.5), breaks = seq(-2,2,2))+
        scale_y_continuous( limits = c(-0.5,6), breaks = seq(-2,7,2))+
        xlab("Ratio protein 1")+
        ylab("Ratio protein 2")+
        my_plot_theme


## Arrange the example plots for define relationships
p_def_examples <- arrangeGrob(pPL, pNL, pEx, pSi, pQu, nrow = 1)
grid.newpage()
grid.draw(p_def_examples)


#### Create PR curves for each relationship ####

my_cols <- c( bicor = "mediumblue", PCC = "steelblue2", RHO = "springgreen", treeClust = "magenta2" )     # Plotting colours

# Positive linear
pos_lin <- rbind( data.table( TCL_pos_linear_pr$curve[,1:2], method = "treeClust" ),
                  data.table( PCC_pos_linear_pr$curve[,1:2], method = "PCC"       ),
                  data.table( RHO_pos_linear_pr$curve[,1:2], method = "RHO"       ),
                  data.table( BIC_pos_linear_pr$curve[,1:2], method = "bicor"     ))
names(pos_lin)[1:2] <- c("Recall", "Precision")
       
pPR_PL <- ggplot(pos_lin, aes(Recall, Precision, colour = method))+
            geom_path( size = 0.25 )+
            scale_colour_manual(values = my_cols)+
            scale_x_continuous( breaks = seq(0,1,0.2))+
            scale_y_continuous( breaks = seq(0,1,0.2))+
            annotate("text", label = format( round( TCL_pos_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.20, size = 2.1, colour = my_cols[ "treeClust"] )+
            annotate("text", label = format( round( PCC_pos_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.17, size = 2.1, colour = my_cols[ "PCC"])+
            annotate("text", label = format( round( RHO_pos_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.14, size = 2.1, colour = my_cols[ "RHO"])+
            annotate("text", label = format( round( BIC_pos_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.11, size = 2.1, colour = my_cols[ "bicor"])+
            my_plot_theme

# Negative linear
neg_lin <- rbind( data.table( TCL_neg_linear_pr$curve[,1:2], method = "treeClust" ),
                  data.table( PCC_neg_linear_pr$curve[,1:2], method = "PCC"       ),
                  data.table( RHO_neg_linear_pr$curve[,1:2], method = "RHO"       ),
                  data.table( BIC_neg_linear_pr$curve[,1:2], method = "bicor"     ))
names(neg_lin)[1:2] <- c("Recall", "Precision")

pPR_NL <- ggplot(neg_lin, aes(Recall, Precision, colour = method))+
            geom_path( size = 0.25 )+
            scale_colour_manual(values = my_cols)+
            scale_x_continuous( breaks = seq(0,1,0.2))+
            scale_y_continuous( breaks = seq(0,1,0.2))+
            annotate("text", label = format( round( TCL_neg_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.20, size = 2.1, colour = my_cols[ "treeClust"] )+
            annotate("text", label = format( round( PCC_neg_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.17, size = 2.1, colour = my_cols[ "PCC"])+
            annotate("text", label = format( round( RHO_neg_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.14, size = 2.1, colour = my_cols[ "RHO"])+
            annotate("text", label = format( round( BIC_neg_linear_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.11, size = 2.1, colour = my_cols[ "bicor"])+
            my_plot_theme

# Exponential
exponential <- rbind( data.table( TCL_exponential_pr$curve[,1:2], method = "treeClust" ),
                      data.table( PCC_exponential_pr$curve[,1:2], method = "PCC"       ),
                      data.table( RHO_exponential_pr$curve[,1:2], method = "RHO"       ),
                      data.table( BIC_exponential_pr$curve[,1:2], method = "bicor"     ))
names(exponential)[1:2] <- c("Recall", "Precision")

pPR_Ex <- ggplot(exponential, aes(Recall, Precision, colour = method))+
            geom_path( size = 0.25 )+
            scale_colour_manual(values = my_cols)+
            scale_x_continuous( breaks = seq(0,1,0.2))+
            scale_y_continuous( breaks = seq(0,1,0.2))+
            annotate("text", label = format( round( TCL_exponential_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.20, size = 2.1, colour = my_cols[ "treeClust"] )+
            annotate("text", label = format( round( PCC_exponential_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.17, size = 2.1, colour = my_cols[ "PCC"])+
            annotate("text", label = format( round( RHO_exponential_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.14, size = 2.1, colour = my_cols[ "RHO"])+
            annotate("text", label = format( round( BIC_exponential_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.11, size = 2.1, colour = my_cols[ "bicor"])+
            my_plot_theme

# Sigmoid
sigmoid <- rbind( data.table( TCL_sigmoid_pr$curve[,1:2], method = "treeClust" ),
                  data.table( PCC_sigmoid_pr$curve[,1:2], method = "PCC"       ),
                  data.table( RHO_sigmoid_pr$curve[,1:2], method = "RHO"       ),
                  data.table( BIC_sigmoid_pr$curve[,1:2], method = "bicor"     ))
names(sigmoid)[1:2] <- c("Recall", "Precision")

pPR_Si <- ggplot(sigmoid, aes(Recall, Precision, colour = method))+
            geom_path( size = 0.25 )+
            scale_colour_manual(values = my_cols)+
            scale_x_continuous( breaks = seq(0,1,0.2))+
            scale_y_continuous( breaks = seq(0,1,0.2))+
            annotate("text", label = format( round( TCL_sigmoid_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.20, size = 2.1, colour = my_cols[ "treeClust"] )+
            annotate("text", label = format( round( PCC_sigmoid_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.17, size = 2.1, colour = my_cols[ "PCC"])+
            annotate("text", label = format( round( RHO_sigmoid_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.14, size = 2.1, colour = my_cols[ "RHO"])+
            annotate("text", label = format( round( BIC_sigmoid_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.11, size = 2.1, colour = my_cols[ "bicor"])+
            my_plot_theme
  

# Quadratic
quadratic <- rbind( data.table( TCL_quadratic_pr$curve[,1:2], method = "treeClust" ),
                    data.table( PCC_quadratic_pr$curve[,1:2], method = "PCC"       ),
                    data.table( RHO_quadratic_pr$curve[,1:2], method = "RHO"       ),
                    data.table( BIC_quadratic_pr$curve[,1:2], method = "bicor"     ))
names(quadratic)[1:2] <- c("Recall", "Precision")

pPR_Qu <- ggplot(quadratic, aes(Recall, Precision, colour = method))+
            geom_path( size = 0.25 )+
            scale_colour_manual(values = my_cols)+
            scale_x_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
            scale_y_continuous( limits = c(0,1), breaks = seq(0,1,0.2))+
            annotate("text", label = format( round( TCL_quadratic_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.20, size = 2.1, colour = my_cols[ "treeClust"] )+
            annotate("text", label = format( round( PCC_quadratic_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.17, size = 2.1, colour = my_cols[ "PCC"])+
            annotate("text", label = format( round( RHO_quadratic_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.14, size = 2.1, colour = my_cols[ "RHO"])+
            annotate("text", label = format( round( BIC_quadratic_pr$auc.integral, 2), nsmall = 2), x = 0.1, y = 0.11, size = 2.1, colour = my_cols[ "bicor"])+
            my_plot_theme


## Arrange the PR plots
p_PR_plots <- arrangeGrob(pPR_PL, pPR_NL, pPR_Ex, pPR_Si, pPR_Qu, nrow = 1)
grid.newpage()
grid.draw(p_PR_plots)


#### Create composite figure ####

# Create the top part
p_illustration <- arrangeGrob( pRd, pPL, pH_PCC, pH_TCL, pH_TCL_Zoom, nrow = 1)

# Create the whole plot
p_composite <- arrangeGrob(pAns, p_illustration, p_def_examples, p_PR_plots, nrow = 4)

# Print it
grid.newpage()
grid.draw(p_composite)

# Save it
ggsave("treeClust_relationships_benchmark.pdf", p_composite, width = 18, height = 14.4, units = c( "cm" ))
