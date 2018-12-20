# Script to test impact of outliers, goodness-of-fit and non-linear relationships on ProHD data,
# Comparing treeClust with Pearson's, Spearman's and biweight midcorrelation
# 2018 P. Grabowski & G. Kustatscher, Rappsilber lab

# Load the required libraries
library(WGCNA); library(data.table); library(ggplot2); library(grid); library(gridExtra);
library(scales); library(treeClust); library(readxl); library(cowplot)

# Set seed to make execution reproducible
set.seed(42)

#### Read in and prepare the ProHD dataset #### 

# Read in ProteomeHD and pre-prep
ProHD <- read.csv("ProteomeHD_v1_1.csv")                        # Get the data
rownames(ProHD) <- ProHD$Majority_protein_IDs                   # Assign rownames
ProHD <- ProHD[, grep("Ratio", colnames(ProHD))]                # Extract SILAC ratio columns
feature_count <- apply(ProHD, 1, function(x){ sum(!is.na(x))})  # Count number of SILAC ratios per protein
ProHD <- ProHD[ feature_count >= 95 ,]                          # Discard proteins detected in fewer than 95 experiments
medians <- apply(ProHD, 2, median, na.rm=T)
ProHD <- t(apply(ProHD, 1, function(x) x - medians))


#### treeClust analysis ####

# Run treeClust
tc <- treeClust.dist( data.frame( ProHD),
                      d.num = 2,
                      verbose = FALSE,
                      rcontrol = rpart.control(cp = 0.105),
                      control = treeClust.control(serule = 1.8))

tc <- as.data.table( melt( as.matrix(tc)))     # Convert it to a long data table
tc <- tc[, .( Protein_1 = as.character(Var1),  # Re-name columns and turn dist into similarity
              Protein_2 = as.character(Var2), 
              tC_sim = (1-value) )]
tc <- tc[ Protein_1 > Protein_2 ]              # Remove self-comparisons and duplicate pairs
tc <- tc[ complete.cases(tc) ]                 # Remove any arising NAs


# Transpose the values for correlation metrics
ProHD <- t( ProHD )


#### Annotate with gold standard ####

# Read in Reactome gold standard
TPFP <- fread("Reactome_TP_FP.csv")                             # Read in gold standard
names(TPFP) <- c("SimpleID_1", "SimpleID_2", "Class")           # Change colnames for merging

# Merge treeClust scores with TP/FP annotation, to have one combined table for downstream processes
tc[, SimpleID_1 := gsub(";.+", "", Protein_1) ][, SimpleID_1 := gsub("-.+", "", SimpleID_1) ]  # Simplify protein IDs
tc[, SimpleID_2 := gsub(";.+", "", Protein_2) ][, SimpleID_2 := gsub("-.+", "", SimpleID_2) ]  # Simplify protein IDs
DT <- merge(tc, TPFP, by = c("SimpleID_1", "SimpleID_2"), all.x = TRUE)                          # Merge data with gold standard annotation

# Clear workspace and memory
rm( list = setdiff( ls(), c("DT", "ProHD") ))
gc()


#### treeClust vs PCC ####

# Calculate PCC
ProHD_cor <- stats::cor(ProHD, use = "pairwise.complete.obs", method = c("pearson"))  # Get correlation matrix
ProHD_cor <- as.data.table( melt( ProHD_cor))                                         # Convert it to a long data table
ProHD_cor <- ProHD_cor[, .( Protein_1 = as.character(Var1),
                        Protein_2 = as.character(Var2),
                        PCC = value ) ]
ProHD_cor <- ProHD_cor[ Protein_1 > Protein_2 ]                                       # Remove self-comparisons and duplicate pairs
ProHD_cor <- ProHD_cor[ complete.cases(ProHD_cor) ]                                     # Drop the (few) pairs which didn't yield a PCC

# Combine with treeClust scores
DT <- merge( DT, ProHD_cor, by = c("Protein_1", "Protein_2"))

# Clear the large ProHD_cor object and clear memory 
rm( "ProHD_cor" )
gc()

# Rank protein pairs by treeClust similarity or PCC
DT[, coreg_score_rank := frank( -tC_sim ) ]
DT[,         PCC_rank := frank( -PCC ) ]

# Define "high" and "low" thresholds as a function of `rank percentage`
high_scoring <- DT[,.N]/100 * 0.1      # The top-ranking 0.1% of pairs
low_scoring <- DT[,.N]/100 * 0.5       # Not in the top 0.5% of pairs 

# Define protein pairs that were either detected by treeClust or by PCC
trC_only <- DT[ coreg_score_rank < high_scoring & PCC_rank > low_scoring  ]
PCC_only <- DT[ coreg_score_rank > low_scoring  & PCC_rank < high_scoring ]

# Combine into one test set
trC_only[, group := "trC_only" ]
PCC_only[, group := "PCC_only" ]
trC_vs_PCC <- rbind( trC_only, PCC_only )


## Model fitting & outlier removal ##

pb <- txtProgressBar(min = 0, max = trC_vs_PCC[,.N], style = 3)  # Initiate progress bar 
res <- data.table()                                              # Initiate result table

for(i in 1:trC_vs_PCC[,.N] ){
  my_pair <- trC_vs_PCC[i, c(Protein_1, Protein_2)]            # Define protein pair to assess in current iteration
  temp_dt <- na.omit( as.data.table( ProHD[,  my_pair ] ))       # Define those protein pair's data        
  names( temp_dt ) <- c("A", "B")                              # Replace protein IDs as colnames
  temp_dt <- temp_dt[, lapply(.SD, rescale )]                  # Rescale for curve-fitting
  
  if( temp_dt[,.N] < 50 ){                                     # If there are fewer than 50 common data points, skip the pair                          
    res <- rbind(res,                                          # ... and rather than fitting a model just add NAs to the results
                 data.table( Protein_1 = my_pair[1],           # See below for the definition of these parameters
                             Protein_2 = my_pair[2],
                             lm_all_r2 = NA,
                             lm_noStuRes_r2 = NA,
                             lm_noMahDis_r2 = NA,
                             mae_lm_mod = NA,
                             RSS_lm_mod = NA,
                             RSS_lin_nls = NA,
                             RSS_lin_nls_noStuRes = NA,
                             RSS_lin_nls_noMahDis = NA,
                             RSS_exp_nls = NA,
                             RSS_exp_nls_noStuRes = NA, 
                             RSS_exp_nls_noMahDis = NA,
                             RSS_sig_nls = NA,
                             RSS_sig_nls_noStuRes = NA,
                             RSS_sig_nls_noMahDis = NA,
                             PCC_all = NA,
                             PCC_noStuRes = NA,
                             PCC_noMahDis = NA,
                             rho_all = NA,
                             rho_noStuRes = NA,
                             rho_noMahDis = NA,
                             bic_all = NA,
                             bic_noStuRes = NA,
                             bic_noMahDis = NA,
                             N_StuRes_outlier = NA,
                             N_MahDis_outlier = NA,
                             N = temp_dt[,.N] ))
    
  } else {                                                    # If there are more than 50 common data, proceed as follows
    
    ## Fit linear regression model, identify outliers (via studentized residuals or mahalanobis distance), re-fit model without outliers
    
    lm_mod <- lm(A ~ B, data = temp_dt)                       # Fit a simple linear regression model
    temp_dt$StuRes <- rstudent(lm_mod)                        # Extract the studentized residuals of the model
    temp_dt$MahDis <- mahalanobis( temp_dt[,.(A,B)],          # Calculate the Mahalanobis distance for each point
                                   center = colMeans( temp_dt[,.(A,B)] ),
                                   cov( temp_dt[,.(A,B)] ))
    
    temp_dt[, StuRes_outlier := abs( StuRes ) > 2 ]           # Mark data points that are outliers according to studentized residuals
    temp_dt[, MahDis_outlier := MahDis > 2 ]                  # Mark data points that are outliers according to the Mahalanobis distance
    
    lm_noStuRes <- lm(A ~ B, data = temp_dt[ StuRes_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    lm_noMahDis <- lm(A ~ B, data = temp_dt[ MahDis_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    
    
    ## Fit non-linear models (exponential and sigmoid, as well as linear for comparison), output NaN if fitting does not converge
    
    # For all data
    lin_nls <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt, start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without StuRes outliers
    lin_nls_noStuRes <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noStuRes <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noStuRes <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without MahDis outliers
    lin_nls_noMahDis <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noMahDis <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noMahDis <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # Extract residual sum of squares (RSS) for pairs where model was built successfully
    RSS_lin_nls          <- if( class(lin_nls) == "nls" ) {          lin_nls$m$deviance()          } else { NaN }
    RSS_lin_nls_noStuRes <- if( class(lin_nls_noStuRes) == "nls" ) { lin_nls_noStuRes$m$deviance() } else { NaN }
    RSS_lin_nls_noMahDis <- if( class(lin_nls_noMahDis) == "nls" ) { lin_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_exp_nls          <- if( class(exp_nls) == "nls" ) {          exp_nls$m$deviance()          } else { NaN }
    RSS_exp_nls_noStuRes <- if( class(exp_nls_noStuRes) == "nls" ) { exp_nls_noStuRes$m$deviance() } else { NaN }
    RSS_exp_nls_noMahDis <- if( class(exp_nls_noMahDis) == "nls" ) { exp_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_sig_nls          <- if( class(sig_nls) == "nls" ) {          sig_nls$m$deviance()          } else { NaN }
    RSS_sig_nls_noStuRes <- if( class(sig_nls_noStuRes) == "nls" ) { sig_nls_noStuRes$m$deviance() } else { NaN }
    RSS_sig_nls_noMahDis <- if( class(sig_nls_noMahDis) == "nls" ) { sig_nls_noMahDis$m$deviance() } else { NaN }
    
    
    ## Collect results in a data.table
    res <- rbind(res,                                                            
                 data.table( Protein_1 = my_pair[1],                               # Name of protein 1
                             Protein_2 = my_pair[2],                               # Name of protein 2
                             lm_all_r2 =      summary( lm_mod      )$r.squared,    # R2 of standard linear model
                             lm_noStuRes_r2 = summary( lm_noStuRes )$r.squared,    # R2 of linear model without StuRes outliers
                             lm_noMahDis_r2 = summary( lm_noMahDis )$r.squared,    # R2 of linear model without MahDis outliers
                             mae_lm_mod = mean(abs(lm_mod$residuals)),             # Mean absolute error of the standard linear model
                             RSS_lm_mod = sum( lm_mod$residuals^2 ),               # Residual sum of squares of standard linear model
                             RSS_lin_nls,                                          # Residual sum of squares from the various non-linear fittings
                             RSS_lin_nls_noStuRes,
                             RSS_lin_nls_noMahDis,
                             RSS_exp_nls,
                             RSS_exp_nls_noStuRes, 
                             RSS_exp_nls_noMahDis,
                             RSS_sig_nls,
                             RSS_sig_nls_noStuRes,
                             RSS_sig_nls_noMahDis,
                             PCC_all = temp_dt[, stats::cor(A,B) ],                                  # PCC of standard data
                             PCC_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B) ],    # PCC of data without StuRes outliers
                             PCC_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B) ],    # PCC of data without MahDis outliers
                             rho_all = temp_dt[, stats::cor(A,B, method = "spearman") ],                                 # rho of standard data
                             rho_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # rho of data without StuRes outliers
                             rho_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # rho of data without MahDis outliers
                             bic_all = temp_dt[, as.numeric( bicor(A,B)) ],                                  # bic of standard data
                             bic_noStuRes = temp_dt[ StuRes_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without StuRes outliers
                             bic_noMahDis = temp_dt[ MahDis_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without MahDis outliers
                             N_StuRes_outlier = temp_dt[ StuRes_outlier == TRUE , .N ],               # Number of StuRes outliers
                             N_MahDis_outlier = temp_dt[ MahDis_outlier == TRUE , .N ],               # Number of MahDis outliers
                             N = temp_dt[,.N] ))                                                      # Total number of common measurements
  }
  setTxtProgressBar(pb, i)   # Update loop progress
}

trC_vs_PCC <- merge( trC_vs_PCC, res, by = c("Protein_1", "Protein_2") )      # Append the results to the trC_vs_PCC table
trC_vs_PCC <- trC_vs_PCC[ N >= 50 ]                                           # Remove pairs for which no models were fitted

# Clear workspace and memory
rm( list = setdiff( ls(), c("ProHD", "trC_vs_PCC", "DT") ))
gc()


#### treeClust vs RHO ####

# Calculate RHO
ProHD_cor <- stats::cor(ProHD, use = "pairwise.complete.obs", method = c("spearman")) # Get correlation matrix
ProHD_cor <- as.data.table( melt( ProHD_cor))                                         # Convert it to a long data table
ProHD_cor <- ProHD_cor[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), RHO = value ) ]
ProHD_cor <- ProHD_cor[ Protein_1 > Protein_2 ]                                       # Remove self-comparisons and duplicate pairs
ProHD_cor <- ProHD_cor[ complete.cases(ProHD_cor) ]                                     # Drop the (few) pairs which didn't yield a RHO

# Combine with treeClust scores
DT <- merge( DT, ProHD_cor, by = c("Protein_1", "Protein_2"))

# Clear the large ProHD_cor object and clear memory 
rm( "ProHD_cor" )
gc()

# Rank protein pairs by co-regulation score or RHO
DT[, RHO_rank := frank( -RHO ) ]

# Define "high" and "low" thresholds as a function of `rank percentage`
high_scoring <- DT[,.N]/100 * 0.1       # The top-ranking 0.1% of pairs
low_scoring <- DT[,.N]/100 * 0.5       # Not in the top 0.5% of pairs 

# Define protein pairs that were either detected by treeClust or by RHO
trC_only <- DT[ coreg_score_rank < high_scoring & RHO_rank > low_scoring  ]
RHO_only <- DT[ coreg_score_rank > low_scoring  & RHO_rank < high_scoring ]

# Combine into one test set
trC_only[, group := "trC_only" ]
RHO_only[, group := "RHO_only" ]
trC_vs_RHO <- rbind( trC_only, RHO_only )


## Model fitting & outlier removal ##

pb <- txtProgressBar(min = 0, max = trC_vs_RHO[,.N], style = 3)  # Initiate progress bar 
res <- data.table()                                              # Initiate result table

for(i in 1:trC_vs_RHO[,.N] ){
  my_pair <- trC_vs_RHO[i, c(Protein_1, Protein_2)]            # Define protein pair to assess in current iteration
  temp_dt <- na.omit( as.data.table( ProHD[,  my_pair ] ))       # Define those protein pair's data        
  names( temp_dt ) <- c("A", "B")                              # Replace protein IDs as colnames
  temp_dt <- temp_dt[, lapply(.SD, rescale )]                  # Rescale for curve-fitting
  
  if( temp_dt[,.N] < 50 ){                                     # If there are fewer than 50 common data points, skip the pair                          
    res <- rbind(res,                                          # ... and rather than fitting a model just add NAs to the results
                 data.table( Protein_1 = my_pair[1],           # See below for the definition of these parameters
                             Protein_2 = my_pair[2],
                             lm_all_r2 = NA,
                             lm_noStuRes_r2 = NA,
                             lm_noMahDis_r2 = NA,
                             mae_lm_mod = NA,
                             RSS_lm_mod = NA,
                             RSS_lin_nls = NA,
                             RSS_lin_nls_noStuRes = NA,
                             RSS_lin_nls_noMahDis = NA,
                             RSS_exp_nls = NA,
                             RSS_exp_nls_noStuRes = NA, 
                             RSS_exp_nls_noMahDis = NA,
                             RSS_sig_nls = NA,
                             RSS_sig_nls_noStuRes = NA,
                             RSS_sig_nls_noMahDis = NA,
                             PCC_all = NA,
                             PCC_noStuRes = NA,
                             PCC_noMahDis = NA,
                             rho_all = NA,
                             rho_noStuRes = NA,
                             rho_noMahDis = NA,
                             bic_all = NA,
                             bic_noStuRes = NA,
                             bic_noMahDis = NA,
                             N_StuRes_outlier = NA,
                             N_MahDis_outlier = NA,
                             N = temp_dt[,.N] ))
    
  } else {                                                    # If there are more than 50 common data, proceed as follows
    
    ## Fit linear regression model, identify outliers (via studentized residuals or mahalanobis distance), re-fit model without outliers
    
    lm_mod <- lm(A ~ B, data = temp_dt)                       # Fit a simple linear regression model
    temp_dt$StuRes <- rstudent(lm_mod)                        # Extract the studentized residuals of the model
    temp_dt$MahDis <- mahalanobis( temp_dt[,.(A,B)],          # Calculate the Mahalanobis distance for each point
                                   center = colMeans( temp_dt[,.(A,B)] ),
                                   cov( temp_dt[,.(A,B)] ))
    
    temp_dt[, StuRes_outlier := abs( StuRes ) > 2 ]           # Mark data points that are outliers according to studentized residuals
    temp_dt[, MahDis_outlier := MahDis > 2 ]                  # Mark data points that are outliers according to the Mahalanobis distance
    
    lm_noStuRes <- lm(A ~ B, data = temp_dt[ StuRes_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    lm_noMahDis <- lm(A ~ B, data = temp_dt[ MahDis_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    
    
    ## Fit non-linear models (exponential and sigmoid, as well as linear for comparison), output NaN if fitting does not converge
    
    # For all data
    lin_nls <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt, start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without StuRes outliers
    lin_nls_noStuRes <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noStuRes <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noStuRes <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without MahDis outliers
    lin_nls_noMahDis <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noMahDis <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noMahDis <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # Extract residual sum of squares (RSS) for pairs where model was built successfully
    RSS_lin_nls          <- if( class(lin_nls) == "nls" ) {          lin_nls$m$deviance()          } else { NaN }
    RSS_lin_nls_noStuRes <- if( class(lin_nls_noStuRes) == "nls" ) { lin_nls_noStuRes$m$deviance() } else { NaN }
    RSS_lin_nls_noMahDis <- if( class(lin_nls_noMahDis) == "nls" ) { lin_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_exp_nls          <- if( class(exp_nls) == "nls" ) {          exp_nls$m$deviance()          } else { NaN }
    RSS_exp_nls_noStuRes <- if( class(exp_nls_noStuRes) == "nls" ) { exp_nls_noStuRes$m$deviance() } else { NaN }
    RSS_exp_nls_noMahDis <- if( class(exp_nls_noMahDis) == "nls" ) { exp_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_sig_nls          <- if( class(sig_nls) == "nls" ) {          sig_nls$m$deviance()          } else { NaN }
    RSS_sig_nls_noStuRes <- if( class(sig_nls_noStuRes) == "nls" ) { sig_nls_noStuRes$m$deviance() } else { NaN }
    RSS_sig_nls_noMahDis <- if( class(sig_nls_noMahDis) == "nls" ) { sig_nls_noMahDis$m$deviance() } else { NaN }
    
    
    ## Collect results in a data.table
    res <- rbind(res,                                                            
                 data.table( Protein_1 = my_pair[1],                               # Name of protein 1
                             Protein_2 = my_pair[2],                               # Name of protein 2
                             lm_all_r2 =      summary( lm_mod      )$r.squared,    # R2 of standard linear model
                             lm_noStuRes_r2 = summary( lm_noStuRes )$r.squared,    # R2 of linear model without StuRes outliers
                             lm_noMahDis_r2 = summary( lm_noMahDis )$r.squared,    # R2 of linear model without MahDis outliers
                             mae_lm_mod = mean(abs(lm_mod$residuals)),             # Mean absolute error of the standard linear model
                             RSS_lm_mod = sum( lm_mod$residuals^2 ),               # Residual sum of squares of standard linear model
                             RSS_lin_nls,                                          # Residual sum of squares from the various non-linear fittings
                             RSS_lin_nls_noStuRes,
                             RSS_lin_nls_noMahDis,
                             RSS_exp_nls,
                             RSS_exp_nls_noStuRes, 
                             RSS_exp_nls_noMahDis,
                             RSS_sig_nls,
                             RSS_sig_nls_noStuRes,
                             RSS_sig_nls_noMahDis,
                             PCC_all = temp_dt[, stats::cor(A,B) ],                                  # RHO of standard data
                             PCC_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B) ],    # RHO of data without StuRes outliers
                             PCC_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B) ],    # RHO of data without MahDis outliers
                             rho_all = temp_dt[, stats::cor(A,B, method = "spearman") ],                                 # rho of standard data
                             rho_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # rho of data without StuRes outliers
                             rho_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # rho of data without MahDis outliers
                             bic_all = temp_dt[, as.numeric( bicor(A,B)) ],                                  # bic of standard data
                             bic_noStuRes = temp_dt[ StuRes_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without StuRes outliers
                             bic_noMahDis = temp_dt[ MahDis_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without MahDis outliers
                             N_StuRes_outlier = temp_dt[ StuRes_outlier == TRUE , .N ],               # Number of StuRes outliers
                             N_MahDis_outlier = temp_dt[ MahDis_outlier == TRUE , .N ],               # Number of MahDis outliers
                             N = temp_dt[,.N] ))                                                      # Total number of common measurements
  }
  setTxtProgressBar(pb, i)   # Update loop progress
}

trC_vs_RHO <- merge( trC_vs_RHO, res, by = c("Protein_1", "Protein_2") )      # Append the results to the trC_vs_RHO table
trC_vs_RHO <- trC_vs_RHO[ N >= 50 ]                                           # Remove pairs for which no models were fitted

# Clear workspace and memory
rm( list = setdiff( ls(), c("ProHD", "trC_vs_PCC", "trC_vs_RHO", "DT") ))
gc()


#### treeClust vs BIC ####

# Calculate BIC
ProHD_cor <- bicor(ProHD, use = "pairwise.complete.obs" )                             # Get correlation matrix
ProHD_cor <- as.data.table( melt( ProHD_cor))                                         # Convert it to a long data table
ProHD_cor <- ProHD_cor[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), BIC = value ) ]
ProHD_cor <- ProHD_cor[ Protein_1 > Protein_2 ]                                       # Remove self-comparisons and duplicate pairs
ProHD_cor <- ProHD_cor[ complete.cases(ProHD_cor) ]                                     # Drop the (few) pairs which didn't yield a BIC

# Combine with treeClust scores
DT <- merge( DT, ProHD_cor, by = c("Protein_1", "Protein_2"))

# Clear the large ProHD_cor object and clear memory 
rm( "ProHD_cor" )
gc()

# Rank protein pairs by co-regulation score or BIC
DT[, BIC_rank := frank( -BIC ) ]

# Define "high" and "low" thresholds as a function of `rank percentage`
high_scoring <- DT[,.N]/100 * 0.1       # The top-ranking 0.1% of pairs
low_scoring <- DT[,.N]/100 * 0.5       # Not in the top 0.5% of pairs 

# Define protein pairs that were either detected by treeClust or by BIC
trC_only <- DT[ coreg_score_rank < high_scoring & BIC_rank > low_scoring  ]
BIC_only <- DT[ coreg_score_rank > low_scoring  & BIC_rank < high_scoring ]

# Combine into one test set
trC_only[, group := "trC_only" ]
BIC_only[, group := "BIC_only" ]
trC_vs_BIC <- rbind( trC_only, BIC_only )


## Model fitting & outlier removal ##

pb <- txtProgressBar(min = 0, max = trC_vs_BIC[,.N], style = 3)  # Initiate progress bar 
res <- data.table()                                              # Initiate result table

for(i in 1:trC_vs_BIC[,.N] ){
  my_pair <- trC_vs_BIC[i, c(Protein_1, Protein_2)]            # Define protein pair to assess in current iteration
  temp_dt <- na.omit( as.data.table( ProHD[,  my_pair ] ))       # Define those protein pair's data        
  names( temp_dt ) <- c("A", "B")                              # Replace protein IDs as colnames
  temp_dt <- temp_dt[, lapply(.SD, rescale )]                  # Rescale for curve-fitting
  
  if( temp_dt[,.N] < 50 ){                                     # If there are fewer than 50 common data points, skip the pair                          
    res <- rbind(res,                                          # ... and rather than fitting a model just add NAs to the results
                 data.table( Protein_1 = my_pair[1],           # See below for the definition of these parameters
                             Protein_2 = my_pair[2],
                             lm_all_r2 = NA,
                             lm_noStuRes_r2 = NA,
                             lm_noMahDis_r2 = NA,
                             mae_lm_mod = NA,
                             RSS_lm_mod = NA,
                             RSS_lin_nls = NA,
                             RSS_lin_nls_noStuRes = NA,
                             RSS_lin_nls_noMahDis = NA,
                             RSS_exp_nls = NA,
                             RSS_exp_nls_noStuRes = NA, 
                             RSS_exp_nls_noMahDis = NA,
                             RSS_sig_nls = NA,
                             RSS_sig_nls_noStuRes = NA,
                             RSS_sig_nls_noMahDis = NA,
                             PCC_all = NA,
                             PCC_noStuRes = NA,
                             PCC_noMahDis = NA,
                             rho_all = NA,
                             rho_noStuRes = NA,
                             rho_noMahDis = NA,
                             bic_all = NA,
                             bic_noStuRes = NA,
                             bic_noMahDis = NA,
                             N_StuRes_outlier = NA,
                             N_MahDis_outlier = NA,
                             N = temp_dt[,.N] ))
    
  } else {                                                    # If there are more than 50 common data, proceed as follows
    
    ## Fit linear regression model, identify outliers (via studentized residuals or mahalanobis distance), re-fit model without outliers
    
    lm_mod <- lm(A ~ B, data = temp_dt)                       # Fit a simple linear regression model
    temp_dt$StuRes <- rstudent(lm_mod)                        # Extract the studentized residuals of the model
    temp_dt$MahDis <- mahalanobis( temp_dt[,.(A,B)],          # Calculate the Mahalanobis distance for each point
                                   center = colMeans( temp_dt[,.(A,B)] ),
                                   cov( temp_dt[,.(A,B)] ))
    
    temp_dt[, StuRes_outlier := abs( StuRes ) > 2 ]           # Mark data points that are outliers according to studentized residuals
    temp_dt[, MahDis_outlier := MahDis > 2 ]                  # Mark data points that are outliers according to the Mahalanobis distance
    
    lm_noStuRes <- lm(A ~ B, data = temp_dt[ StuRes_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    lm_noMahDis <- lm(A ~ B, data = temp_dt[ MahDis_outlier == FALSE ])    # Re-fit model without studentized residual outliers
    
    
    ## Fit non-linear models (exponential and sigmoid, as well as linear for comparison), output NaN if fitting does not converge
    
    # For all data
    lin_nls <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt, start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt, start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without StuRes outliers
    lin_nls_noStuRes <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noStuRes <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noStuRes <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ StuRes_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # For data without MahDis outliers
    lin_nls_noMahDis <- tryCatch( nls( A ~ (a + b*B                    ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    exp_nls_noMahDis <- tryCatch( nls( A ~ (a + exp(b)^B               ), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 0, b = 1),          control = nls.control( maxiter = 100 )), error = function(e) NaN )
    sig_nls_noMahDis <- tryCatch( nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = temp_dt[ MahDis_outlier == FALSE ], start = list(a = 1, b = 5, c = 0.5), control = nls.control( maxiter = 100 )), error = function(e) NaN )
    
    # Extract residual sum of squares (RSS) for pairs where model was built successfully
    RSS_lin_nls          <- if( class(lin_nls) == "nls" ) {          lin_nls$m$deviance()          } else { NaN }
    RSS_lin_nls_noStuRes <- if( class(lin_nls_noStuRes) == "nls" ) { lin_nls_noStuRes$m$deviance() } else { NaN }
    RSS_lin_nls_noMahDis <- if( class(lin_nls_noMahDis) == "nls" ) { lin_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_exp_nls          <- if( class(exp_nls) == "nls" ) {          exp_nls$m$deviance()          } else { NaN }
    RSS_exp_nls_noStuRes <- if( class(exp_nls_noStuRes) == "nls" ) { exp_nls_noStuRes$m$deviance() } else { NaN }
    RSS_exp_nls_noMahDis <- if( class(exp_nls_noMahDis) == "nls" ) { exp_nls_noMahDis$m$deviance() } else { NaN }
    
    RSS_sig_nls          <- if( class(sig_nls) == "nls" ) {          sig_nls$m$deviance()          } else { NaN }
    RSS_sig_nls_noStuRes <- if( class(sig_nls_noStuRes) == "nls" ) { sig_nls_noStuRes$m$deviance() } else { NaN }
    RSS_sig_nls_noMahDis <- if( class(sig_nls_noMahDis) == "nls" ) { sig_nls_noMahDis$m$deviance() } else { NaN }
    
    
    ## Collect results in a data.table
    res <- rbind(res,                                                            
                 data.table( Protein_1 = my_pair[1],                               # Name of protein 1
                             Protein_2 = my_pair[2],                               # Name of protein 2
                             lm_all_r2 =      summary( lm_mod      )$r.squared,    # R2 of standard linear model
                             lm_noStuRes_r2 = summary( lm_noStuRes )$r.squared,    # R2 of linear model without StuRes outliers
                             lm_noMahDis_r2 = summary( lm_noMahDis )$r.squared,    # R2 of linear model without MahDis outliers
                             mae_lm_mod = mean(abs(lm_mod$residuals)),             # Mean absolute error of the standard linear model
                             RSS_lm_mod = sum( lm_mod$residuals^2 ),               # Residual sum of squares of standard linear model
                             RSS_lin_nls,                                          # Residual sum of squares from the various non-linear fittings
                             RSS_lin_nls_noStuRes,
                             RSS_lin_nls_noMahDis,
                             RSS_exp_nls,
                             RSS_exp_nls_noStuRes, 
                             RSS_exp_nls_noMahDis,
                             RSS_sig_nls,
                             RSS_sig_nls_noStuRes,
                             RSS_sig_nls_noMahDis,
                             PCC_all = temp_dt[, stats::cor(A,B) ],                                  # BIC of standard data
                             PCC_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B) ],    # BIC of data without StuRes outliers
                             PCC_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B) ],    # BIC of data without MahDis outliers
                             rho_all = temp_dt[, stats::cor(A,B, method = "spearman") ],                                 # BIC of standard data
                             rho_noStuRes = temp_dt[ StuRes_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # BIC of data without StuRes outliers
                             rho_noMahDis = temp_dt[ MahDis_outlier == FALSE , stats::cor(A,B, method = "spearman") ],   # BIC of data without MahDis outliers
                             bic_all = temp_dt[, as.numeric( bicor(A,B)) ],                                  # bic of standard data
                             bic_noStuRes = temp_dt[ StuRes_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without StuRes outliers
                             bic_noMahDis = temp_dt[ MahDis_outlier == FALSE , as.numeric( bicor(A,B)) ],    # bic of data without MahDis outliers
                             N_StuRes_outlier = temp_dt[ StuRes_outlier == TRUE , .N ],               # Number of StuRes outliers
                             N_MahDis_outlier = temp_dt[ MahDis_outlier == TRUE , .N ],               # Number of MahDis outliers
                             N = temp_dt[,.N] ))                                                      # Total number of common measurements
  }
  setTxtProgressBar(pb, i)   # Update loop progress
}

trC_vs_BIC <- merge( trC_vs_BIC, res, by = c("Protein_1", "Protein_2") )      # Append the results to the trC_vs_BIC table
trC_vs_BIC <- trC_vs_BIC[ N >= 50 ]                                           # Remove pairs for which no models were fitted

# Clear workspace and memory
rm( list = setdiff( ls(), c("ProHD", "trC_vs_PCC", "trC_vs_RHO", "trC_vs_BIC", "DT") ))
gc()


#### Prep result data for plotting ####

# Plot formatting settings
my_plot_theme <-  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25),
                        legend.position = "none")
my_cols <- c( bicor = "mediumblue", PCC = "steelblue2", RHO = "springgreen", treeClust = "magenta2" )


## Merge the revelant data into one table for easy plotting

# Assign identifier for after merging
trC_vs_PCC[, type := "trC_vs_PCC" ]    
trC_vs_RHO[, type := "trC_vs_RHO" ]
trC_vs_BIC[, type := "trC_vs_BIC" ]

# Define subset of columns relevant for plotting
relevant_cols <- c("type", "group", "Class",
                   "mae_lm_mod", "RSS_lm_mod",
                   "RSS_lin_nls", "RSS_lin_nls_noMahDis", "RSS_exp_nls", "RSS_exp_nls_noMahDis", "RSS_sig_nls", "RSS_sig_nls_noMahDis",
                   "PCC_all", "PCC_noStuRes", "PCC_noMahDis",
                   "rho_all", "rho_noStuRes", "rho_noMahDis",
                   "bic_all", "bic_noStuRes", "bic_noMahDis",
                   "N_StuRes_outlier", "N_MahDis_outlier", "N")

# Merge
plot_dt <- rbind( trC_vs_PCC[, relevant_cols, with = FALSE ],
                  trC_vs_RHO[, relevant_cols, with = FALSE ],
                  trC_vs_BIC[, relevant_cols, with = FALSE ])

# Assign plotting order for "type"
plot_dt[, type := factor(type, levels = c("trC_vs_PCC", "trC_vs_RHO", "trC_vs_BIC"))]

# Simplify "group" and assign plotting order
plot_dt[, group := gsub("PCC", "cor", group) ]
plot_dt[, group := gsub("RHO", "cor", group) ]
plot_dt[, group := gsub("BIC", "cor", group) ]
plot_dt[, group := factor(group, levels = c("trC_only", "cor_only")) ]


#### Plot outlier example scatterplots #### 

## Example scatterplots to illustrate regression vs Mahalanobis outliers

# Select a false-positive example pair and find outliers as above
example_pair <- c("Q14344;Q14344-2", "P61006;P61006-2")
example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))
names( example_pair ) <- c("A", "B")
example_pair <- example_pair[, lapply(.SD, rescale )] # Rescale for curve-fitting
example_pair_lm <- lm(A ~ B, data = example_pair)     # Fit a simple linear regression model
StuRes <- rstudent(example_pair_lm)                   # Extract the studentized residuals of the model
MahDis <- mahalanobis( example_pair[,.(A,B)] ,        # Calculate the Mahalanobis distance for each point
                       center = colMeans( example_pair[,.(A,B)] ),
                       cov( example_pair[,.(A,B)] ))
example_pair[, StuRes := abs( StuRes ) > 2 ]          # Mark data points that are outliers according to studentized residuals
example_pair[, MahDis := MahDis > 2 ]                 # Mark data points that are outliers according to the Mahalanobis distance

pFP_scatter <- ggplot( example_pair, aes(x = A, y = B))+
        geom_point( size = 0.1 , colour = "grey50" )+
        geom_smooth( method = "lm", fullrange = TRUE, size = 0.25, se = FALSE, colour = "black", linetype = "dashed")+
        geom_point( data = example_pair[ MahDis == TRUE ] , colour = "deeppink2", size = 0.5 )+
        xlab("G-protein subunit alpha-13")+
        ylab("Ras-related protein Rab-8A")+
        annotate("text", label = "false positive", x = 0.8, y = 0.2, size = 2.1 )+
        xlim(0,1)+
        ylim(0,1)+
        my_plot_theme

# Select a true-positive example pair and find outliers as above
example_pair <- c("P68431", "P62805")
example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))
names( example_pair ) <- c("A", "B")
example_pair <- example_pair[, lapply(.SD, rescale )] # Rescale for curve-fitting
example_pair_lm <- lm(A ~ B, data = example_pair)     # Fit a simple linear regression model
StuRes <- rstudent(example_pair_lm)                   # Extract the studentized residuals of the model
MahDis <- mahalanobis( example_pair[,.(A,B)] ,        # Calculate the Mahalanobis distance for each point
                       center = colMeans( example_pair[,.(A,B)] ),
                       cov( example_pair[,.(A,B)] ))
example_pair[, StuRes := abs( StuRes ) > 2 ]          # Mark data points that are outliers according to studentized residuals
example_pair[, MahDis := MahDis > 2 ]                 # Mark data points that are outliers according to the Mahalanobis distance

pTP_scatter <- ggplot( example_pair, aes(x = A, y = B))+
        geom_point( size = 0.1 , colour = "grey50" )+
        geom_smooth( method = "lm", fullrange = TRUE, size = 0.25, se = FALSE, colour = "black")+
        geom_point( data = example_pair[ StuRes == TRUE ] , size = 0.5, shape = 21, fill = NA, colour = "deeppink2")+
        xlab("Histone H3")+
        ylab("Histone H4")+
        annotate("text", label = "true positive", x = 0.8, y = 0.2, size = 2.1 )+      
        xlim(0,1)+
        ylim(0,1)+
        my_plot_theme


#### Plot true and false positive stats (pie charts) ####

# Calculate TP percentages
pcTP <- plot_dt[, round( sum( Class == "TP", na.rm = TRUE ) / sum( !is.na(Class) ) *100, 2), .(type, group) 
                ][, .(type, group, pcTP_label = paste(V1, "% TP", sep = ""), Class = "TP") ]

# Get total protein numbers
totN <-  plot_dt[, .N, .(type, group) ][, .(type, group, totN_label = paste("n = ", N, sep = ""), Class = NA) ]

# Create the pie charts
pPie <- ggplot( plot_dt[ !is.na(Class) , .N, .(type, group, Class)], aes(x = factor(1), y = N , fill = Class))+
          geom_bar( stat = "identity", position = "fill") + coord_polar(theta = "y")+
          facet_grid( ~ type + group )+
          scale_fill_manual( values = c( TP = "dodgerblue", FP = "salmon" ))+
          geom_text( data = pcTP, aes(label = pcTP_label, x = 1, y = 1), size = 2)+
          geom_text( data = totN, aes(label = totN_label, x = 1, y = 0.5), size = 2)+
          my_plot_theme + theme( strip.background = element_rect(fill = NA, colour = "black", size = 0.25),
                                 strip.text = element_text(size = 6),
                                 axis.text = element_blank(), axis.ticks = element_blank(),
                                 axis.title = element_blank(), legend.position = "top")


#### Plot outlier boxplots #### 

# Regression outlier numbers
pOut_1 <- ggplot(plot_dt, aes( x = group, y = (N_StuRes_outlier/N*100) ))+
  facet_grid(.~type)+
  geom_boxplot( notch = TRUE , fill = "deeppink2", outlier.color = NA, lwd = 0.25)+
  scale_y_continuous( breaks = seq(0,16,2))+
  coord_cartesian( ylim = c(0,12))+
  ylab("Regression outliers [%]")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major.y = element_line(colour = "grey50", linetype = "dashed", size = 0.25),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Mahalanobis outlier numbers
pOut_2 <- ggplot(plot_dt, aes( x = group, y = (N_MahDis_outlier/N*100) ))+
  facet_grid(.~type)+
  geom_boxplot( notch = TRUE , fill = "lightseagreen", outlier.color = NA, lwd = 0.25)+
  scale_y_continuous( breaks = seq(0,50,10))+
  coord_cartesian( ylim = c(0,50))+
  ylab("Mahalanobis outliers [%]")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major.y = element_line(colour = "grey50", linetype = "dashed", size = 0.25),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())


## Outlier removal: Can't do this in one plot because I'm plotting either PCC, RHO or BIC
## So I'll remove the appropriate panels in Inkscape later

# Regression outlier removal (PCC)
pOut_3a <- ggplot(plot_dt, aes( x = group, y = (PCC_noStuRes - PCC_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-0.12, 0.35) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "deeppink2", lwd = 0.25)+
  ylab("PCC w/o outliers - PCC with outliers")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Regression outlier removal (RHO)
pOut_3b <- ggplot(plot_dt, aes( x = group, y = (rho_noStuRes - rho_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-0.12, 0.35) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "deeppink2", lwd = 0.25)+
  ylab("RHO w/o outliers - RHO with outliers")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Regression outlier removal (BIC)
pOut_3c <- ggplot(plot_dt, aes( x = group, y = (bic_noStuRes - bic_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-0.12, 0.35) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "deeppink2", lwd = 0.25)+
  ylab("BIC w/o outliers - BIC with outliers")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())


# Mahalanobis outlier removal (PCC)
pOut_4a <- ggplot(plot_dt, aes( x = group, y = (PCC_noMahDis - PCC_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-1, 0.5) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "lightseagreen", lwd = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  ylab("PCC w/o outliers - PCC with outliers")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Mahalanobis outlier removal (RHO)
pOut_4b <- ggplot(plot_dt, aes( x = group, y = (rho_noMahDis - rho_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-1, 0.5) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "lightseagreen", lwd = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  ylab("RHO w/o outliers - RHO with outliers")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Mahalanobis outlier removal (RHO)
pOut_4c <- ggplot(plot_dt, aes( x = group, y = (bic_noMahDis - bic_all)))+
  facet_grid(.~type)+
  coord_cartesian( ylim = c(-1, 0.5) )+
  geom_boxplot(notch = TRUE , outlier.colour = NA , fill = "lightseagreen", lwd = 0.25)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+
  ylab("BIC w/o outliers - BIC with outliers")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())


#### Plot dispersion example scatterplots #### 

## Example scatterplots to illustrate different MAEs

# Rescale the similarity / correlation values to make them comparable
DT[, tC_sim_scaled := rescale( tC_sim ) ]
DT[, PCC_scaled := rescale( PCC ) ]
DT[, RHO_scaled := rescale( RHO ) ]
DT[, BIC_scaled := rescale( BIC ) ]

# Select a true-positive example pair with low MAE
example_pair <- c("P62847;P62847-2;P62847-3", "P61254")
TP_pair_tC_score <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , tC_sim]  # Find the treeClust similarity of this pair
TP_pair_BIC <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , BIC]          # Find the correlation values of this pair
TP_pair_PCC <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , PCC]
TP_pair_RHO <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , RHO]

 TP_pair_tC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , tC_sim_scaled]  # Find the corresponding rescaled values
TP_pair_BIC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , BIC_scaled]         
TP_pair_PCC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , PCC_scaled]
TP_pair_RHO_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , RHO_scaled]

example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))       # Get SILAC ratios
names( example_pair ) <- c("A", "B")                                    # Rename columns for model fitting
example_pair <- example_pair[, lapply(.SD, rescale )]                   # Rescale for curve-fitting
example_pair_lm <- lm(A ~ B, data = example_pair)                       # Fit a simple linear regression model
example_pair$residuals <- example_pair_lm$residuals                     # Extract the residuals
example_pair_mae <- mean(abs(example_pair$residuals))                   # Calculate MAE

pD_TP <- ggplot( example_pair, aes(x = B, y = A))+
          geom_smooth( method = "lm", fullrange = TRUE, size = 0.25, se = FALSE, colour = "black")+
          geom_point( size = 0.1 )+
          geom_segment( aes(x = B, y = A, xend = B, yend = A-residuals ), size = 0.25, colour = "orange")+
          xlab("Ribosomal protein L26")+
          ylab("Ribosomal protein S24")+
          annotate("text", label = "true positive", x = 0.8, y = 0.2, size = 2.1 )+
          coord_cartesian( xlim = c(0,1), ylim = c(0,1))+
          annotate("text", x = 0.2, y = 0.9, label = paste( "MAE", round( example_pair_mae, 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.7, label = paste( "BIC", round( TP_pair_BIC     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.6, label = paste( "PCC", round( TP_pair_PCC     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.5, label = paste( "RHO", round( TP_pair_RHO     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.4, label = paste( "treeClust", round(TP_pair_tC_score, 2)), size = 2.1)+
          my_plot_theme


# Select a false-positive example pair with high MAE  (see previous example for code annotation)
example_pair <- c("Q08380", "P02656")

FP_pair_tC_score <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , tC_sim]  # Find the treeClust similarity of this pair
FP_pair_BIC <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , BIC]          # Find the correlation values of this pair
FP_pair_PCC <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , PCC]
FP_pair_RHO <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , RHO]

 FP_pair_tC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , tC_sim_scaled]  # Find the corresponding rescaled values
FP_pair_BIC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , BIC_scaled]         
FP_pair_PCC_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , PCC_scaled]
FP_pair_RHO_scaled <- DT[ Protein_1 == example_pair[1] & Protein_2 == example_pair[2] , RHO_scaled]

example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))
names( example_pair ) <- c("A", "B")
example_pair <- example_pair[, lapply(.SD, rescale )] 
example_pair_lm <- lm(A ~ B, data = example_pair) 
example_pair$residuals <- example_pair_lm$residuals 
example_pair_mae <- mean(abs(example_pair$residuals))

pD_FP <- ggplot( example_pair, aes(x = B, y = A))+
          geom_smooth( method = "lm", fullrange = TRUE, size = 0.25, se = FALSE, colour = "black")+
          geom_point( size = 0.1 )+
          geom_segment( aes(x = B, y = A, xend = B, yend = A-residuals ), size = 0.25, colour = "orange")+
          xlab("Galectin-3-binding protein")+
          ylab("Apolipoprotein C-III")+
          annotate("text", label = "false positive", x = 0.8, y = 0.2, size = 2.1 )+
          coord_cartesian( xlim = c(0,1), ylim = c(0,1))+
          annotate("text", x = 0.2, y = 0.9, label = paste( "MAE", round( example_pair_mae, 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.7, label = paste( "BIC", round( FP_pair_BIC     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.6, label = paste( "PCC", round( FP_pair_PCC     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.5, label = paste( "RHO", round( FP_pair_RHO     , 2)), size = 2.1)+
          annotate("text", x = 0.2, y = 0.4, label = paste( "treeClust", round(FP_pair_tC_score, 2)), size = 2.1)+
          my_plot_theme


## Example barplot to show score differences

# Plot the rescaled values
TP_scaled <- c(treeClust = TP_pair_tC_scaled, PCC = TP_pair_PCC_scaled, RHO = TP_pair_RHO_scaled, BIC = TP_pair_BIC_scaled)
FP_scaled <- c(treeClust = FP_pair_tC_scaled, PCC = FP_pair_PCC_scaled, RHO = FP_pair_RHO_scaled, BIC = FP_pair_BIC_scaled)
both_scaled <- data.table( TP_scaled , FP_scaled, metric = names(TP_scaled) )
both_scaled <- melt(both_scaled, id.vars = "metric")
both_scaled[, metric := factor(metric, levels = c("PCC", "RHO", "BIC", "treeClust")) ]                 # Set plotting order

pEx_bar <- ggplot(both_scaled, aes( x = variable, y = value, colour = variable, group = metric))+
            facet_wrap(~metric, nrow = 1)+
            geom_line(colour = "black", size = 0.25)+
            geom_point(size = 2)+
            ylim(0,1)+
            ylab("Correlation coefficient / similarity\n(scaled to [0,1] for comparison)")+
            scale_colour_manual( values = c("red", "green"))+
            theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
                  axis.text=element_text(size=5), axis.title.y=element_text(size=6), axis.ticks.y = element_line(size=0.25),
                  axis.ticks.x = element_blank(), axis.title.x = element_blank())

pEx_bar_v2 <- ggplot(both_scaled, aes( x = metric, y = value, fill = variable))+
                geom_bar(stat = "identity", position = "dodge", colour = "black", size = 0.25)+
                ylab("Correlation coefficient / similarity\n(scaled to [0,1] for comparison)")+
                scale_fill_manual( values = c("white", "black"))+
                theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
                      axis.text=element_text(size=5), axis.title.y=element_text(size=6), axis.ticks.y = element_line(size=0.25),
                      axis.ticks.x = element_blank(), axis.title.x = element_blank())

#### Plot impact of data dispersion #### 

pMAEbx <- ggplot( plot_dt, aes( x = group, y = mae_lm_mod ))+
  geom_boxplot( notch = TRUE , outlier.colour = "black" )+
  coord_cartesian( ylim = c(0, 0.2) )+
  facet_grid( ~ type, scales = "free" )+
  ylab("Mean absolute error (MAE")+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major.y = element_line(colour = "grey50", linetype = "dashed", size = 0.25),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())


pMAE <- ggplot( plot_dt, aes( fill = group, x = mae_lm_mod ))+
  geom_histogram( position = "dodge", binwidth = 0.01, boundary = 0)+
  scale_x_continuous( limits = c(0, 0.2), expand = c(0,0))+
  facet_grid( ~ type )+
  xlab("Mean absolute error (MAE)")+
  ylab("Number of protein pairs")+
  scale_fill_manual( values = c("darkslateblue", "magenta", "darkturquoise"))+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks = element_line(size=0.25))


#### Non-linear relationships: Example scatterplots ####

# Prepare an example where the exponential curve fits better than the linear model
example_pair <- c("P08567", "O75558")
example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))
names( example_pair ) <- c("A", "B")
example_pair <- example_pair[, lapply(.SD, rescale )]      # Rescale for curve-fitting

lin_nls <- nls( A ~ (a + b*B                    ), data = example_pair, start = list(a = 0, b = 1)         )
exp_nls <- nls( A ~ (a + exp(b)^B               ), data = example_pair, start = list(a = 0, b = 1)         )

lin_RSS <- lin_nls$m$deviance()
exp_RSS <- exp_nls$m$deviance()

lin_nls_1 <- lin_nls$m$getPars()  # Extract the fitted parameters (_1 is necessary to distinguish the fit from the next example)
exp_nls <- exp_nls$m$getPars()    # Extract the fitted parameters

pExp <- ggplot( example_pair, aes(x = B, y = A))+
          geom_point( size = 0.1 )+
          stat_function( fun = function(x){ lin_nls_1["a"] + lin_nls_1["b"]*x },    geom = "line", colour = "black", size = 0.25)+
          stat_function( fun = function(x){ exp_nls["a"] + exp( exp_nls["b"] )^x }, geom = "line", colour = "mediumblue", size = 0.25)+
          annotate("text", x = 0.1, y = 0.9, colour = "black",    size = 2.1, label = paste("Linear model: RSS", round( lin_RSS, 2)))+
          annotate("text", x = 0.1, y = 0.8, colour = "mediumblue", size = 2.1, label = paste("Exponential model: RSS", round( exp_RSS, 2)))+ 
          xlab("Pleckstrin")+
          ylab("Syntaxin-11")+
          scale_x_continuous( limits = c(0,1), breaks = c(0.0, 0.5, 1.0))+
          scale_y_continuous( limits = c(0,1), breaks = c(0.0, 0.5, 1.0))+
          my_plot_theme


# Prepare an example where the logistic curve fits better than the linear model
example_pair <- c("P61224;P61224-3;A6NIZ1;P61224-2", "P04040")
example_pair <- na.omit( as.data.table( ProHD[, example_pair ] ))
names( example_pair ) <- c("A", "B")
example_pair <- example_pair[, lapply(.SD, rescale )]

lin_nls <- nls( A ~ (a + b*B                    ), data = example_pair, start = list(a = 0, b = 1)         )
sig_nls <- nls( A ~ (a / (1 + exp(1)^(-b*(B-c)))), data = example_pair, start = list(a = 1, b = 5, c = 0.5))

lin_RSS <- lin_nls$m$deviance()
sig_RSS <- sig_nls$m$deviance()

lin_nls_2 <- lin_nls$m$getPars()  # Extract the fitted parameters
sig_nls <- sig_nls$m$getPars()    # Extract the fitted parameters


pSig <- ggplot( example_pair, aes(x = B, y = A))+
          geom_point( size = 0.1 )+
          stat_function( fun = function(x){ lin_nls_2["a"] + lin_nls_2["b"]*x },        geom = "line", colour = "black", size = 0.25)+
          stat_function( fun = function(x){ sig_nls["a"] / (1 + exp(1)^(-sig_nls["b"]*(x-sig_nls["c"]))) }, geom = "line", colour = "red", size = 0.25)+
          annotate("text", x = 0.1, y = 0.9,   colour = "black",    size = 2.1, label = paste("Linear model: RSS", round( lin_RSS, 2)))+
          annotate("text", x = 0.1, y = 0.8, colour = "red", size = 2.1, label = paste("Logistic model: RSS", round( sig_RSS, 2)))+ 
          xlab("Ras-related protein Rap-1b")+
          ylab("Catalase")+
          scale_x_continuous( limits = c(0,1), breaks = c(0.0, 0.5, 1.0))+
          scale_y_continuous( limits = c(0,1), breaks = c(0.0, 0.5, 1.0))+
          my_plot_theme


#### Non-linear relationships: Barplot ####

## Deal with missing values in fitted models
# RSS = NaN indicates that no curve could be fitted, therefore assign those pairs the maximum RSS
maxRSS <- plot_dt[, max(RSS_lin_nls, na.rm = TRUE) ]

plot_dt[ is.na( RSS_lin_nls )          , RSS_lin_nls          := maxRSS ]
plot_dt[ is.na( RSS_lin_nls_noMahDis ) , RSS_lin_nls_noMahDis := maxRSS ]
plot_dt[ is.na( RSS_exp_nls )          , RSS_exp_nls          := maxRSS ]
plot_dt[ is.na( RSS_exp_nls_noMahDis ) , RSS_exp_nls_noMahDis := maxRSS ]
plot_dt[ is.na( RSS_sig_nls )          , RSS_sig_nls          := maxRSS ]
plot_dt[ is.na( RSS_sig_nls_noMahDis ) , RSS_sig_nls_noMahDis := maxRSS ]

# Which % of pairs is better explained by exponential or sigmoid models?
# "Better explained" is defined as a >= 10% reduction in RSS
dt_nls <- plot_dt[, .(pc_exp_better_than_lin = sum( (RSS_exp_nls / RSS_lin_nls) <= 0.9 ) / .N * 100 ,
                      pc_sig_better_than_lin = sum( (RSS_sig_nls / RSS_lin_nls) <= 0.9 ) / .N * 100 ,
                      pc_exp_better_than_lin_noMahOut = sum( (RSS_exp_nls_noMahDis / RSS_lin_nls_noMahDis) <= 0.9 ) / .N * 100 ,
                      pc_sig_better_than_lin_noMahout = sum( (RSS_sig_nls_noMahDis / RSS_lin_nls_noMahDis) <= 0.9 ) / .N * 100 ),
                  by = .(type, group) ]
dt_nls <- melt(dt_nls, id.vars = c("type", "group"))

# Change plotting order
dt_nls[, variable := factor(variable, levels = c("pc_exp_better_than_lin", "pc_exp_better_than_lin_noMahOut", "pc_sig_better_than_lin", "pc_sig_better_than_lin_noMahout"))]

# Create a barplot
pNLS <- ggplot( dt_nls, aes( x = group, y = value, fill = variable, colour = variable))+
  facet_wrap(~type)+
  geom_bar( stat="identity", position = "dodge")+                   
  ylab("Protein pairs[%]")+
  scale_fill_manual(   values = c("darkred", "red", "turquoise4", "turquoise2") )+
  scale_colour_manual( values = c("darkred", "red", "turquoise4", "turquoise2") )+
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA, colour="black", size=0.25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top",
        axis.text=element_text(size=5), axis.title=element_text(size=6), axis.ticks.y = element_line(size=0.25),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())


#### Plots for Figure 3 ######

Fig3 <- ggdraw()+
          draw_plot( pTP_scatter, x = 0,    y = 0.75, width = 0.25, height = 0.25)+
          draw_plot( pFP_scatter, x = 0.25, y = 0.75, width = 0.25, height = 0.25)+
          draw_plot( pPie       , x = 0.0,  y = 0.50, width = 0.5,  height = 0.25)+
          draw_plot( pOut_1     , x = 0.5,  y = 0.75, width = 0.25, height = 0.25)+
          draw_plot( pOut_2     , x = 0.75, y = 0.75, width = 0.25, height = 0.25)+
          draw_plot( pOut_3a    , x = 0.5 , y = 0.50, width = 0.25, height = 0.25)+
          draw_plot( pOut_4a    , x = 0.75, y = 0.50, width = 0.25, height = 0.25)+
          draw_plot( pOut_3b    , x = 0.5 , y = 0.25, width = 0.25, height = 0.25)+
          draw_plot( pOut_4b    , x = 0.75, y = 0.25, width = 0.25, height = 0.25)+
          draw_plot( pOut_3c    , x = 0.5 , y = 0.00, width = 0.25, height = 0.25)+
          draw_plot( pOut_4c    , x = 0.75, y = 0.00, width = 0.25, height = 0.25)

save_plot("Figure_3.pdf", Fig3, ncol = 4, nrow = 4, base_height = 1.8, base_aspect_ratio = 1)


#### Plots for Figure 4 ######

Fig4 <- ggdraw()+
          draw_plot( pD_TP  , x = 0   , y = 0.5, width = 0.33, height = 0.5)+
          draw_plot( pD_FP  , x = 0.33, y = 0.5, width = 0.33, height = 0.5)+
          draw_plot( pEx_bar, x = 0.66, y = 0.5, width = 0.33, height = 0.5)+
          draw_plot( pMAE   , x = 0   , y = 0.0, width = 1.00, height = 0.5)

save_plot("Figure_4.pdf", Fig4, ncol = 3, nrow = 2, base_height = 1.8, base_aspect_ratio = 1.04)

  
#### Plots for Supplementary Figure S1 ######

FigS1 <- ggdraw()+
          draw_plot( pExp  , x = 0   , y = 0.5, width = 0.5, height = 0.5)+
          draw_plot( pSig  , x = 0.5 , y = 0.5, width = 0.5, height = 0.5)+
          draw_plot( pNLS  , x = 0   , y = 0  , width = 1  , height = 0.5)

save_plot("Figure_S1.pdf", FigS1, ncol = 2, nrow = 2, base_height = 1.8, base_aspect_ratio = 1)


