
cat("|| ============================= || \n")
cat("Starting permutation test... \n")

# Load libraries
library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv
library(tidyverse)
library(parallel)
library(doParallel)

# Read in command line arguments
args<-commandArgs(trailingOnly=TRUE)

cat("the arguments are:")
print(args)

n_snps_values <- as.numeric(args[1])
data_file <- as.character(args[2])

# Load data
load(data_file)



## Parallelize calculations across main calculations and folds
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

cat("Using this number of cores:", cores, " ... \n")

# Set up model terms
fixed_effects_resids <- paste0(paste0("rgr_resids ~ s(gebv_by_count, bs=\"cr\") + 
                                               s(tmax_sum_dif, by = gebv_by_count, bs=\"cr\")"))

fixed_effects_resids_permuted <- paste0(paste0("rgr_resids_permuted ~ s(gebv_by_count, bs=\"cr\") + 
                                                 s(tmax_sum_dif, by = gebv_by_count, bs=\"cr\")"))


# Calculate polygenic risk score for each maternal genotype
gen_dat_moms <- dplyr::select(gen_dat_clim, id, accession, snp_col_names)
dim(gen_dat_moms)


# Define function that calculates GEBV 
calc_gebv <- function(fold_input = NA, parallel = TRUE){
  
  # Choose SNPs
  n_snps_values <- n_snps_values

  # If not in parallel
  if(parallel == FALSE){
    fold = fold_input
  }
  
  # Initialize dataframe for output
  gebv_df_summary <- tibble(n_snps = n_snps_values,
                            fold = fold,
                            r2 = NA,
                            dev = NA,
                            high_mid = NA,
                            high_low = NA)
  gebv_values <- tibble(accession = gen_dat_moms$accession)
  
  # Order snps based on lowest P values
  if(fold == 0 ){ # If main model with full dataset
    snp_list_ordered <- sum_df$snp[order(sum_df$p_val_gen_int)]
  } else {
    snp_list_ordered <- sum_df$snp[order(sum_df[paste0("p_val_gen_int_fold", fold)])]
  }
  
  # Randomizing order of P values as suggested by reviewers
  # if(fold == 0 ){ # If main model with full dataset
  #   snp_list_ordered <- sum_df$snp[sample(1:length(sum_df$p_val_gen_int))]
  # } else {
  #   snp_list_ordered <- sum_df$snp[sample(1:length(sum_df$p_val_gen_int))]
  # }
  
  # Initialize
  gebv <- NULL
  gebv_count <- NULL
  
  # Initialize counter
  x = 1
  
  for(snp_name in snp_list_ordered){
    
    if(x %% 500 == 0){
      cat("Working on snp number: ", x, "... \n")
    }
    
    # Subset to genotypes
    sub <- gen_dat_moms[ , snp_name]
    
    # If main model with full dataset
    if(fold == 0 ){ 
      preds_sub <- dplyr::filter(sum_df_long_w_preds, snp == snp_name) %>%
        dplyr::select(genotype, height_change_warmer_absolute)
    } else {
      preds_sub <- dplyr::filter(sum_df_long_w_preds, snp == snp_name) %>%
        dplyr::select(genotype, paste0("height_change_warmer_absolute_fold", fold)) %>%
        dplyr::rename(height_change_warmer_absolute =  
                        paste0("height_change_warmer_absolute_fold", fold))
    }
    
    # Remove NAs
    preds_sub <- preds_sub %>%
      filter(!is.na(height_change_warmer_absolute))
    
    
    counted <- preds_sub$genotype
    
    # Add in median value for missing genotypes
    if(!"0" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "0",
                                               height_change_warmer_absolute = 0))
    }
    if(!"1" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "1",
                                               height_change_warmer_absolute = 0))
    }
    if(!"2" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "2",
                                               height_change_warmer_absolute = 0))
    }
    
    gebv_temp <-  case_when(
      sub == "0" ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "0"],
      sub == "1"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "1"],
      sub == "2"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "2"],
      #  is.na(sub) ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "missing"],
      TRUE ~ 0 # When individual has missing data
    )
    
    gebv_count_temp <- ifelse(as.data.frame(sub)[, 1] %in% counted, yes = 1, no = 0)
    
    if(x == 1){
      gebv <- gebv_temp
      gebv_count <- gebv_count_temp
    } else {
      gebv <- gebv + gebv_temp
      gebv_count <- gebv_count + gebv_count_temp
    }
    
    # Error checking
    if(any(is.na(gebv))){
      stop("NAs detected in GEBV calculation")
    }
    
    if(x %in% n_snps_values | x == length(snp_list_ordered)){
      cat("\n||--------------------------------||\n")
      cat("Running model for:", x, " SNPS... \n")
      
      # Adjust for missing data
      gebv_by_count <- gebv / gebv_count
      # gebv_by_count <- gebv
      
      # Save to output DF
      gebv_values[, paste0("gebv_fold",fold, "_", x)] <- gebv_by_count
      
      ### Join to Testing dataframe
      if(fold == 0 ){
        dat_snp_testing <- dat_snp_all[, c("accession" ,"rgr_resids",
                                           "tmax_sum_dif")]
      } else {
        dat_snp_testing <- dat_snp_all[-training_folds[[fold]], c("accession" ,"rgr_resids",
                                                                  "tmax_sum_dif")]
      }
      
      # Convert accession to numeric
      dat_snp_testing$accession <- as.numeric(as.character(dat_snp_testing$accession))
      
      # Join gebv to testing dataframe
      dat_snp_testing <- dplyr::left_join(dat_snp_testing, 
                                          cbind(gen_dat_moms[ ,"accession"], gebv_by_count))  
      
      # Scale GEBV
      dat_snp_testing$gebv_by_count <- c(scale(dat_snp_testing$gebv_by_count))
      
      # Gam with interaction effect
      gam_test_resids = bam(formula = formula(fixed_effects_resids),
                            data = dat_snp_testing,
                            discrete = FALSE, 
                            nthreads = 8,
                            method = "fREML")
      
      print(summary(gam_test_resids))
      
      ## Save results   
      gebv_df_summary$r2[gebv_df_summary$n_snps == x] <- summary(gam_test_resids)$r.sq
      gebv_df_summary$dev[gebv_df_summary$n_snps == x] <- summary(gam_test_resids)$dev.expl
      
      cat(":::::::::::::::::::::::::::::\n")
      cat("Working on permutations... \n ")
      ## Permutation test
      for(perm in 1:500){
        
        if(perm %% 50 == 0){
          cat("Working on permutation: ", perm, "... \n")
        }
        
        # Permute phenotype
        dat_snp_testing$rgr_resids_permuted <- sample(dat_snp_testing$rgr_resids)
        
        # Gam with interaction effect
        gam_test_resids_permuted = bam(formula = formula(fixed_effects_resids_permuted),
                                       data = dat_snp_testing,
                                       discrete = FALSE, 
                                       nthreads = 8,
                                       method = "fREML")
        
        gebv_df_summary[gebv_df_summary$n_snps == x, paste0("r2_permuted_", perm)] <- summary(gam_test_resids_permuted)$r.sq
        
      }
      
      
    } # End model loop
    
    x = x + 1 # Move onto next SNP
    
  } # End snp loop
  
  return(list(gebv_df_summary, gebv_values)) # Return summary df as output of loop
  
} # End calc_gebv function

# Run parallelized funtion to calculate gebv for each fold and 
   gebv_df_summary_list <- foreach(fold=0:10, 
                                   .packages = c("tidyverse", "mgcv"),
                                   .combine = rbind,
                                   .export = "n_snps_values") %dopar% {
                                     calc_gebv() # call function
                                   }

  ## Non parallel version
 # for(fold in 0:10){
 #   cat("Working on fold: ", fold, " ... \n")

 #   if(fold == 0){
 #     gebv_df_summary_list <- calc_gebv(fold = fold, parallel = FALSE)
 #   } else {
 #     gebv_df_summary_list <- rbind(gebv_df_summary_list, 
 #                             calc_gebv(fold = fold, parallel = FALSE))
 #   }

 # }


   
   # Quick helper function
   bind_tibbles <- function(list, mode){
     for(i in 1:length(list)){
       if(i == 1){
         out <- list[[i]]
       } else {
         if(mode == "rbind"){
           out <- rbind(out, list[[i]])
         } else if(mode == "cbind"){
           out <- cbind(out, list[[i]])
         }
       }
     }
     return(out)
   }
   
   gebv_df_summary <- bind_tibbles(gebv_df_summary_list[, 1], mode = "rbind")
   
   # Save output
   write_csv(gebv_df_summary, paste0("gebv_permutations_", n_snps_values,"_", Sys.Date(), ".csv"))
