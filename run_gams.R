# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

	task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])
  climate_var_dif <- as.character(args[3])

# Load libraries
	library(mgcv)
	library(tidyverse)
  library(visreg) # Installed on home directory

# Load data
 load("../gam_cluster_2018-05-31.Rdata")


 # Initialize variables
  gam_list <- list()
  x = 1
  
  ## Convert fixed and random parts of formula to actual formula objetc
  # No idea why we need to call formula twice, but it works
  make_formula <- function(fixed, random){
    return(formula(formula(enquote(paste0(fixed, random)))))
  }
  
   # Formula for random effects
  random_effects <-  '+ s(accession, bs = "re") + s(section_block, bs = "re")'
  



  # Loop through snps by index based on task id
  for(snp_index in task_id:(task_id + interval - 1)){
        
    # To avoid going past number of snps
    if(snp_index > length(snp_col_names)){
      next
    }
    
    # For timing loops
    start_time <- Sys.time()
    
    # Choose snp
    snp <- paste0("snp_", snp_col_names[snp_index])
    
    # Make sure there's at least 3 levels
      if(length(table(dat_snp[, snp])) < 3){
        cat("Skipping because not at least 3 levels")
        next
      }

    # Choose numbers of knots
    k = length(table(pull(dat_snp, snp)))
    
    # Convert to factor
    # dat_snp[, snp] <- as.factor(pull(dat_snp, snp))
    
  
      
    
  # For SNPS as factors    
     
    # Stack overflow on adding m = 1 to by= smooths - https://stats.stackexchange.com/questions/32730/how-to-include-an-interaction-term-in-gam
     
    # fixed_effects_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by =",snp,", bs=\"cr\", m = 1) + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    # fixed_effects_no_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
  # For SNPS as continuous  
    fixed_effects_int <- paste0(paste0("height_2017 ~ site + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\") + ti(", climate_var_dif,", ",snp,", bs=\"cr\", k =", k," ) + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    fixed_effects_no_int <- paste0(paste0("height_2017 ~ site + s(height_2014, bs=\"cr\") + te(", snp, ", bs=\"cr\", k = ", k, ") +  te(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
  # Run gams
  
   cat("Working on: ", snp, "... number: ", x, " ...\n" )

    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp, 
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML", family = "tw")
    
    gam_snp_no_int =  bam(formula = make_formula(fixed_effects_no_int,
                                                 random_effects),
                          data = dat_snp, 
                          discrete = TRUE,
                          nthreads = 8,
                          method = "fREML", family = "tw")
    
    # Compare AICs
    gam_snp_int$delta_aic <- gam_snp_int$aic - gam_snp_no_int$aic



   # Visualization plots
    
    # For snps as factors
      # pdf(paste0("./model_plots/", snp,".pdf"), width = 8, height = 5)
      # visreg(gam_snp_int, xvar = climate_var_dif, 
      #        by = snp, scale = "response")
      # dev.off()
    
    # For SNPS as continuous
        pdf(paste0("./model_plots/", snp,".pdf"), width = 5, height = 5)
        visreg2d(gam_snp_int, 
                 xvar = climate_var_dif, yvar = snp, 
                 scale = "response", 
                 plot.type = "persp", theta = 35)
        dev.off()


  # Make predictions

   # When SNPs are continuous, make own genotypes
      gen_temp <- data.frame(accession = "1", 
                 genotype = c(0, 1, 2),             
               genotype_scaled = c((0 - scaled_snps_means[snp]) / scaled_snps_sds[snp],
                              (1 - scaled_snps_means[snp]) / scaled_snps_sds[snp],
                              (2 - scaled_snps_means[snp]) / scaled_snps_sds[snp]))
                 
      
      predictions_3deg_sub <- left_join(predictions_3deg, gen_temp)
        
        
    # Make sure genotypes are represented and rename column      
    predictions_3deg_sub <-  predictions_3deg_sub  %>%
      dplyr::filter(genotype_scaled %in% pull(dat_snp, snp)) %>% 
      dplyr::mutate(!!snp := genotype_scaled)
    
  # For SNPs as factors  
    # Change accession number if genetic data is missing so we don't get errors in prediction
    # while(any(is.na(dat_snp[as.character(dat_snp$accession) %in% 
    #                         as.character(predictions_3deg_sub$accession), snp]))){
    #   predictions_3deg_sub$accession = as.numeric(predictions_3deg_sub$accession) + 1
    # }
    
      
    # Subset predictions based on 
    predictions_3deg_sub$pred <- predict(gam_snp_int, 
                                         newdata = predictions_3deg_sub,
                                     type = "response")
    
    
    climate_min <- min(predictions_3deg_sub$tmax_sum_dif)
    climate_max <- max(predictions_3deg_sub$tmax_sum_dif)
    
    for(genotype in c(0,1,2)){
      
      if(!genotype %in% predictions_3deg_sub$genotype){
        gam_snp_int[paste0("height_change_gen_", genotype)] <- NA
        next
      }
      
      pred_sub <- predictions_3deg_sub[predictions_3deg_sub$genotype == genotype, ]
      
      height_change <- (pred_sub$pred[pred_sub$tmax_sum_dif == climate_max] - 
        pred_sub$pred[pred_sub$tmax_sum_dif == climate_min]) / pred_sub$pred[pred_sub$tmax_sum_dif == climate_max] * 100 
      
      gam_snp_int[paste0("height_change_gen_", genotype)] <- height_change
      
      }

    # Attach data so that we can use visreg later if we need
    # gam_snp$data <- dat_snp
    

    # Save into list
    gam_list[[x]] <- gam_snp_int
    names(gam_list)[x] <- snp
    
    x = x + 1
    
    
    # Save individual .Rdata file for each model
    cat("Saving gam model to file... \n")

    gam_name <- paste0("gam_", climate_var_dif, "_", snp)
    assign(gam_name, gam_snp_int) # Assign model the name

    # save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))
    
    # Timing of loop
    end_time <- Sys.time()
    cat("This loop took: ")
    print(end_time - start_time)   
    
  } # End SNP loop


  
  
  # Run summary for each model
  gam_list_summary <- lapply(gam_list, summary)
  
  # Save into dataframe 
  gam_mod_summary_df <- data.frame(
    
    # climate var
    climate_var = climate_var_dif,
    
    # Snp name
    snp = names(gam_list),
    
    # Sample size
    n = unlist(lapply(gam_list_summary, function(x) x$n)),

    # Converged?
    converged = unlist(lapply(gam_list, function(x) x$mgcv.conv)),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # Estimated degrees of freedom
    edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
    
    # AIC
    aic = unlist(lapply(gam_list, function(x) x$aic)),
    
    # Delta AIC with model without interaction
    delta_aic = unlist(lapply(gam_list, function(x) x$delta_aic)),

    # Predictions of height
    height_change_gen_0 = unlist(lapply(gam_list, function(x) x$height_change_gen_0)),
    height_change_gen_1 = unlist(lapply(gam_list, function(x) x$height_change_gen_1)),
    height_change_gen_2 = unlist(lapply(gam_list, function(x) x$height_change_gen_2)),

    # P values for interaction terms
    p_val_int = ifelse(unlist(lapply(gam_list_summary, function(x)
      length(x$s.table[grep(pattern = "^ti", 
                            rownames(x$s.table)), "p-value"]) == 0)), 
      yes = NA, 
      no = unlist(lapply(gam_list_summary, function(x)
        x$s.table[grep(pattern = "^ti",
                       rownames(x$s.table)), "p-value"]))),
    
    # Use grep to look at the last number of the smooth summary table - 
    # Should correspond to the genotype
    # Returns NA if no match
    p_val_gen_0 = ifelse(unlist(lapply(gam_list_summary, function(x)
                                length(x$s.table[grep(pattern = "0$",
                                               rownames(x$s.table)), "p-value"]) == 0)), 
                         yes = NA, 
                         no = unlist(lapply(gam_list_summary, function(x)
                                                    x$s.table[grep(pattern = "0$",
                                                   rownames(x$s.table)), "p-value"]))),
    p_val_gen_1 = ifelse(unlist(lapply(gam_list_summary, function(x)
                          length(x$s.table[grep(pattern = "1$",
                            rownames(x$s.table)), "p-value"]) == 0)), 
                              yes = NA, 
                           no = unlist(lapply(gam_list_summary, function(x)
                                  x$s.table[grep(pattern = "1$",
                                  rownames(x$s.table)), "p-value"]))),
    
    p_val_gen_2 = ifelse(unlist(lapply(gam_list_summary, function(x)
                           length(x$s.table[grep(pattern = "2$",
                            rownames(x$s.table)), "p-value"]) == 0)), 
                          yes = NA, 
                          no = unlist(lapply(gam_list_summary, function(x)
                          x$s.table[grep(pattern = "2$",
                                         rownames(x$s.table)), "p-value"])))
    
      )
  
  head(gam_mod_summary_df)
  
 
 write.csv(gam_mod_summary_df, file = paste0("./model_summaries/gam_summaries_", task_id, ".csv"),
                               row.names = FALSE)


###### LAGNIAPPE








