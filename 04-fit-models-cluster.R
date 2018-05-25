

# Load libraries ----------------------------------------------------------

# install_github("jdstorey/qvalue")
library(qvalue)
library(MuMIn)
library(vegan)


# Goal is to run an individual gam for each SNP, like a fancy GWAS, and then correct for multiple testing after

# Format and save to .Rdata for cluster -----------------------------------

  # Joining climate and genetic data
  # Do an inner join to avoid having GBS moms with climate data but poor genotyping data
  gen_dat_clim <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim)
  
  # Getting genetic data
  snp_col_names <- colnames(gen_dat[, -c(1, 2)])
  
  # Calculate percent missing data  
  missing <- gen_dat_clim %>%
    dplyr::select(id, accession) %>%
    mutate(missing = apply(gen_dat_clim[, snp_col_names], 1, function(x) sum(is.na(x))/length(x)))
  
  missing %>%
    arrange(desc(missing))

  # Scale data and impute missing before converting 012 to -101
  bglr_gen_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
  bglr_gen_scaled[1:10, 1:10]
  summary(bglr_gen_scaled[, 1:10])
  
  
  ## Calculate PCA on genetic data
  pca_gen = prcomp(bglr_gen_scaled, center = FALSE, scale = FALSE)
  summary(pca_gen)
  screeplot(pca_gen, bstick = TRUE)
  
  pcs <- as.data.frame(pca_gen$x[, 1:10])
  colnames(pcs) <- paste0(colnames(pcs), "_gen")
  
  bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

  
  # Combine with genetic data 
    dat_snp = left_join(dat_gbs_only_scaled, 
                        bind_cols(gen_dat_clim[, "accession"],
                                           bglr_gen_scaled), 
                        by = "accession")
    
    dat_snp$accession <- factor(dat_snp$accession)
    dat_snp$section_block <- factor(dat_snp$section_block)
    dat_snp$section <- factor(dat_snp$section)
    
    
    # Adding snp_ to column names - useful for gams
    colnames(dat_snp)[which(colnames(dat_snp) %in% snp_col_names)] <- paste0("snp_", colnames(dat_snp[, snp_col_names]))
 
# Set a prediction data frame
  
  # function to transform from raw to scaled variables
 #  forward_transform <- function(x, var, means, sds){
 #   ( x - means[var]) / sds[var]
 #  }
 # 
 # # Prediction frame with a 3 degree increase 
 #  # Set up each climate variable individually
 #  predictions_3deg <-  expand.grid(height_2014 = 0,
 #                  site = "Chico",
 #                  accession = "1",
 #                  section_block = "Block1_1",
 #                  tmax_sum_dif = forward_transform(c(0, 3), 
 #                                                   var = "tmax_sum_dif", 
 #                                                   means = scaled_var_means_gbs_only,
 #                                                   sds = scaled_var_sds_gbs_only),
 #                  tmin_winter_dif = forward_transform(c(0, 3), 
 #                                                      var = "tmin_winter_dif", 
 #                                                      means = scaled_var_means_gbs_only,
 #                                                      sds = scaled_var_sds_gbs_only))
 #  
 #  predictions_3deg

  
# Save data to file that will be uploaded to cluster  
save(dat_snp, snp_col_names,
  file = paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))

  
# Run baseline models -----------------------------------------------------


  base_gams <- list()
  x = 1
  
  ## Convert fixed and random parts of formula to actual formula objetc
  # No idea why we need to call formula twice, but it works
  make_formula <- function(fixed, random){
    return(formula(formula(enquote(paste0(fixed, random)))))
  }
  
  for(climate_var_dif in c("tmax_sum_dif", "tmin_winter_dif")){
 
    # Formula for fixed effects
    fixed_effects <- paste0(paste0("height_2017 ~ site + s(height_2014) + s(", 
                                   climate_var_dif,")"))
    
    # Formula for random effects
    random_effects <-  '+ s(accession, bs = "re") + s(section_block, bs = "re")'
    
    
    gam_base <- bam(formula = make_formula(fixed_effects, random_effects),
                              data = dat_snp, 
                              discrete = TRUE,
                              nthreads = 8,
                              method = "fREML", family = "tw")
    
    base_gams[[x]] <- gam_base
    names(base_gams)[x] <- climate_var_dif
    
    x = x + 1
  
  }
  
# Run summary for each model
  base_gams_summary <- lapply(base_gams, summary)
  
  # Save into dataframe 
  gam_mod_summary_df <- data.frame(
    # climate var
    climate_var = names(base_gams),
    
    # Snp name
    snp = "base",
    
    # Deviance explained
    dev_explained = unlist(lapply(base_gams_summary, function(x) x$dev.expl)),
    
    # Estimated degrees of freedom
    edf = unlist(lapply(base_gams, function(x) sum(x$edf))),
    
    # AIC
    aic = unlist(lapply(base_gams, function(x) x$aic)),
    
    # P value of interaction term
    p_val_int = NA
    
    # Percent change in height
    # height_change_3_deg = unlist(lapply(lapply(base_gams, function(x) predict(x, 
    #                                              newdata = predictions_3deg  %>% 
    #                                                dplyr::distinct(get(climate_var_dif), 
    #                                                                .keep_all = TRUE), 
    #                                              type = "response")),
    #        function(x) ((x[2] - x[1]) / x[1])*100))
    )

  head(gam_mod_summary_df)
  
  write.csv(gam_mod_summary_df, file = paste0("./output/model_summaries_test/gam_summaries_base_models.csv"),
            row.names = FALSE)
  
  
  

# Code from cluster for testing -------------------------------------------

  
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
    
  
    # Formula for fixed effects
    fixed_effects <- fixed_effects <- paste0(paste0("height_2017 ~ site + 
                                    s(height_2014) + te(", 
                                  climate_var_dif,") + te(",snp,",k = 3) + ti(",
                                  climate_var_dif, ",", snp,", k = 3) +
                                  s(PC1_gen) + s(PC2_gen) + s(PC3_gen)"))
  
    
    cat("Working on: ", snp, "... number: ", x, " ...\n" )
    
    gam_snp = bam(formula = make_formula(fixed_effects, random_effects),
                  data = dat_snp, 
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML", family = "tw")
    
    summary(gam_snp)
    visreg(gam_snp, scale = "response")

    # Attach data so that we can use visreg later if we need
    # gam_snp$data <- dat_snp
    
    # Save into list
    gam_list[[x]] <- gam_snp
    names(gam_list)[x] <- snp
    
    x = x + 1
    
    
    # Save individual .Rdata file for each model
    cat("Saving gam model to file... \n")
    
    gam_name <- paste0("gam_", climate_var_dif, "_", snp)
    assign(gam_name, gam_snp) # Assign model the name
    
    save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))
    
    # Timing of loop
    end_time <- Sys.time()
    
    cat("This loop took: ")
    print(end_time - start_time)   
    
  }
  
  
  
  
  # Run summary for each model
  gam_list_summary <- lapply(gam_list, summary)
  
  # Save into dataframe 
  gam_mod_summary_df <- data.frame(
    # climate var
    climate_var = climate_var_dif,
    
    # Snp name
    snp = names(gam_list),
    
    # Deviance explained
    dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # Estimated degrees of freedom
    edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
    
    # AIC
    aic = unlist(lapply(gam_list, function(x) x$aic)),
    
    # P value of interaction term
    
    ## IF ORDER OF VARIABLES IN MODEL EVER CHANGES - THIS NEEDS TO CHANGE TOO!!
    
    p_val_int = unlist(lapply(gam_list_summary, function(x) x$s.table[4, "p-value"]))
    
    # Percent change in height with 3 degree increase
    # height_change_3_deg = unlist(lapply(lapply(gam_list, function(x) predict(x, 
    #                                                                          newdata = predictions_3deg  %>% 
    #                                                                            dplyr::distinct(get(climate_var_dif), 
    #                                                                                            .keep_all = TRUE), 
    #                                                                          type = "response")),
    #                                     function(x) ((x[2] - x[1]) / x[1])*100))
    
  )
  
  head(gam_mod_summary_df)

# Read in cluster run output ----------------------------------------------

  path_to_summaries <- "./output/model_summaries_tmax_sum/"  
  
  sum_df <- ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df)
  sum_df
  
  summary(sum_df)

# Sort and arrange data ---------------------------------------------------


# Calculate AIC weights
  sum_df <- sum_df %>%
    dplyr::group_by(climate_var) %>%
    dplyr::mutate(aic_weight = Weights(aic))

  
# Look at top SNPs
    
  # Sort by lowest AIC  
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(aic)
    
  # Sort by lowest p value  
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(p_val_int)
    
  # Sort by highest deviance explained
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(desc(dev_explained))
    
    summary(sum_df$dev_explained)
    
  
## Calculating q value  
  # Vignette - https://bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
  
  hist(sum_df$p_val_int, nclass = 20)
  
  qvals <- qvalue(p = sum_df %>%
                    dplyr::filter(climate_var == "tmax_sum_dif") %>%
                    dplyr::pull(p_val_int))
  
  sum_df <- sum_df %>%
    dplyr::filter(climate_var == "tmax_sum_dif") %>%
    mutate(q_val = qvals$qvalues)
  
  qvals$qvalues
  qvals$lfdr
  qvals$pi0 # overall proportion of true null hypotheses
  
  max(qvals$qvalues[qvals$pvalues <= 0.01]) # FDR for p values below threshold
  
  summary(qvals)
  
  hist(qvals)
  
  plot(qvals)
  
  # Setting an alpha threshold 
  qobj_fdrlevel <- qvalue(p = sum_df$p_val_int, fdr.level = 0.05)
  qobj_fdrlevel$significant
  table(qobj_fdrlevel$significant)
  
  summary(qobj_fdrlevel)
  
  # Sort by lowest q value  
  sum_df %>%
    dplyr::group_by(climate_var) %>%
    dplyr::arrange(q_val)
  
  

# Manhattan plots ---------------------------------------------------------


sum_df_sub <-  sum_df %>%
  dplyr::filter(climate_var == "tmax_sum_dif")    
    
  plot(-log10(sum_df_sub$p_val_int), pch = 20,
       xlab = "SNP", ylab = "-log10(P)", 
       main = "tmax_sum_dif", las = 1,
       col = factor(snp_pos$chrom %% 2), # To alternate colors for pseudo-chromosomes
       ylim = c(0, 9))
  abline(h = -log10(5e-8), lty = 2, col = "red") # Genome-wide significance line
  abline(h = -log10(1e-5)) # Suggestive significance line  
  
  

# Look at top models ------------------------------------------------------

  top_snp = "snp_7_39668055"
  
  plot(pull(dat_snp, top_snp), dat_snp$height_2017)
  plot(pull(dat_snp, tmax_sum_dif),pull(dat_snp, top_snp))
  
  summary(gam_tmax_sum_dif_snp_7_39668055)
  visreg(gam_tmax_sum_dif_snp_7_39668055, scale = "response")
  visreg2d(gam_tmax_sum_dif_snp_7_39668055, xvar = "tmax_sum_dif", 
           yvar = "snp_7_39668055", 
           scale = "response", plot.type = "persp", theta = 45, phi = 15)
  
  pairs.panels(dplyr::select(dat_snp, top_snp, latitude, longitude, climate_vars ))
  
