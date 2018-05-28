

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
    # dat_snp = left_join(dat_gbs_only_scaled, 
    #                     bind_cols(gen_dat_clim[, "accession"],
    #                                        bglr_gen_scaled), 
    #                     by = "accession")
  
  # For using snps as factors
  dat_snp = left_join(dat_gbs_only_scaled, 
                    bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                       pcs),
                    by = "accession")
  
    # dat_snp <- dat_snp %>%
    #   dplyr::mutate_at(snp_col_names, funs(factor(.)))
    #   
    dat_snp$accession <- factor(dat_snp$accession)
    dat_snp$section_block <- factor(dat_snp$section_block)
    dat_snp$section <- factor(dat_snp$section)
    
    
    # Adding snp_ to column names - useful for gams
    colnames(dat_snp)[which(colnames(dat_snp) %in% snp_col_names)] <- paste0("snp_", colnames(dat_snp[, snp_col_names]))
 
# Set a prediction data frame
  
  # function to transform from raw to scaled variables
 #                                                      sds = scaled_var_sds_gbs_only))
  forward_transform <- function(x, var, means, sds){
   ( x - means[var]) / sds[var]
  }

 # Prediction frame with a 3 degree increase
  # Set up each climate variable individually
  predictions_3deg <-  expand.grid(height_2014 = 0,
                  site = "Chico",
                  accession = "1",
                  section_block = "Block1_1",
                  PC1_gen = 0, PC2_gen = 0, PC3_gen = 0,
                  genotype = c(0, 1, 2),
                  tmax_sum_dif = forward_transform(c(0, 3),
                                                   var = "tmax_sum_dif",
                                                   means = scaled_var_means_gbs_only,
                                                   sds = scaled_var_sds_gbs_only))

  predictions_3deg

  
# Save data to file that will be uploaded to cluster  
save(dat_snp, snp_col_names, predictions_3deg,
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
                                   climate_var_dif,") + s(PC1_gen) + s(PC2_gen) + s(PC3_gen)"))
    
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
    
    # Rsquared adjusted
    r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
    
    # AIC
    aic = unlist(lapply(base_gams, function(x) AIC(x))),
    
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
    # if(length(table(dat_snp[, snp])) < 3){
    #   cat("Skipping because not at least 3 levels")
    #   next
    # }
    
    # Convert to factor
    dat_snp[, snp] <- as.factor(pull(dat_snp, snp))
    
  
    # Formula for fixed effects
    fixed_effects_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\") + s(", climate_var_dif,", by =",snp,", bs=\"cr\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
    
    fixed_effects_no_int <- paste0(paste0("height_2017 ~ site + ", snp, " + s(height_2014, bs=\"cr\") + s(", climate_var_dif,", bs=\"cr\")  + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\")"))
  
    
    cat("Working on: ", snp, "... number: ", x, " ...\n" )
    
    gam_snp_int = bam(formula = make_formula(fixed_effects_int, random_effects),
                  data = dat_snp, 
                  discrete = TRUE,
                  nthreads = 8,
                  method = "fREML", family = "tw", 
                  control = list(trace = FALSE))
    
    gam_snp_no_int =  bam(formula = make_formula(fixed_effects_no_int,
                                                 random_effects),
                          data = dat_snp, 
                          discrete = TRUE,
                          nthreads = 8,
                          method = "fREML", family = "tw")
    
    # Compare AICs
    gam_snp_int$delta_aic <- gam_snp_int$aic - gam_snp_no_int$aic
    
    # Make predictions
    predictions_3deg <-  predictions_3deg  %>%
      dplyr::mutate(!!snp := genotype)
    
    predictions_3deg$pred <- predict(gam_snp_int, newdata = predictions_3deg,
                                     type = "response")
    
    
    climate_min <- min(predictions_3deg$tmax_sum_dif)
    climate_max <- max(predictions_3deg$tmax_sum_dif)
    
    for(genotype in c(0,1,2)){
      
      pred_sub <- predictions_3deg[predictions_3deg$genotype == genotype, ]
      
      height_change <- (pred_sub$pred[pred_sub$tmax_sum_dif == climate_max] - 
        pred_sub$pred[pred_sub$tmax_sum_dif == climate_min]) / pred_sub$pred[pred_sub$tmax_sum_dif == climate_max] * 100 
      
      gam_snp_int[paste0("height_change_gen_", genotype)] <- height_change
      
      }

   
    # summary(gam_snp_int)
    # visreg(gam_snp_int, scale = "response")
    # visreg(gam_snp_int, xvar = "tmax_sum_dif", by = snp, scale = "response")
    # visreg2d(gam_snp, xvar = "tmax_sum_dif", 
    #          yvar = snp, 
    #          scale = "response", 
    #          plot.type = "persp", theta = 45, phi = 15)

    # Attach data so that we can use visreg later if we need
    # gam_snp$data <- dat_snp
    
    # Save into list
    gam_list[[x]] <- gam_snp_int
    names(gam_list)[x] <- snp
    
    x = x + 1
    
    
    # Save individual .Rdata file for each model
    # cat("Saving gam model to file... \n")
    # 
    # gam_name <- paste0("gam_", climate_var_dif, "_", snp)
    # assign(gam_name, gam_snp_int) # Assign model the name
    # 
    # save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))
    
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
    height_change_gen_2 = unlist(lapply(gam_list, function(x) x$height_change_gen_2))
  )
  
  head(gam_mod_summary_df)
  
  
  
  

# Read in cluster run output ----------------------------------------------

  path_to_summaries <- "./output/model_summaries_tmax_sum/"  
  
  sum_df <- ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df)
  sum_df
  
  summary(sum_df)

# Sort and arrange data ---------------------------------------------------

  # How many converged?
  table(sum_df$converged)
  
  # Sort by delta AIC  
  test = sum_df %>%
    dplyr::group_by(climate_var) %>%
    dplyr::filter(converged == TRUE) %>%
    dplyr::filter(delta_aic > -1000) %>% 
    dplyr::filter(delta_aic < 500) %>%
    dplyr::arrange(delta_aic)
  
  
  hist(test$delta_aic, breaks = 30)
  
  
  
  
  
  # Take snps with lowest delta aic 
  # Filter out ones that did not converge
  # - run in same model to generate super model?
  # - filter down to snps that aren't super correlated with each other?
  # - try to figure out more about those SNPs with gradient forest / RDA, etc?
  
  
  
  
  
  
  
  

# Calculate AIC weights
  sum_df <- sum_df %>%
    dplyr::group_by(climate_var) %>%
    dplyr::mutate(aic_weight = Weights(aic))
  
  # sum_df <- sum_df %>%
  #   dplyr::mutate(aic_delta = aic - 28941.18) # AIC calculated from base model
  # 
  # 
# Look at top SNPs
    
  # Sort by lowest AIC  
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(aic)
    
  # Sort by delta AIC  
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(delta_aic)
 
  # Sort by lowest p value - have the strongest interaction effect, but maybe not best AIC
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(p_val_int)
    
  # Sort by highest deviance explained
    sum_df %>%
      dplyr::group_by(climate_var) %>%
      dplyr::arrange(desc(dev_explained))
    
    summary(sum_df$dev_explained)
    
  # Sort by Rsquared  - base model is 0.5620117
    sum_df %>%
      dplyr::arrange(desc(r_sq_adj))
    
    
  
## Calculating q value  
  # Vignette - https://bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
  
  hist(sum_df$p_val_int, nclass = 20)
  
  qvals <- qvalue(p = sum_df %>%
                    dplyr::filter(climate_var == "tmax_sum_dif") %>%
                    dplyr::pull(p_val_int),
                  fdr.level = 0.05)
  
  
  qvals$qvalues == sum_df$p_val_int
  
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
  
 plot(qobj_fdrlevel$qvalues ~ factor(as.numeric(qobj_fdrlevel$significant)))
 summary(qobj_fdrlevel$qvalues[qobj_fdrlevel$significant])
 summary(qobj_fdrlevel)
 
 
  summary(qobj_fdrlevel)
  
  # Sort by lowest q value  
  sum_df %>%
    dplyr::group_by(climate_var) %>%
    dplyr::arrange(q_val)
  
  

# Manhattan plots ---------------------------------------------------------


sum_df_sub <-  sum_df %>%
  dplyr::filter(climate_var == "tmax_sum_dif")    
    
  plot(-log10(sum_df_sub$q_val), pch = 20,
       xlab = "SNP", ylab = "-log10(P)", 
       main = "tmax_sum_dif", las = 1,
       col = factor(snp_pos$chrom %% 2), # To alternate colors for pseudo-chromosomes
       ylim = c(0, 2))
  abline(h = -log10(5e-8), lty = 2, col = "red") # Genome-wide significance line
  abline(h = -log10(1e-5)) # Suggestive significance line  
  abline(h = -log10(.05))
  
  
  

# RDA on top snps ---------------------------------------------------------

number = 100
  
top_snps <- sum_df %>%  
    dplyr::filter(climate_var == "tmax_sum_dif")  %>%
    dplyr::top_n(., -number, q_val) %>%
    dplyr::select(snp) %>%
    dplyr::mutate(snp = gsub("snp_", "", snp)) %>%
    dplyr::filter(snp %in% colnames(gen_dat_clim)) # TEMP FILTERING

# Based on delta AIC
top_snps <- sum_df %>%  
  dplyr::filter(climate_var == "tmax_sum_dif")  %>%
  dplyr::filter(aic_delta < -2) %>%
  dplyr::select(snp) %>%
  dplyr::mutate(snp = gsub("snp_", "", snp)) %>%
  dplyr::filter(snp %in% colnames(gen_dat_clim)) # TEMP FILTERING



bglr_gen_scaled[, top_snps$snp]


rda_full <- vegan::rda(bglr_gen_scaled[, top_snps$snp] ~ latitude + 
                         longitude + tmax_sum + tmin_winter + cwd + bioclim_04 + bioclim_15 + bioclim_18 + bioclim_19,
                       data = gen_dat_clim)



rda_full <- vegan::rda(bglr_gen_scaled[, top_snps$snp] ~ tmax_sum + Condition(latitude + longitude),
                       data = gen_dat_clim)

rda_full
summary(rda_full)

anova(rda_full)

anova(rda_full, by="axis", permutations = 100)

anova(rda_full, 
      by = "margin",
      permutations = 100)

round(summary(rda_full)$biplot[, c(1,2)], 2)

plot(rda_full)


plot(gen_dat_clim$tmax_sum, bglr_gen_scaled[, snp])
    cor.test(gen_dat_clim$tmax_sum, bglr_gen_scaled[, snp])

    summary(gam(bglr_gen_scaled[, snp] ~  s(tmax_sum), data = gen_dat_clim ))  
  
    

# Gradient forest test ----------------------------------------------------
    
    number = 50
    
    # Based on delta AIC
    top_snps = sum_df %>%
      dplyr::filter(climate_var == "tmax_sum_dif")  %>%
      dplyr::filter(converged == TRUE) %>%
      dplyr::filter(delta_aic > -1000) %>% 
      dplyr::filter(delta_aic < 500) %>%
      dplyr::top_n(., -number, delta_aic) %>%
      dplyr::arrange(delta_aic) %>%
      dplyr::select(snp) %>%
      dplyr::mutate(snp = gsub("snp_", "", snp))
    
   
  # Based on q value   
    top_snps <- sum_df %>%  
      dplyr::filter(climate_var == "tmax_sum_dif")  %>%
      dplyr::top_n(., -number, q_val) %>%
      dplyr::select(snp) %>%
      dplyr::mutate(snp = gsub("snp_", "", snp)) %>%
      dplyr::filter(snp %in% colnames(gen_dat_clim)) # TEMP FILTERING
    
    
    
    
    bglr_gen_scaled[, top_snps$snp]
    
    library(gradientForest)
    
    climate_vars_gf <- c(climate_vars, "latitude", "longitude")
    
    length(top_snps$snp)
    top_snps$snp
    
    gf <- gradientForest(cbind(bglr_gen_scaled, gen_dat_clim),
                         predictor.vars = climate_vars_gf, 
                        response.vars = top_snps$snp,
                      # response.vars = colnames(bglr_gen_scaled)[1:1000],
                         ntree = 500, transform = NULL, compact = T,
                         nbin = 201)
    
    gf$result # Rsquared of positive loci
    gf$species.pos.rsq
    
    mean(gf$result)
    
    plot(gf, plot.type = "O")
    
    most_important <- names(importance(gf))
    
    pred <- expand.grid(latitude = seq(35, 41, by = .1))
    pred$pred <- predict(gf, newdata = pred)
    
    pred
    
    plot(pred$latitude, pred$pred$latitude)
    
    library(randomForest)
    library(forestFloor)
    
    snp = "1_62307519"
    
    gen_dat_clim_sub <- gen_dat_clim %>%
                      dplyr::filter(!is.na(pull(., snp)))
    
    rfo <- randomForest(y = as.factor(pull(gen_dat_clim_sub, snp)),
                       x = gen_dat_clim_sub[, climate_vars_gf], 
                       keep.inbag = TRUE,  # mandatory,
                       importance = TRUE)  # recommended, else ordering by giniImpurity (unstable))
    rfo
    plot(rfo)
    
    rfo$importance
    
    ff = forestFloor(
      rf.fit = rfo,       # mandatory
      X = gen_dat_clim_sub[, climate_vars_gf],              # mandatory
      calc_np = FALSE,    # TRUE or FALSE both works, makes no difference
      binary_reg = FALSE  # takes no effect here when rfo$type="regression"
    )
    
    plot(ff)
    
    plot(ff,                       # forestFloor object
         plot_seq = 1:6,           # optional sequence of features to plot
         orderByImportance=TRUE    # if TRUE index sequence by importance, else by X column  
    )
    
    rf$importance
    
    pred <- expand.grid(latitude = seq(33, 42, by = .1),
                        tmax_sum = 0, tmin_winter = 0, random = 0, bioclim_04 = 0, PC1 = 0, PC2 = 0, PC3 = 0, longitude = seq(-123.5, -118.5, by = .1))
    pred$pred <- predict(rfo, newdata = pred)
    
    pred
    
    table(pred$pred)
    
    
    plot(pred$latitude, pred$pred)
    
    
    
    
    plot(gf, plot.type = "S", imp.vars = most_important,
          leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
          cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
          0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
    
    plot(gf, plot.type = "C", imp.vars = most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,
     0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, + 0.3, 0, 0)))
    
    plot(gf, plot.type = "C", imp.vars = most_important,
         show.species = F, common.scale = T, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
       + 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,  0.3, 0, 0)))
    
    plot(gf, plot.type = "P", show.names = T, horizontal = F,
         cex.axis = 1, cex.labels = 0.7, line = 2.5)
    
    
    

# Look at top models ------------------------------------------------------

  top_snp = "snp_4_93755225"
  
  plot(pull(dat_snp, top_snp), dat_snp$height_2017)
  plot(pull(dat_snp, tmax_sum_dif),pull(dat_snp, top_snp))
  
  summary(gam_tmax_sum_dif_snp_7_39668055)
  visreg(gam_tmax_sum_dif_snp_7_39668055, scale = "response")
  visreg2d(gam_tmax_sum_dif_snp_7_39668055, xvar = "tmax_sum_dif", 
           yvar = "snp_7_39668055", 
           scale = "response", plot.type = "persp", theta = 45, phi = 15)
  
  pairs.panels(dplyr::select(dat_snp, top_snp, latitude, longitude, climate_vars ))
  
