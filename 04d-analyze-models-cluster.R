
# Run model setup ---------------------------------------------------------

  source("./04b-setup-models-cluster.R")

# Read in cluster run output ----------------------------------------------

  # Path to results
  # Used in first submission to PNAS - "run_475547_tmax_sum_dif_training_set_resids_7030_2019_02_26"
  # path <- "run_475547_tmax_sum_dif_training_set_resids_7030_2019_02_26"

## 10 fold cross validation - 1st revision August 2019

   # Not across families - used for plotting main results in text
   path <- "run_642454_tmax_sum_dif_rgr_10fold_notacross_seed300"
   
   # Across families
   path <- "run_648224_tmax_sum_dif_rgr_10fold_acrossfams_seed300"
   
   
   ## Testing Tmax as main effect
   path <- "run_1041948_tmax_sum_dif_rgr_10fold_notacross_seed300_tmax_main"
  
   
   
  path_to_summaries <- paste0("./output/", path, "/model_summaries/")
  path_to_predictions <- paste0("./output/", path, "/model_predictions/")
  
  ## Load in data
  load(grep("*.Rdata", list.files(paste0("./output/", path), full = TRUE), value = TRUE))

  
 # Read in files
    sum_df_raw <- plyr::ldply(list.files(path_to_summaries, full = TRUE), read_csv)
    
    dim(sum_df_raw)
    sum_df_raw
    
    head(sum_df_raw)
    summary(sum_df_raw)
    str(sum_df_raw)
    
    which(table(sum_df_raw$snp) > 1) # Should return nothing
  
  
  ## Missing snps
  # Find if there indexes if needed to re-run manually or on the cluster
    which(!snp_col_names %in% sum_df_raw$snp)
    snp_col_names[which(!snp_col_names %in% sum_df_raw$snp)]
    length( which(!snp_col_names %in% sum_df_raw$snp))
    
    ## snp_chr5_22319422 doesn't have 2 levels
    
    ## Write snps that were run to file to rerun missing snps on cluster
      # write_csv(data.frame(snp = snp_col_names[which(snp_col_names %in% sum_df_raw$snp)]),
      #           path = "./output/skip_these_SNPs_2019_08_21_acrossfams.csv")

    
  # For backup
  sum_df <- sum_df_raw
  
  # How many acccession used for prediction
  table(sum_df_raw$acc_pred)
  
  
## Read in model predictions
  pred_df_raw <- plyr::ldply(list.files(path_to_predictions, full = TRUE), 
                             read_csv)
  pred_df_raw$genotype <- as.character(pred_df_raw$genotype)
  
  str(pred_df_raw)
  head(pred_df_raw)
  dim(pred_df_raw)
  summary(pred_df_raw)
  
  # Set missing data (NA) to 'missing'
  pred_df_raw$genotype[is.na(pred_df_raw$genotype)] <- "missing"
  
  # check that each SNP has the same number of predictions
  all(table(pred_df_raw$snp) == median(table(pred_df_raw$snp)))
  which(table(pred_df_raw$snp) != median(table(pred_df_raw$snp)))
  

  
# Convert to long format and set filters ---------------------------------------------------

  ## Convert to from wide to long format
  n_gen_long <- sum_df %>%
    dplyr::select(snp, n_gen0, n_gen1, n_gen2, n_missing, acc_pred, n) %>%
    gather(key = "genotype", value = "n_gen", -snp, -acc_pred, -n) %>%
    mutate(genotype = gsub("n_gen", "", genotype)) %>%
    mutate(genotype = gsub("n_", "", genotype)) # To reformat missing genotypes
  
  head(n_gen_long)
  
## Join all long dataframes together
  sum_df_long <- n_gen_long
  
  head(sum_df_long) 
  dim(sum_df_long)
  str(sum_df_long)
  
table(sum_df_long$genotype)

  
## Calculate height change in warmer conditions  
    pred_df_long =  pred_df_raw %>%
      dplyr::group_by(snp, genotype) %>%
      do(data.frame( 
     height_change_warmer_absolute = mean(.$pred[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold1 = mean(.$pred_fold1[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold2 = mean(.$pred_fold2[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold3 = mean(.$pred_fold3[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold4 = mean(.$pred_fold4[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold5 = mean(.$pred_fold5[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold6 = mean(.$pred_fold6[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold7 = mean(.$pred_fold7[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold8 = mean(.$pred_fold8[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold9 = mean(.$pred_fold9[.$tmax_sum_dif_unscaled > 0]),
     height_change_warmer_absolute_fold10 = mean(.$pred_fold10[.$tmax_sum_dif_unscaled > 0])
     ))
    
    head(pred_df_long)
    dim(pred_df_long)
    
    pred_df_long <- pred_df_long %>%
      ungroup()
    
    summary(pred_df_long$height_change_warmer_absolute)

    plot(pred_df_long$height_change_warmer_absolute, 
         pred_df_long$height_change_warmer_absolute_fold1)
    cor.test(pred_df_long$height_change_warmer_absolute, 
             pred_df_long$height_change_warmer_absolute_fold1)
  
 # Join prediction dataframe to summary dataframe
    sum_df_long_w_preds <- left_join(sum_df_long,
                                     pred_df_long)
    
    sum(is.na(sum_df_long_w_preds$height_change_warmer_absolute))
    
    table(sum_df_long_w_preds$genotype)
    
  
# Correct for multiple testing

  library(fdrtool)

  # F values of interaction term
  hist(sum_df$f_val_gen_int, breaks = 50)

  fdr_fvals = fdrtool(c(sum_df$p_val_gen_int),
                      statistic = "pvalue", plot = FALSE)

  # P values
  summary(fdr_fvals$pval)
  hist(fdr_fvals$pval, breaks = 40)
  sum(fdr_fvals$pval < 0.05)
  sort(fdr_fvals$pval)[1:50]
    
  # Q values
  hist(fdr_fvals$qval, breaks = 40) # Each q-value is the expected proportion of false positives among all SNP effects that are at least as extreme as that observed for the current SNP.
  sum(fdr_fvals$qval < 0.05) 
  summary(fdr_fvals$qval)
  sort(fdr_fvals$qval, decreasing = F)[1:100]
  
  
  
  # Make a panel plot for supplementary info with p values and q values
  library(qqman)
  
  par(mfrow = c(1,2))
  
  # P values
  hist(fdr_fvals$pval, breaks = 40, col = "steelblue2",
       xlab = "p-value", main = "Uncorrected p-values", las = 1)
  qq(fdr_fvals$pval, main = "Uncorrected p-values", las = 1)
  
  # Q values
  # hist(fdr_fvals$qval, breaks = 40, col = "steelblue2",
  #      xlab = "q-value", main = "q-values", las = 1)
  # qq(fdr_fvals$qval, main = "q-values", las = 1)
  # sort(fdr_fvals$qval, decreasing = T)
  
  # Plot figure
  #Need to hand-edit the axis labesl for q values to say Q instead of P
  # dev.copy(pdf, paste0("./figs_tables/Figure S9 - histogram and qqplots of p values ", Sys.Date(), ".pdf"),
  #          width = 7, height = 3.5, useDingbats = FALSE)
  # dev.off()

  par(mfrow = c(1,1))
  
  
  ## Combine p values of main and interaction terms
  # library(metap)
  # 
  # for(row in 1:nrow(sum_df)){
  # 
  #   if(row %% 500 == 0){
  #     cat("Working on row... ", row, " \n")
  #   }
  # 
  #   # For full dataset
  #   sum_df$p_val_gen_combined[row] <- metap::sumlog(c(sum_df$p_val_gen_int[row],
  #                                                sum_df$p_val_gen_main[row]))$p
  #   # Loop through folds
  #   for(fold in 1:10){
  #     pval_main <- sum_df[row, paste0("p_val_gen_main_fold", fold)]
  #     pval_int <- sum_df[row, paste0("p_val_gen_int_fold", fold)]
  # 
  #     # Check to see if pval of main effect is NaN
  #     if(is.na(pval_main)){
  #       sum_df[row, paste0("p_val_gen_combined_fold", fold)] <- pval_int
  #       next
  #     }
  # 
  #     # Assign new p value for fold
  #     sum_df[row, paste0("p_val_gen_combined_fold", fold)] <-
  #           metap::sumlog(c(pval_main, pval_int))$p
  #   }
  # }
  # 
  # hist(sum_df$p_val_gen_combined)
  # summary(sum_df$p_val_gen_combined)
  
  
# Join back to main dataframe
  sum_df$q_val <- fdr_fvals$qval
  sum_df_long_w_preds <- dplyr::left_join(sum_df_long_w_preds,
                                          dplyr::select(sum_df, snp, q_val),
                                          by = "snp")
  summary(sum_df_long_w_preds$q_val)
  
  

# Figure 2 - Manhattan plot of Q values ---------------------------------------------------------
  
   man_df1 <- data.frame(SNP = snp_col_names, 
                         snp_pos)
   
   man_df1$index <- 1:nrow(man_df1) ## For x axis plotting
   
   man_df <- man_df1
   
   man_df <- left_join(man_df1, sum_df, by = c("SNP" = "snp"))
  
   head(man_df)
   dim(man_df)
   
   ## Clean up chromosome names - clear scaffold names
   man_df <- man_df %>%
     mutate(chrom = ifelse(grepl("chr", chrom), chrom, ""))
  
   # For labels on x axis 
   x_axis_breaks <- man_df %>%
     dplyr::select(index, chrom) %>%
     group_by(chrom) %>%
     summarise_all(mean)
   
   # For vertical lines showing chromosome sections
   chrom_breaks <- man_df %>%
     dplyr::select(index, chrom) %>%
     group_by(chrom) %>%
     summarise_all(max)
   
   # Set colors for chromosomes
   # Discrete color cateogies - https://twitter.com/debivort/status/994583058031546369
   odd_chroms <- paste0("chr", seq(1, 12, 2))
   even_chroms <- paste0("chr", seq(2, 12, 2))
   
   man_df$color <- NA
   man_df$color[man_df$chrom %in% odd_chroms] <- "#1f78b4"
   man_df$color[man_df$chrom %in% even_chroms] <- "#a6cee3"
   
   # Start plot
   ggplot(man_df, aes(x = index, 
                    y = -log10(p_val_gen_int),
                    label = SNP)) +
  
     # Vertical lines
     geom_vline(xintercept = chrom_breaks$index, col = "grey50", alpha = 0.8, size = .2) + 
     
     # Horizontal line
   #  geom_hline(yintercept = -log10(q_val_thresh), lty = 1, size = .25) +

     # Add points
     geom_point(pch = 16, 
                col = man_df$color,
                alpha = 1,
                size = 1.75) +
     
     # # Overplot top and bottom SNPS
     # geom_point(data = man_df[man_df$top_snp | man_df$bottom_snp, ],
     #            aes(x = index,
     #                y = -log10(q_val)),
     #                bg = man_df[man_df$top_snp | man_df$bottom_snp, ]$color,
     #                col = "black",
     #                pch = 21,
     #                size = .65 + 1) +  
     
     # Theme options
     scale_x_continuous(expand = c(0.00, 0.00), 
                        breaks = x_axis_breaks$index,
                        labels = x_axis_breaks$chrom) + 
     scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) + 
     ylab("-log10(p value)") + 
     xlab("Chromosome") + 
     theme_bw(8) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
          # legend.position = "none") +
             )+
     # axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
     NULL
   
   # # Save as file - PNG
     # ggsave(paste0("./figs_tables/Figure S9 - manhattan plot interaction term_",
     #                      Sys.Date(), ".png"),
     #          width = 19, height = 4, units = "cm", dpi = 600)


   
   
   
# Figure 2 - Manhattan plot of growth predictions  --------------------------------------------
   
   man_df1 <- data.frame(SNP = snp_col_names, 
                         snp_pos)
   
   man_df1$index <- 1:nrow(man_df1) ## For x axis plotting
   
   man_df <- left_join(man_df1, sum_df_long_w_preds, by = c("SNP" = "snp"))
   
   head(man_df)
   
   # Set alpha for q values
  # man_df$q_val_alpha <- ifelse(man_df$q_val < q_val_thresh, 1, .35)
   man_df$q_val_alpha <- 1
   
   # Set up shapes
   man_df <- man_df %>%
     mutate(pch = ifelse(genotype == 0 | genotype == 2, 16, 16))
  
   table(man_df$pch)
   
   ## Clean up chromosome names - clear scaffold names
   man_df <- man_df %>%
     mutate(chrom = ifelse(grepl("chr", chrom), chrom, ""))
   
   # For labels on x axis 
   x_axis_breaks <- man_df %>%
     dplyr::select(index, chrom) %>%
     group_by(chrom) %>%
     summarise_all(mean)
   
   # For vertical lines showing chromosome sections
   chrom_breaks <- man_df %>%
     dplyr::select(index, chrom) %>%
     group_by(chrom) %>%
     summarise_all(max)
   
   # Start plot
   ggplot(man_df, aes(x = index, 
                      y = height_change_warmer_absolute,
                      col = height_change_warmer_absolute,
                      label = SNP)) +
     
     # Vertical lines
     geom_vline(xintercept = chrom_breaks$index, col = "grey50", alpha = 0.8, size = .2) + 
     
    # Add points
     geom_point(pch = man_df$pch, 
                alpha = 1, 
                size = 1.25) +
     
     # Horizontal line
     # geom_hline(yintercept = mean(man_df$height_change_warmer_absolute, na.rm = TRUE), 
     #            lty = 2, size = .25) +
     
     geom_hline(yintercept = 0, 
                lty = 2, size = .25) +
     
    
   # Scale colors
   scale_colour_distiller(type = "seq", palette = "PRGn", 
                          direction = 1, 
                          aesthetics = c("color", "fill"),
                          name = "% Change in RGR \n with 4.8°C increase") + 
     
     # Theme options
     scale_x_continuous(expand = c(0.00, 0.00), 
                        breaks = x_axis_breaks$index,
                        labels = x_axis_breaks$chrom) + 
   #  scale_y_continuous(breaks = seq(-30, 20, 5)) + 
     ylab("-log10(q-value)") + ## Need same y axis as interaction plot so dimensions are exactly the same
     xlab("Chromosome") + 
     theme_bw(8) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           legend.position = "none") +
     #)+
     # axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
     NULL
   
   # # Save as file - PNG
     # ggsave(paste0("./figs_tables/fig2/Figure 2 - manhattan plot growth predictions",
     #                      Sys.Date(), ".png"),
     #          width = 19, height = 4, units = "cm", dpi = 600)


# Calculate polygenic risk score ------------------------------------------
   
   
   # Save workspace to upload to hoffman
  # save.image(file = "./output/gebv_calculation_acrossfamily.Rdata")
   
   
   ## Parallelize calculations across main calculations and folds
   cores=detectCores()
   cl <- makeCluster(cores[1]-1) #not to overload your computer
   registerDoParallel(cl)
   
   # Set up model terms
   fixed_effects_resids <- paste0(paste0("rgr_resids ~ s(gebv_by_count, bs=\"cr\") + 
                                               s(tmax_sum_dif, by = gebv_by_count, bs=\"cr\")"))
   
   fixed_effects_resids_permuted <- paste0(paste0("rgr_resids_permuted ~ s(gebv_by_count, bs=\"cr\") + 
                                                 s(tmax_sum_dif, by = gebv_by_count, bs=\"cr\")"))
   
   
   # Calculate polygenic risk score for each maternal genotype
   gen_dat_moms <- dplyr::select(gen_dat_clim, id, accession, snp_col_names)
   dim(gen_dat_moms)
   

  # Define function that calculates GEBV 
   calc_gebv <- function(){
     
     # Choose SNPs
     n_snps_values <- c(25, 50, seq(100, 999, by = 100), seq(1000, 12000, by = 500), nrow(sum_df))

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
         # for(perm in 1:500){
         #   
         #   if(perm %% 10 == 0){
         #     cat("Working on permutation: ", perm, "... \n")
         #   }
         #   
         #   # Permute phenotype
         #   dat_snp_testing$rgr_resids_permuted <- sample(dat_snp_testing$rgr_resids)
         #   
         #   # Gam with interaction effect
         #   gam_test_resids_permuted = bam(formula = formula(fixed_effects_resids_permuted),
         #                         data = dat_snp_testing,
         #                         discrete = FALSE, 
         #                         nthreads = 8,
         #                         method = "fREML")
         # 
         # gebv_df_summary[gebv_df_summary$n_snps == x, paste0("r2_permuted_", perm)] <- summary(gam_test_resids_permuted)$r.sq
         # 
         # }
         
         ## Make prediction of effect size
           # newdata <-  expand.grid(section_block = "Block1_1",
           #                         locality = "FHL",
           #                         height_2014 = 0,
           #                         accession = "6",
           #                         tmax_sum_dif = c(forward_transform(c(0, 4.8), "tmax_sum_dif",
           #                                                            scaled_var_means_gbs_all,
           #                                                            scaled_var_sds_gbs_all)),
           #                         PC1_gen = 0, PC2_gen = 0,
           #                         PC1_clim = 0, PC2_clim = 0,
           #                         gebv_by_count = c(-1, 0, 1)) # SD away from average
           # 
           # # Make predictions
           # newdata$pred <- predict(gam_test_resids,
           #                         newdata = newdata,
           #                         se.fit = FALSE, type = "response")
           # 
           # # Add on baseline predictions from model without genetic information
           # newdata$pred <- newdata$pred + predict(gam_snp_all,
           #                                        newdata = newdata,
           #                                        se.fit = TRUE, type = "response")$fit
           # 
           # ## Calculating differences in growth rates
           # growth_dif =  newdata %>%
           #   dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
           #                                                        means = scaled_var_means_gbs_all,
           #                                                        sds = scaled_var_sds_gbs_all))
           # 
           # # Comparing high to mid
           # high_mid <- (growth_dif$pred[growth_dif$gebv_by_count == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
           #     growth_dif$pred[growth_dif$gebv_by_count == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
           #   abs(growth_dif$pred[growth_dif$gebv_by_count == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100
           # 
           # # Compare high to low
           # high_low <-  (growth_dif$pred[growth_dif$gebv_by_count == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
           #     growth_dif$pred[growth_dif$gebv_by_count == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
           #   abs(growth_dif$pred[growth_dif$gebv_by_count == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100
           # 
           # gebv_df_summary$high_mid[gebv_df_summary$n_snps == x] <- high_mid
           # gebv_df_summary$high_low[gebv_df_summary$n_snps == x] <- high_low
         
       } # End model loop
       
       x = x + 1 # Move onto next SNP
       
     } # End snp loop
     
     return(list(gebv_df_summary, gebv_values)) # Return summary df as output of loop
     
   } # End calc_gebv function
   
   ## Run parallelized funtion to calculate gebv for each fold and 
   gebv_df_summary_list <- foreach(fold=0:10, 
                   .packages = c("tidyverse", "mgcv"),
                   .combine = rbind) %dopar% {
     calc_gebv() # call function
                   }
   
   beepr::beep(4)
  
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
   
   gebv_df_summary %>%
     arrange(desc(r2))
   
   gebv_df_summary %>%
     group_by(fold) %>%
     summarise_all(max)
   
   gebv_df_summary %>%
     group_by(fold) %>%
     summarise_all(max) %>%
     filter(fold != 0) %>%
     ungroup() %>%
     summarise_all(mean)
   
   gebv_df_summary %>%
     filter(fold != 0) %>%
     group_by(fold) %>%
     filter(r2 == max(r2)) %>%
     ungroup() %>%
     summarise_all(mean)
   
   # How many snps at max r2?
   gebv_df_summary %>%
     group_by(fold) %>%
     filter(r2 == max(r2))
  
   
   ## Plot relationship between # of SNPS used and variance explained in testing set
   ggplot(gebv_df_summary, aes(x = n_snps, y = r2, col = factor(fold))) + 
     geom_point(size = 2) + 
     geom_line() + 
     xlab("Number of SNPs in GEBV estimation") +
     #ylab("R2 adjusted in testing set") +
     ylab(bquote('R'^2[adj]~'training test')) +
     scale_x_continuous(breaks = gebv_df_summary$n_snps) +
     # scale_y_continuous(breaks = seq(0, 0.02, by = 0.0025), limits = c(-0.0015, 0.015)) + 
     theme_bw(7) + 
     geom_hline(yintercept = 0, lty = 2) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           plot.margin = (margin(0, 0, 0, 0, "cm")),
           axis.text.x = element_text(angle = 60, hjust = 1)) +
     NULL
   
   
   
   ## Plot relationship between # of SNPS used and variance explained in full model
   gebv_df_fold0 <-  gebv_df_summary %>%
     filter(fold == 0)
   
   
   ggplot(gebv_df_fold0, aes(x = n_snps, y = r2)) + 
     geom_point(size = 1) + 
     geom_line() + 
     xlab("Number of SNPs in GEBV estimation") +
     #ylab("R2 adjusted in testing set") +
     ylab(bquote('R'^2[adj])) +
     scale_x_continuous(breaks = gebv_df_fold0$n_snps[c(-1, -3, -5, -7, -9, -11)]) +
     # scale_y_continuous(breaks = seq(0, 0.02, by = 0.0025), limits = c(-0.0015, 0.015)) + 
     theme_bw(5) + 
     geom_hline(yintercept = 0, lty = 2) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           plot.margin = (margin(0, 0, 0, 0, "cm")),
           axis.text.x = element_text(angle = 60, hjust = 1)) +
     NULL
   
   
   
   ### Read in GEBV cross-validation results
   cv_files <- list.files("./output/cv_permutations_across", full.names = T)
   
   cv_out = plyr::ldply(cv_files, read_csv)
   
   
   # Calculate average R2 and 95% confidence interval
   cv_out2 <- cv_out %>%
    mutate(r2_permuted = rowMeans(dplyr::select(cv_out, contains("r2_permuted")))) %>%
     dplyr::select(n_snps, fold, r2, r2_permuted)

   cv_out2 %>%
     filter(fold != 0) %>%
     group_by(fold) %>%
     filter(r2 == max(r2)) %>%
     ungroup() %>%
     summarise_all(mean)
   
   
   cv_out2 %>%
     filter(fold != 0) %>%
     group_by(fold) %>%
     filter(r2 == max(r2)) %>%
     ungroup() %>%
     summarise_all(sd)
   
   
   # Save as file
   # ggsave(paste0("./figs_tables/Figure S4 - testing variance explained by number of snps",
   #                      Sys.Date(), ".png"),
   #          width = 12, height = 4, units = "cm", dpi = 300)
   
   
   
   
   ## Get GEBVs
   gebvs <- bind_tibbles(gebv_df_summary_list[, 2], mode = "cbind")
   head(gebvs)

   
  # Set top gebv
  top_r2_index <- which(gebv_df_summary$r2[gebv_df_summary$fold == 0] == 
                          max(gebv_df_summary$r2[gebv_df_summary$fold == 0],
                                                          na.rm = T))
  top_r2_n_snps <- gebv_df_summary$n_snps[gebv_df_summary$fold == 0][top_r2_index]
  
  cat("Top number of SNPs is: ", top_r2_n_snps, " ... \n")
  
  # Set top gebv
  gen_dat_moms$gebv_best <- gebvs[ , paste0("gebv_fold0_", top_r2_n_snps)]
 
  # Scale gebv scores
  gebv_best_mean <- mean(gen_dat_moms$gebv_best)
  gebv_best_sd <- sd(gen_dat_moms$gebv_best)
  gen_dat_moms$gebv_best_unscaled <- gen_dat_moms$gebv_best
  gen_dat_moms$gebv_best_scaled <- (gen_dat_moms$gebv_best - gebv_best_mean) / gebv_best_sd
  
  summary(gen_dat_moms$gebv_best_scaled); hist(gen_dat_moms$gebv_best_scaled, breaks = 50)
  
  
# Plot predictions of PRS -------------------------------------------------

    
    ## Calculate number of top and bottom SNPs 
    dat_snp_gebv <- dat_snp_all %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
      dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
                    rgr, rgr_resids, 
                    section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
      dplyr::left_join(., dplyr::select(gen_dat_moms, accession, gebv_best, gebv_best_scaled)) %>%
      dplyr::mutate(accession = factor(accession))
    
    head(dat_snp_gebv)
    dim(dat_snp_gebv)
    
    summary(dat_snp_gebv$gebv_best)
    summary(dat_snp_gebv$gebv_best_scaled)
    
    dat_snp_gebv_unscaled <- dat_snp_gebv
    
    
    # Set formula
     fixed_effects <- paste0(paste0("rgr_resids ~ s(gebv_best_scaled, bs=\"cr\") + 
                                           s(tmax_sum_dif, by = gebv_best_scaled, bs=\"cr\")"))
 
    # Gam with interaction effect
    gam_gebv_all = bam(formula = formula(fixed_effects),
                  data = dat_snp_gebv,
                  discrete = FALSE, 
                  nthreads = 8,
                  method = "fREML")
    
    summary(gam_gebv_all)
    
    ## Save output
    
   # Save gam summary to file
   #  sink(file = paste0("./figs_tables/Table S2 - gebv_all_gam_model_summary_",
   #                     Sys.Date(), ".txt"))
   #  summary(gam_gebv_all)
   # # anova(gam_gebv_all)
   #  sink()
    
    
    
## Plot predictions

    # Set up dataframe for prediction
    newdata <-  expand.grid(section_block = "Block1_1",
                            locality = "FHL",
                            height_2014 = 0,
                            accession = "6",
                            tmax_sum_dif = c(seq(min(dat_snp_all_unscaled$tmax_sum_dif),
                                                 max(dat_snp_all_unscaled$tmax_sum_dif),
                                                 length.out = 150), 
                                             forward_transform(4.8, "tmax_sum_dif",
                                                               scaled_var_means_gbs_all,
                                                               scaled_var_sds_gbs_all)),
                            PC1_gen = 0, PC2_gen = 0,
                            PC1_clim = 0, PC2_clim = 0,
                            gebv_best_scaled = c(-1, 0, 1)) # SD away from average
    
    # Make predictions
    newdata$pred <- predict(gam_gebv_all,
                            newdata = newdata,
                            se.fit = TRUE, type = "response")$fit
    
    newdata$se <- predict(gam_gebv_all,
                          newdata = newdata,
                          se.fit = TRUE, type = "response")$se
    
    
    # Add on baseline predictions from model without genetic information
    newdata$pred <- newdata$pred + predict(gam_snp_all,
                                            newdata = newdata,
                                            se.fit = TRUE, type = "response")$fit

    # Plot predictions
    newdata %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all)) %>%
      ggplot(., aes(x = tmax_sum_dif_unscaled, y = pred)) +
      geom_vline(aes(xintercept = 0), lty = 2, size = .25) +
      geom_vline(xintercept = 4.8, size = .25) +
      geom_vline(xintercept = -5.2, size = .25) +
      
      
      geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96 * se, 
                      fill = factor(gebv_best_scaled)), alpha = .35) +
      geom_line(aes(col = factor(gebv_best_scaled))) + 
      scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
    #  scale_y_continuous(breaks = seq(0, 1.75, by = 0.25)) + 
      ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
      xlab("Tmax transfer distance (°C)") +
      scale_color_manual(values = c("#af8dc3", "#227eb8", "#7fbf7b")) +
      scale_fill_manual(values =  c("#af8dc3", "#227eb8", "#7fbf7b")) +
      
      theme_bw(8) + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      NULL
    
    
    # Save plot output
    # ggsave(filename = paste0("./figs_tables/fig2/Figure 2 - gebv growth ", Sys.Date(), ".pdf"),
    #        units = "cm", width = 8, height = 6)

    
    
    ## Calculating differences in growth rates
    growth_dif =  newdata %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                         means = scaled_var_means_gbs_all,
                                                         sds = scaled_var_sds_gbs_all))

    # Comparing high to mid
    (growth_dif$pred[growth_dif$gebv_best_scaled == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
        growth_dif$pred[growth_dif$gebv_best_scaled == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
      abs(growth_dif$pred[growth_dif$gebv_best_scaled == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100

    # Compare high to low
    (growth_dif$pred[growth_dif$gebv_best_scaled == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
        growth_dif$pred[growth_dif$gebv_best_scaled == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
     abs(growth_dif$pred[growth_dif$gebv_best_scaled == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100

    
    
    

# Calculate PRS based on BLUP ---------------------------------------------

    ## Quadratic model    
    gam_acc_quad <- lmer(rgr ~
                           height_2014 + I(height_2014^2) + I(tmax_sum_dif^2) + 
                           section_block  + PC1_gen + PC2_gen + PC1_clim + PC2_clim +  
                           (tmax_sum_dif | accession),
                         data = dat_snp_all)
    
    summary(gam_acc_quad)
    
    
    ## Calculate RGR at tmax > 0 based on individual accession curves
    # Make dataframe to save 'true' GEBV
    acc_rgr_absolute <- data.frame(accession = NA, 
                                   height_change_warmer_quad = NA) 
    
    plot = FALSE
    x = 1
    for(acc in levels(dat_snp_all$accession)){
      
      cat("Working on accession:", acc, "... \n")
      
      # Build prediction dataframe
      pred_acc <-  expand.grid(height_2014 = 0,
                               accession = acc,
                               locality = "FHL",
                               section_block = "Block1_1",
                               PC1_clim =0, PC2_clim = 0,
                               PC1_gen = 0, PC2_gen = 0)
      
      pred_acc <- left_join(pred_acc, expand.grid(accession = acc,
                                            tmax_sum_dif= c(seq(min(dat_snp_all$tmax_sum_dif),
                                                                  max(dat_snp_all$tmax_sum_dif),
                                                                  length.out = 100))))
      
      pred_acc$tmax_sum_dif_unscaled = back_transform(pred_acc$tmax_sum_dif,
                                                      var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all)
      
      pred_acc$preds_quad <- predict(gam_acc_quad, newdata = pred_acc)
      
      ## Scaled plot
        # plot(pred_acc$tmax_sum_dif_unscaled,
        #      pred_acc$preds,
        #      type  = "p", lwd = 1.5, las = 1,
        #      ylab = "Scaled RGR effect", xlab = "Tmax dif",
        #      main = paste0("Accession: ", acc))
        # abline(v = 0, lty = 2)
        # 
        # sub_progeny <- dat_snp_all[dat_snp_all$accession == acc, ]
        # sub_progeny$tmax_sum_dif_unscaled = back_transform( sub_progeny$tmax_sum_dif,
        #                                                 var = "tmax_sum_dif",
        #                                                 means = scaled_var_means_gbs_all,
        #                                                 sds = scaled_var_sds_gbs_all)
        # sub_progeny$y <- mean(pred_acc$preds, na.rm = T)
        # points(sub_progeny$tmax_sum_dif_unscaled,
        #        sub_progeny$y,
        #        pch = 19, col = "blue")
      
      acc_rgr_absolute[x, "accession"] <- acc
      acc_rgr_absolute[x, "height_change_warmer_quad"] <- mean(pred_acc$preds_quad[pred_acc$tmax_sum_dif_unscaled > 0])
      
      x = x + 1
    }
    
    acc_rgr_absolute
    summary(acc_rgr_absolute)
    hist(acc_rgr_absolute$height_change_warmer_quad)
    
    # Comapre with estimates of PRS from gams
    acc_rgr_absolute$accession <- as.numeric(as.character(acc_rgr_absolute$accession))
    acc_rgr_absolute2 <- left_join(acc_rgr_absolute, 
                             dplyr::select(gen_dat_moms, accession, gebv_best_scaled))
    
    plot(acc_rgr_absolute2$gebv_best_scaled, scale(acc_rgr_absolute2$height_change_warmer_quad))
    cor.test(acc_rgr_absolute2$gebv_best_scaled, acc_rgr_absolute2$height_change_warmer_quad)
    
    
    ### Calculate BLUPs based on BGLR
    library(BGLR)
    
    gen_dat_bglr <- gen_dat_moms %>%
      filter(accession %in% acc_rgr_absolute$accession)
    
    dim(gen_dat_bglr)
    
    
    # Order properly
    gen_dat_bglr <-  gen_dat_bglr[match(acc_rgr_absolute$accession, gen_dat_bglr$accession), ]

    all(gen_dat_bglr$accession == acc_rgr_absolute$accession)
    
    
    # Computing the genomic relationship matrix for SNPs
    A.snp.temp <- scale(as.matrix(gen_dat_bglr[, snp_col_names]),center=TRUE,scale=TRUE)
    
    # Impute missing data
    A.snp.temp[is.na(A.snp.temp)] <- 0
    
    A.snp <- tcrossprod(A.snp.temp)/ncol(A.snp.temp)
    
    
    # Fitting the model
    
    # Parameters    
    nIter <- 50000
    burnIn <- 20000
    verbose <- FALSE
    
    # Set response variable - estimated from quadratic model
    y = acc_rgr_absolute$height_change_warmer_quad
    summary(y)
    

    ## Model with SNP data
    mod.snp <- BGLR(y = y,
                    ETA = list(list(K = A.snp, model = 'RKHS')),
                    nIter = nIter, burnIn = burnIn,
                    saveAt='./output/bglr/', 
                    verbose = verbose)

    # Variance explained
    1 - mod.snp$varE/var(y)
    
    
    ## Compare estimated breeding values to true breeding values
    bglr_gebv <- mod.snp$ETA[[1]]$u
    acc_rgr_absolute$bglr_gebv <- c(scale(bglr_gebv))

    
    ## Plot for Supplementary information
      png("./figs_tables/Figure S10 - GEBV correlation gBLUP and GAMs.png", 
          width = 6, height = 5, units = "in", res = 300)
      
      plot(scale(acc_rgr_absolute2$gebv_best_scaled), scale(bglr_gebv), pch = 19,
           xlab = "GEBVs estimated by SNP-by-SNP GAM", ylab = "GEBVs Estimated by gBLUP",
           col = "black")
      abline(a = 0, b = 1, lwd = 2)
      
      cor_overall <- cor(acc_rgr_absolute2$gebv_best_scaled, bglr_gebv)
      
      legend("bottomright", legend = c(paste0("R = ", round(cor_overall, 3))),
             col = c("black"), pch = 19)
      
      dev.off()     
      
    
    
    ## Calculate number of top and bottom SNPs 
    dat_snp_gebv_blup <- dat_snp_all %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
      dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
                    rgr, rgr_resids, 
                    section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
      dplyr::left_join(., dplyr::select(gen_dat_moms, accession, gebv_best, gebv_best_scaled))
    
    head(dat_snp_gebv_blup)
    dim(dat_snp_gebv_blup)

    ## Add GEBV from BLUP
    dat_snp_gebv_blup <- left_join(dat_snp_gebv_blup, acc_rgr_absolute)
    head(dat_snp_gebv_blup)

    dat_snp_gebv_blup_unscaled <- dat_snp_gebv_blup
    
    
    # Set formula
    fixed_effects <- paste0(paste0("rgr_resids ~ s(bglr_gebv, bs=\"cr\") + 
                                   s(tmax_sum_dif, by = bglr_gebv, bs=\"cr\")"))
    
    # Gam with interaction effect
    gam_gebv_all = bam(formula = formula(fixed_effects),
                      data = dat_snp_gebv_blup,
                      discrete = TRUE, 
                      nthreads = 8,
                      method = "fREML",
                      #   family = "tw",
                      control = list(trace = FALSE))
    
    summary(gam_gebv_all)
    
    ## Save output
    
    # Save gam summary to file
    #  sink(file = paste0("./figs_tables/Table S2 - gebv_all_gam_model_summary_",
    #                     Sys.Date(), ".txt"))
    #  summary(gam_gebv_all)
    # # anova(gam_gebv_all)
    #  sink()
    
    

    ## Plot predictions
    
    # Set up dataframe for prediction
    newdata <-  expand.grid(section_block = "Block1_1",
                            locality = "FHL",
                            height_2014 = 0,
                            accession = "6",
                            tmax_sum_dif = c(seq(min(dat_snp_all_unscaled$tmax_sum_dif),
                                                 max(dat_snp_all_unscaled$tmax_sum_dif),
                                                 length.out = 150), 
                                             forward_transform(4.8, "tmax_sum_dif",
                                                               scaled_var_means_gbs_all,
                                                               scaled_var_sds_gbs_all)),
                            PC1_gen = 0, PC2_gen = 0,
                            PC1_clim = 0, PC2_clim = 0,
                            bglr_gebv = c(-1, 0, 1)) # SD away from average
    
    # Make predictions
    newdata$pred <- predict(gam_gebv_all,
                            newdata = newdata,
                            se.fit = TRUE, type = "response")$fit
    
    newdata$se <- predict(gam_gebv_all,
                          newdata = newdata,
                          se.fit = TRUE, type = "response")$se
    
    
    # Add on baseline predictions from model without genetic information
    newdata$pred <- newdata$pred + predict(gam_snp_all,
                                           newdata = newdata,
                                           se.fit = TRUE, type = "response")$fit
    
    # Plot predictions
    newdata %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all)) %>%
      ggplot(., aes(x = tmax_sum_dif_unscaled, y = pred)) +
      geom_vline(aes(xintercept = 0), lty = 2, size = .25) +
      geom_vline(xintercept = 4.8, size = .25) +
    #  geom_vline(xintercept = -5.2, size = .25) +
      
      
      geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96 * se, 
                      fill = factor(bglr_gebv)), alpha = .35) +
      geom_line(aes(col = factor(bglr_gebv))) + 
      scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
     # scale_y_continuous(breaks = seq(0, 1.75, by = 0.25)) + 
      ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
      xlab("Tmax transfer distance (°C)") +
      scale_color_manual(values = c("#af8dc3", "#227eb8", "#7fbf7b")) +
      scale_fill_manual(values =  c("#af8dc3", "#227eb8", "#7fbf7b")) +
      
      theme_bw(8) + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      NULL
    
    
    
    # Save plot output
    ggsave(filename = paste0("./figs_tables/Figure S10 - gebv growth from gBLUP ", Sys.Date(), ".pdf"),
           units = "cm", width = 8, height = 6)

    
    
    
    
    
    
    
    ## Calculate BLUPS based on GAPIT - http://zzlab.net/GAPIT/gapit_help_document.pdf - https://github.com/jiabowang/GAPIT3#gs
    
    # library(multtest)
    # library(gplots)
    # library(LDheatmap)
    # library(genetics)
    # library(ape)
    # library(EMMREML)
    # library(compiler) #this library is already installed in R
    # library("scatterplot3d")
    # source("http://zzlab.net/GAPIT/gapit_functions.txt")
    # 
    # 
    # 
    # ## Output in GEMMA format
    # ### Output SNP information in BIMBAM format
    # snp.sim <- as.data.frame(X.snp.unscaled) %>%
    #   `colnames<-`(paste0("snp_",
    #                       str_pad(as.character(1:ncol(X.snp.unscaled)),
    #                               4, pad = "0"))) %>%
    #   mutate(ID = paste0("ind_",
    #                      str_pad(as.character(1:nrow(X.snp.unscaled)), 4, pad = "0"))) %>%
    #   gather(key = "site", value = "value", -ID) %>%
    #   spread(key = ID, value = value) %>%
    #   mutate(allele_type1 = "C", allele_type2 = "T") %>%
    #   select(site, allele_type1, allele_type2, everything())
    # 
    # head(snp.sim[, 1:20])
    # 
    # dim(snp.sim)
    # 
    # head(X.snp.unscaled[, 1:10])
    # cbind()
    # 
    # # Write to file
    # # write_tsv(snp.sim,
    # #           "./data_simulation/snp.sim.txt",
    # #           col_names = FALSE)
    # 
    # ## Make SNP annotation file
    # snp.sim.annotation <- data.frame(snp.name = snp.sim$site,
    #                                  snp.position = str_split(snp.sim$site,
    #                                                           "_", simplify = T)[,2],
    #                                  snp.chrom = "1")
    # 
    # head(snp.sim.annotation)
    # 
    # myGM <- data.frame(snp = snp_col_names,
    #                    chrom = "1",
    #                    position = 1:length(snp_col_names))
    # head(myGM)
    # 
    # # write_tsv(x = snp.sim.annotation,
    # #           path = "./data_simulation/snp.sim.annotation.txt",
    # #           col_names = FALSE)
    # 
    # myGD <- cbind(gen_dat_moms)
    # phenotype <- 
    # 
    # myGAPIT_gBLUP <- GAPIT(
    #   Y=acc_rgr_absolute[, c("accession", "height_change_warmer")],
    #   GD=gen_dat_moms[, c("accession", snp_col_names)],
    #   GM=myGM,
    #   model="gBLUP",
    #   PCA.total=1,
    #   file.output=T,
    #   SNP.MAF = 0.05
    # )
    # 
    # ## Read in gblup results
    # 
    # gapit_out <- read_csv("./GAPIT.MLM.Pred.result.csv")
    # 
    # # Put in correct order
    #  gapit_out <-  gapit_out[match(acc_prs$accession, gapit_out$Taxa),]
    #  
    #  all(gapit_out$Taxa == acc_prs$accession) # Should return TRUE
    # 
    # 
    # 
    # 
    #  png(paste0("./figures/BLUP GAPIT correlation with true GEBV ", 
    #            n_sites, "_sites.png"), width = 8, height = 5, units = "in", res = 300)
    # 
    # plot(scale(acc_prs$prs_true), scale(gapit_out$BLUP), pch = 19,
    #      xlab = "GEBV True", ylab = "GEBV Estimated (BLUP)", col = col)
    # abline(a = 0, b = 1, lwd = 2)
    # 
    # cor_overall <- cor(acc_prs$prs_true, scale(gapit_out$BLUP))
    # cor_training <- cor(acc_prs$prs_true[acc_prs$accession %in% training_moms], 
    #                     scale(gapit_out$BLUP)[gen_dat_moms$accession %in% training_moms])
    # cor_testing <- cor(acc_prs$prs_true[!acc_prs$accession %in% training_moms], 
    #                    scale(gapit_out$BLUP)[!gen_dat_moms$accession %in% training_moms])
    # 
    # 
    # legend("bottomright", legend = c(paste0("R = ", round(cor_overall, 3)),
    #                                  paste0("R = ", round(cor_training, 3)),
    #                                  paste0("R = ", round(cor_testing, 3))),
    #        col = c("black", "steelblue", "coral"), pch = 19)
    # 
    # dev.off()     
    # 
    
    
    ## Compare BLUP approach to GAM approach
    # png(paste0("./figures/GAM and BLUP correlation with true GEBV ", 
    #            n_sites, "_sites.png"), width = 10, height = 10, units = "in", res = 300)
    # pairs.panels(data.frame(gebv_true = acc_prs$prs_true, 
    #                         gebv_blup_BGLR = bglr_gebv, 
    #                         #  gev_blup_GAPIT = gapit_out$BLUP,
    #                         gebv_gam = gen_dat_moms$prs_best_scaled),
    #              ellipses = F)
    # dev.off()    
    # 
    # 
    
    
    
    
    
    
    
    
    
    
#  Simulate selection based on PRS score like in Hayden & SOD -------------

    
  ### Run full model and save out residuals for cluster analysis
    fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\")  + s(accession, bs = \"re\") + s(locality, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC1_clim, bs=\"cr\") + s(PC2_clim, bs=\"cr\")"))
    
    # Gam with interaction effect
    gam_snp_all2 = bam(formula = formula(fixed_effects),
                      data = dat_snp_all,
                      discrete = FALSE, 
                      nthreads = 8,
                      method = "fREML",
                      family = "tw",
                      control = list(trace = FALSE))
    
    summary(gam_snp_all2)
    
    dat_snp_all2 <- dat_snp_all
    
    dat_snp_all2$rgr_resids <- resid(gam_snp_all2)
    hist(dat_snp_all2$rgr_resids, breaks = 50)
    
    ## Mean center the residuals
    dat_snp_all2$rgr_resids <-  dat_snp_all2$rgr_resids - mean(dat_snp_all2$rgr_resids)
    summary(dat_snp_all2$rgr_resids)
    
    
    dat_snp_gebv2 <- dat_snp_all2 %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
      dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
                    rgr, rgr_resids, 
                    section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
      dplyr::left_join(., dplyr::select(gen_dat_moms, accession, gebv_best, gebv_best_scaled)) %>%
      dplyr::mutate(accession = factor(accession)) %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all))

   # Set parameters
   set.seed(1)
   sample_size = 50
   n_reps <- 10000
   
   # Indices for top and bottom quartiles
   upper_ind <- which(dat_snp_gebv2$gebv_best_scaled >= 1 & 
                        dat_snp_gebv2$tmax_sum_dif_unscaled > 0)
   length(upper_ind)
   
   lower_ind <- which(dat_snp_gebv2$gebv_best_scaled <= 
                        quantile(dat_snp_gebv2$gebv_best_scaled, .25) 
                       & dat_snp_gebv2$tmax_sum_dif_unscaled > 0)
   length(lower_ind)
   
   random_ind <- which(dat_snp_gebv2$tmax_sum_dif_unscaled > 0)
   length(random_ind)
   
   climate_match <- dat_snp_gebv2  %>%
     dplyr::filter(tmax_sum_dif_unscaled <= .5 & tmax_sum_dif_unscaled >= -.5)  
   dim(climate_match)
   
   climate_colder <- dat_snp_gebv2  %>%
     dplyr::filter(tmax_sum_dif_unscaled <= -2)  
   dim(climate_colder)
   
   # Initialize output
   output_df <- tibble(rep = 1:n_reps,
                           `Upper quartile GEBV` = NA,
                           `Lower quartile GEBV` = NA,
                           `None` = NA,
                           `Climate Match` = NA,
                           `Climate Colder` = NA)
   
   
   for(x in 1:n_reps){
     
     if(x %% 500 == 0){
      cat("Working on x: ", x, " .... \n")
     }
     
     upper_sample <- sample(upper_ind, size = sample_size)
     lower_sample <- sample(lower_ind, size = sample_size)
     avg_sample <- sample(random_ind, size = sample_size)
     climate_match_sample <- sample(1:nrow(climate_match), size = sample_size)
     climate_colder_sample <- sample(1:nrow(climate_colder), size = sample_size)
     
     
     output_df$`Upper quartile GEBV`[x] <-     mean(dat_snp_gebv2$rgr_resids[upper_sample])
     output_df$`Lower quartile GEBV`[x] <-     mean(dat_snp_gebv2$rgr_resids[lower_sample])
     output_df$None[x] <-                      mean(dat_snp_gebv2$rgr_resids[avg_sample])
     output_df$`Climate Match`[x] <-            mean(climate_match$rgr_resids[climate_match_sample])
     output_df$`Climate Colder`[x] <-          mean(climate_colder$rgr_resids[climate_colder_sample])

   }
   
   
  # Plot results 
   output_df %>%
     gather(key = type, value = rgr_resid, -rep) %>%
     mutate(type = factor(type, levels = c("None", "Climate Match", "Climate Colder",
                                           "Upper quartile GEBV", "Lower quartile GEBV"))) %>%
     dplyr::filter(type != "Lower quartile GEBV") %>%
     ggplot(., aes(x = type, y = rgr_resid, fill = type)) + 
     geom_hline(yintercept = 0, lty = 2)+
    # geom_violin(draw_quantiles = c(.5)) +
     geom_boxplot(outlier.shape = NA) + 
     xlab("Choice criteria") +
     ylab("Adjusted Relative Growth Rate") + 
     scale_fill_manual(values =  c("#227eb8", "#fdc086", "#ffff99", "#7fbf7b")) + # "#af8dc3")) +
     scale_x_discrete(labels = c("Random", "Climate match \n -0.5° to 0.5°\nTmax transfer", 
                                 "Climate match \n <= -2°\nTmax transfer",
                                 "GEBV \n >= 1 SD")) +
                                 #"GEBV \n Lower quartile")) + 
     scale_y_continuous(breaks = seq(-.25, .25, by = 0.05),
                        limits = c(-.1, .15)) +
     theme_bw(8) + 
     theme(panel.border = element_blank(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           axis.line = element_line(colour = "black"),
           axis.text.x = element_text(size = 8),
           axis.title.x = element_text(size = 10),
           axis.title.y = element_text(size = 10),
           legend.position = "none") +
     NULL
   
   # Save output
     # ggsave(filename = paste0("./figs_tables/fig2/Figure 2 - rgr predictions by choice ",
     #                          Sys.Date(), ".pdf"),
     #        units = "cm",
     #        width = 9, height = 8)

   
   
  ## Comparing different scenarios 
  growth_by_scenario <-  output_df %>%
     gather(key = type, value = rgr_resid, -rep) %>%
     mutate(type = factor(type, levels = c("None", "Climate Match", "Climate Colder",
                                           "Upper quartile GEBV", "Lower quartile GEBV"))) %>%
     dplyr::filter(type != "Lower quartile GEBV") %>%
     group_by(type) %>%
     summarise_all(mean)
  
  growth_by_scenario
  
  growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "Climate Match"]
  
  growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "Climate Colder"]
  
  growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "None"]
  
  
  
  
  
  
  
  
  
   
# Plotting current distribution of GEBV  --------------------------

  library(raster)
  library(factoextra)
  
  # Shape outlines
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range_extended <- readShapePoly("./data/gis/valley_oak_range/qlobata_extended.shp", delete_null_obj = TRUE,
                                         proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Elevation map
  dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
  
  slope = raster::terrain(dem, opt='slope')
  aspect = raster::terrain(dem, opt='aspect')
  hill = hillShade(slope, aspect, 40, 270)
  hill_cropped <- crop(hill, extent(-125, -113.5, 32, 42.5)) # crop based on california outline
  
  # Climate rasters
  tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
  
  tmin_winter <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")
  
  aet <- raster("./data/gis/climate_data/BCM/historical/1951-1980/aet1951_1980_ave_HST_1513103180/aet1951_1980_ave_HST_1513103180.tif")
  cwd <- raster("./data/gis/climate_data/BCM/historical/1951-1980/cwd1951_1980_ave_HST_1513035274/cwd1951_1980_ave_HST_1513035274.tif")
  bioclim_04 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_04.tif")
  bioclim_15 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_15.tif")
  bioclim_18 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_18.tif")
  bioclim_19 <- raster("./data/gis/climate_data/BCM/historical/1951-1980/bioclim_19.tif")
  
  
  # Function to downscale and project raster 
  process_raster <- function(rast){
    rast <- aggregate(rast, fact = 5) # Make smaller by averaging across cells
    rast_latlon <- projectRaster(rast, crs = CRS("+proj=longlat +datum=WGS84"))
    plot(rast_latlon)
    return(rast_latlon)
  }
  
  # Aggregate and project tmax_raster
  tmax_rast   <- process_raster(tmax_rast)
  tmin_winter <- process_raster(tmin_winter)
  aet         <- process_raster(aet)
  cwd         <- process_raster(cwd)
  bioclim_04  <- process_raster(bioclim_04)
  bioclim_15  <- process_raster(bioclim_15)
  bioclim_18  <- process_raster(bioclim_18)
  bioclim_19  <- process_raster(bioclim_19)
  
  # Match up resolution and extent of DEM to climate rasters
  dem <- resample(dem, tmax_rast) 
  bioclim_04 <- resample(bioclim_04, tmax_rast) 
  bioclim_15 <- resample(bioclim_15, tmax_rast) 
  bioclim_18 <- resample(bioclim_18, tmax_rast) 
  bioclim_19 <- resample(bioclim_19, tmax_rast) 
  
  compareRaster(tmax_rast, dem) # Make sure rasters are equal
  
  # Find lat long points
  latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))
  
  gam_dist_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                                 longitude = latlon[, 1],
                                 latitude = latlon[, 2],
                                 tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                                 tmin_winter = raster::extract(tmin_winter, 1:ncell(tmin_winter)),
                                 aet = raster::extract(aet, 1:ncell(aet)),
                                 cwd = raster::extract(cwd, 1:ncell(cwd)),
                                 elevation = raster::extract(dem, 1:ncell(dem)),
                                 bioclim_04 = raster::extract(bioclim_04, 1:ncell(bioclim_04)),
                                 bioclim_15 = raster::extract(bioclim_15, 1:ncell(bioclim_15)),
                                 bioclim_18 = raster::extract(bioclim_18, 1:ncell(bioclim_18)),
                                 bioclim_19 = raster::extract(bioclim_19, 1:ncell(bioclim_19)),
                                 random = rnorm(1:ncell(dem)))
  # gam_dist_rast_df
  
  dim(gam_dist_rast_df)
  
 
  
  
  # Remove NAs from dataframe or random forest will complain
  gam_dist_rast_df_nona <-  gam_dist_rast_df %>%
    dplyr::filter(!is.na(tmax_sum)) %>%
    dplyr::filter(!is.na(elevation)) %>%
    dplyr::filter(!is.na(aet)) %>%
    dplyr::filter(!is.na(bioclim_04))
  
  summary(gam_dist_rast_df_nona)
  dim(gam_dist_rast_df_nona)
  
  
  # Choose climate variables   
  climate_vars_gam_dist <-  c("tmax_sum", 
                              "tmin_winter", 
                              "cwd",
                            #     "aet",
                            #     "bioclim_04",
                              "bioclim_15",
                              "bioclim_18",
                              # "bioclim_19", # R = 0.79 with aet; R = 0.84 with bioclim_18
                              "elevation",
                              "latitude",
                              "longitude") # R = -0.84 with latitude
  #   "random")
  
  # Check variance inflation factor
  car::vif(lm(snp_chr8_30946730 ~ . , 
              data = gen_dat_clim[, c("snp_chr8_30946730", climate_vars_gam_dist)]))
  
  pairs.panels(gen_dat_clim[, climate_vars_gam_dist])
  
  
  # Remove invididuals based on missing data
  # Count % of missing data for full gbs dataset
  missing_all <- data.frame(gen_dat_clim_all[ "id"],
                            missing = apply(dplyr::select(gen_dat_clim_all, snp_col_names),
                                            MARGIN = 1, function(x) sum(is.na(x)))/length(snp_col_names))
  
  # Get list of samples to exclude because of missing data
  exclude_based_on_missing <- missing_all$id[which(missing_all$missing > .20)]
  
  
  
  
  rda_df <-  gen_dat_clim %>%
    dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality),
                     by = c("id" = "gbs_name")) %>%
    dplyr::filter(!(id %in% exclude_based_on_missing))
  
  
  rda_df2 <- left_join(rda_df, dplyr::select(gen_dat_moms, accession, gebv_best_scaled) %>% 
                         distinct(accession, .keep_all = T) %>% 
                         mutate(accession = as.numeric(as.character(accession)))) %>%
     filter(!is.na(accession)) %>%
    filter(!is.na(gebv_best_scaled))
  
  pairs.panels(rda_df2[, climate_vars_gam_dist])
  
  car::vif(lm(snp_chr8_30946730 ~ . , 
              data = rda_df2[, c("snp_chr8_30946730", climate_vars_gam_dist)]))
  
  
  rda_df2$locality <- factor(rda_df2$locality)

  gam_dist <- bam(gebv_best_scaled ~ s(tmax_sum, bs = "cr") +
                    s(tmin_winter, bs = "cr") +
                    s(cwd, bs = "cr") + 
                    s(bioclim_15, bs = "cr") + 
                    s(bioclim_18, bs = "cr") +  
                    s(elevation, bs = "cr") + 
                    s(latitude, bs = "cr") + 
                    s(longitude, bs = "cr") + 
                    ti(latitude, longitude, bs = "cr"),
                  data = rda_df2,
                method = "fREML",
                discrete = FALSE)
  
  summary(gam_dist)
  
  
  # Save gam summary to file
    #  sink(file = paste0("./figs_tables/Table S3 - gam_dist_summary_",
    #                     Sys.Date(), ".txt"))
    #  summary(gam_dist)
    # anova(gam_dist)
    #  sink()

  #  visreg(gam_dist)
  
  gam_dist_rast_df_nona$locality <- "FHL"
  
  ## Predict across a raster and add to stack
  gam_dist_rast_df_nona$pred = predict(gam_dist,
                                       gam_dist_rast_df_nona, 
                                       type = "response")
  
  gam_dist_rast_df_temp <- left_join(gam_dist_rast_df,
                                     gam_dist_rast_df_nona[, c("cell_id", "pred")])
  
  rast <- tmax_rast # Initialize
  values(rast) <- gam_dist_rast_df_temp$pred
  names(rast) <- "test" # Name the raster
  
  
  gebv_rast_gam <- mask(rast, lobata_range_extended)
  plot(gebv_rast_gam)
  
  hist(gebv_rast_gam, breaks = 50)
  

  
  # Choose climate variables
  climate_vars_core <- climate_vars[!grepl(pattern = "random|PC", x = climate_vars)]
  climate_vars_core
  length(climate_vars_core) # Should be 10


## Calculate number of top and bottom SNPs for moms
    dat_snp_gebv_mom <- gen_dat_clim_all %>%
      dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality), 
                       by = c("id" = "gbs_name")) %>%
      dplyr::select(id, accession, locality, longitude, latitude, climate_vars_core) %>%
      dplyr::filter(!id %in% exclude_based_on_missing) %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>% 
      dplyr::left_join(., dplyr::select(gen_dat_moms, id, gebv_best_scaled))
    
    head(dat_snp_gebv_mom)
    
    summary(dat_snp_gebv_mom$gebv_best_scaled)
    hist(dat_snp_gebv_mom$gebv_best_scaled, breaks = 50)
    
    # Average by locality
    dat_snp_gebv_locality_mean <- dat_snp_gebv_mom %>%
      dplyr::group_by(locality) %>% 
      summarise_all(mean)
 
    # Set up dataframe   
    dat_snp_gebv_locality_sp <- SpatialPointsDataFrame(dat_snp_gebv_locality_mean[, c("longitude", "latitude")],
                                                          data = dat_snp_gebv_locality_mean[, c("locality", "gebv_best_scaled", "latitude", "longitude")],
                                                          proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    

    library(raster)
    library(rasterVis)

    ## Read in elevation dem and create a hillshade for mapping
    dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")

    slope = raster::terrain(dem, opt='slope')
    aspect = raster::terrain(dem, opt='aspect')
    hill = hillShade(slope, aspect, 40, 270)

    hill_cropped <- crop(hill, extent(-125, -113.5, 32, 42.5)) # crop based on california outline

    ## Load in california outline
    cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))

    lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))

    lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                        proj4string = CRS("+proj=longlat +datum=WGS84"))

    lobata_range_extended <- readShapePoly("./data/gis/valley_oak_range/qlobata_extended.shp", delete_null_obj = TRUE,
                                           proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    
    
    # Function to downscale and project raster 
    process_raster <- function(rast){
      rast <- aggregate(rast, fact = 5) # Make smaller by averaging across cells
      rast_latlon <- projectRaster(rast, crs = CRS("+proj=longlat +datum=WGS84"))
      rast_latlon <- mask(rast_latlon, lobata_range_extended)
      plot(rast_latlon)
      return(rast_latlon)
    }
    
   
  ## Creates vector with list of file paths to all .tif raster files
    ## Directory path where future climate scenarios are located
    dir_name <- "./data/gis/climate_data/BCM/future/"
    raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    raster_files <- raster_files[grep("*[[:digit:]].tif$", raster_files)] # Only those with .tif extension
    raster_files <- raster_files[grep("jja_ave", raster_files)] # Only tmax sum files

    ## Select only 85 scenario
    raster_files <- raster_files[grep("rcp85", raster_files)]
    
    raster_files
    
    future_stack <- stack()
    x = 1
    
    ## Loop through raster files and add to stack
    for(file in raster_files){
      
      # Load in raster
      raster_temp <- raster(file)
      
      cat("Working on file:", file, "...\n")
      
      ## Aggregate, project to lat long, and crop to area of interest
      raster_temp <- process_raster(raster_temp)
      
      # Check to see if extent is different, if it is - change it
      if(x > 1 & !compareRaster(future_stack, raster_temp, stopiffalse = FALSE)){
        extent(raster_temp) <- extent(future_stack)
      }
      
      future_stack <- raster::stack(future_stack, raster_temp)
      x = x + 1
    }
    
    future_stack
    
    ## Clean up names
    names(future_stack) <-  unlist(lapply(strsplit(x = gsub(pattern = "tmx2070_2099jja_ave_", 
                                                            replacement = "", names(future_stack)),
                                                   split = "_1"), "[", 1))
    
    # Load in contemporary Tmax raster
    tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
    
    tmax_rast <- process_raster(tmax_rast)
    
    
  # Calculate difference between future and current
    future_stack_tmax_dif <- future_stack - tmax_rast
    names(future_stack_tmax_dif) <- names(future_stack) # Reassign name
    
    # Scale values
    future_stack_tmax_dif_scaled <- (future_stack_tmax_dif - 
                                       scaled_var_means_gbs_all["tmax_sum_dif"]) /
                                      scaled_var_sds_gbs_all["tmax_sum_dif"]
   
    
  ## Resample gebv Raster to match with tmax rasters

      # From gam predictions
      gebv_rast_mean <- resample(gebv_rast_gam, tmax_rast) # From gam predictions
      
      compareRaster(tmax_rast, gebv_rast_mean) # Make sure rasters are equal
      
      
      # Set up scenario where we pick the highest gebv score within a XXXX radius of the point

      ## Change to equal area projection and resolution to 1km
        gebv_rast_mean2 <- projectRaster(from = gebv_rast_mean,
                                        crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m", 
                                        res = 1000) # Resolution of 1km
        gebv_rast_mean2
        
        
        buffer = focalWeight(gebv_rast_mean2, d =  50000, # 25 km radius
                             type = "circle")
      
        buffer[buffer > 0] <- 1 # Set all non 0 buffer weights to 1
        buffer
        dim(buffer)
        
        # Run focal statistics to find max values within radius
          gebv_rast_mean_focal = focal(gebv_rast_mean2, w = buffer, fun = max, na.rm = T)
          
          plot( gebv_rast_mean_focal)
          
        # Project back to lat long
          gebv_rast_mean_focal <- projectRaster(from =  gebv_rast_mean_focal, 
                                to = gebv_rast_mean)
    
          gebv_rast_mean_focal <- mask( gebv_rast_mean_focal, lobata_range_extended)
          plot( gebv_rast_mean_focal)
        
        # Check to make sure all positive values
          hist( gebv_rast_mean_focal -  gebv_rast_mean)
          summary(values( gebv_rast_mean_focal -  gebv_rast_mean))
        
        # Check to make sure max values aren't above
          summary(values( gebv_rast_mean))
          summary(values( gebv_rast_mean_focal))
          
      # Save raster    
        gebv_rast_max <-  gebv_rast_mean_focal
  
        
   ## Get null prediction of height with no change in tmax

     future_null <- tmax_rast
     names(future_null) <- "tmax_sum_dif"
     
     # Need to scale tmax when actual transfer distance is 0
     future_null[!is.na(future_null)] <- forward_transform(0, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all)
     
     future_null_stack <- stack(future_null, gebv_rast_mean) #
     
     
     names(future_null_stack) <- c("tmax_sum_dif", "gebv_best_scaled") # Rename so predict function works
     future_null_stack
     
     future_null_height_gebv <- predict(future_null_stack,
                                       model = gam_gebv_all,
                                       progress = "text",
                                       type = "response")
     
     future_null_height_gebv
     plot(future_null_height_gebv)
     summary(future_null_height_gebv)
     
     ## Add baseline data
     
     # Make data frame for prediction - one for top and one for bottom
     newdata_rast_const <- data.frame(section_block = "Block1_1",
                                      height_2014 = 0,
                                      accession = "6",
                                      PC1_gen = 0, PC2_gen = 0, 
                                      PC1_clim = 0, PC2_clim = 0,
                                      locality = "FHL")
     
     
     future_null_height_gebv <- future_null_height_gebv +  predict(future_null_stack,
             model = gam_snp_all,
             const = newdata_rast_const,
             progress = "text",
             type = "response")
     

   # Initialize stack for height changes
       future_stack_height_change_gebv <- future_stack_tmax_dif
       future_stack_height_change_gebv_selection <- future_stack_tmax_dif
   
  
    # Loop over scenarios
    for(scenario in 1:nlayers(future_stack_tmax_dif)){
      
      cat("Working on scenario:", names(future_stack_tmax_dif)[scenario], "... \n")
      
      ## Use scaled raster
        rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
        rast_temp <- stack(rast_temp, gebv_rast_mean)
        
        names(rast_temp) <- c("tmax_sum_dif", "gebv_best_scaled") # Rename so predict function works
        rast_temp
        
        # Calculate baseline Rgr
        baseline_temp <- predict(rast_temp,
                                 model = gam_snp_all,
                                 const = newdata_rast_const,
                                 progress = "text",
                                 type = "response")
        
        
        
        rast_temp_height_gebv <- predict(rast_temp,
                                        model = gam_gebv_all,
                                        progress = "text",
                                        type = "response")
        
        ## Add baseline RGR
        rast_temp_height_gebv <- rast_temp_height_gebv + baseline_temp
        
      
      ## Use selection raster
        rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
        rast_temp <- stack(rast_temp, gebv_rast_max)
        
        names(rast_temp) <- c("tmax_sum_dif", "gebv_best_scaled") # Rename so predict function works
        rast_temp
        
        rast_temp_height_gebv_selection <- predict(rast_temp,
                                        model = gam_gebv_all,
                                        progress = "text",
                                        type = "response")
        # Add baseline RGR
        rast_temp_height_gebv_selection <- rast_temp_height_gebv_selection + baseline_temp

      # Dif bw top and null 
     rast_temp_height_change_gebv <- (rast_temp_height_gebv - future_null_height_gebv) / abs(future_null_height_gebv) * 100
     # rast_temp_height_change_gebv <- rast_temp_height_gebv - future_null_height_gebv 
     #   rast_temp_height_change_gebv <- rast_temp_height_gebv
      plot(rast_temp_height_change_gebv, main = scenario)
      future_stack_height_change_gebv[[scenario]] <- rast_temp_height_change_gebv
      names(future_stack_height_change_gebv)[scenario] <- names(future_stack_tmax_dif)[scenario] 
      
      
      # Selection
   #  rast_temp_height_change_gebv_selection <- rast_temp_height_gebv_selection - future_null_height_gebv
   #   rast_temp_height_change_gebv_selection <- rast_temp_height_gebv_selection
      rast_temp_height_change_gebv_selection <- (rast_temp_height_gebv_selection - future_null_height_gebv) / abs(future_null_height_gebv) * 100
      plot(rast_temp_height_change_gebv_selection,
           main = paste0(scenario, " | selection"))
      future_stack_height_change_gebv_selection[[scenario]] <- rast_temp_height_change_gebv_selection
      names(future_stack_height_change_gebv_selection)[scenario] <- names(future_stack_tmax_dif)[scenario] 
      

    } # End scenario loop
     
    
    # Average across scenarios
     future_stack_height_change_gebv_avg <- mean(future_stack_height_change_gebv, na.rm = TRUE)
     
     future_stack_height_change_gebv_selection_avg <- mean(future_stack_height_change_gebv_selection,
                                                          na.rm = TRUE)
     
     hist(future_stack_height_change_gebv_avg)
     hist(future_stack_height_change_gebv_selection_avg)
     
     
     height_change_dif = future_stack_height_change_gebv_selection_avg - future_stack_height_change_gebv_avg
     plot(height_change_dif)
     hist(height_change_dif)
     
     # Find min across both rasters
     rast_min <- min(c(min(values(future_stack_height_change_gebv_avg), na.rm = T), 
                    min(values(future_stack_height_change_gebv_selection_avg), na.rm = T)))
     rast_max <- max(c(max(values(future_stack_height_change_gebv_avg), na.rm = T), 
                       max(values(future_stack_height_change_gebv_selection_avg), na.rm = T)))
     
     
     hill_cropped <- crop(hill, extent(-125, -118, 33, 42.1)) # crop based on california outline
   
     
     
     library(tmap)
     
     
  # Plot of baseline gebv scores  
    # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +

       # Add kriged raster
     #  tm_shape(future_stack_height_change_gebv_selection_avg) +
       tm_shape(gebv_rast_mean) +
       tm_raster(palette = "PRGn",
             #    breaks = seq(rast_min, rast_max, length.out = 10),
                breaks = seq(-2, 2, by = 0.5),
                 alpha = .95,
              midpoint = NA) +

       tm_shape(lobata_range_extended) +
       tm_borders("black") +
       tm_shape(lobata_range) + # Fine scale range
       tm_borders("black", lwd = .15) +

       # # Bubble plots
       tm_shape(dat_snp_gebv_locality_sp) +
       tm_bubbles(size = .025,
                  shape = 21,
                  col = "black",
                  border.alpha = 0,
                  midpoint = NA)
     
     tm

     tmap_save(tm, filename = paste0("./figs_tables/fig3/Figure 3 - gebv map ", Sys.Date(), ".pdf"),
               width = 12, height =12, units = "cm", useDingbats=FALSE)
     
    
     
     
      
   ## Plot of predicted growth rates in 2080 based on current gebv distribution
     
   # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +
       
       # Add kriged raster
       tm_shape(future_stack_height_change_gebv_avg) +
       tm_raster(palette = "YlGn",
               #      breaks = seq(rast_min, rast_max, length.out = 20),
                 breaks = seq(-10, 10, by = 2),
                 alpha = .95,
                 midpoint = NA) +
       
       tm_shape(lobata_range_extended) +
       tm_borders("black") +
       tm_shape(lobata_range) + # Fine scale range
       tm_borders("black", lwd = .15)  +
       
       # # # Bubble plots
       tm_shape(dat_snp_gebv_locality_sp) +
       tm_bubbles(size = .025,
                  shape = 21,
                  col = "black",
                  border.alpha = 0,
                  midpoint = NA)
     tm
     
     
     # tmaptools::palette_explorer()
     
     # Save plot
     tmap_save(tm, filename = paste0("./figs_tables/fig3/Figure 3 - growth 2080 map ", Sys.Date(), ".pdf"),
               width = 12, height =12, units = "cm", useDingbats=FALSE)
     
     
     
     
     ## Plot of predicted growth rates in 2080 based on MAX gebv
     
     # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +
       
       # Add kriged raster
       tm_shape(future_stack_height_change_gebv_selection_avg) +
       tm_raster(palette = "YlGn",
                 #      breaks = seq(rast_min, rast_max, length.out = 20),
                 breaks = seq(-10, 10, by = 2),
                 alpha = .95,
                 midpoint = NA) +
       
       tm_shape(lobata_range_extended) +
       tm_borders("black") +
       tm_shape(lobata_range) + # Fine scale range
       tm_borders("black", lwd = .15)  +
       
       # # # Bubble plots
       tm_shape(dat_snp_gebv_locality_sp) +
       tm_bubbles(size = .025,
                  shape = 21,
                  col = "black",
                  border.alpha = 0,
                  midpoint = NA)
     tm
     
     
     # tmaptools::palette_explorer()
     
     # Save plot
     tmap_save(tm, filename = paste0("./figs_tables/fig3/Figure 3 - growth gebv selection 2080 map ", Sys.Date(), ".pdf"),
               width = 12, height =12, units = "cm", useDingbats=FALSE)
     

  

