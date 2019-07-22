
# Run model setup ---------------------------------------------------------

  source("./04b-setup-models-cluster.R")

# Read in cluster run output ----------------------------------------------

  # Path to results
  # Used in first submission to PNAS - "run_475547_tmax_sum_dif_training_set_resids_7030_2019_02_26"
#  path <- "run_475547_tmax_sum_dif_training_set_resids_7030_2019_02_26"

  # With siblings in both training and testing sets July 16 2019
   path <- "run_389456_tmax_sum_dif_rgr_family_acrossCV_reduced_model"
  
  path_to_summaries <- paste0("./output/", path, "/model_summaries/")
  path_to_predictions <- paste0("./output/", path, "/model_predictions/")

  
 # Read in files
    sum_df_raw <- plyr::ldply(list.files(path_to_summaries, full = TRUE), read_csv)
    
    dim(sum_df_raw)
    sum_df_raw
    
    summary(sum_df_raw)
  
  
  ## Missing snps
  # Find if there indexes if needed to re-run manually or on the cluster
    which(!snp_col_names %in% sum_df_raw$snp)
    snp_col_names[which(!snp_col_names %in% sum_df_raw$snp)]
    length( which(!snp_col_names %in% sum_df_raw$snp))
    
    ## snp_chr5_22319422 doesn't have 2 levels
    
    ## Write snps that were run to file to rerun missing snps on cluster
      # write_csv(data.frame(snp = snp_col_names[which(snp_col_names %in% sum_df_raw$snp)]),
      #           path = "./output/skip_these_SNPs_2018_12_28.csv")
    
  # Make sure columns are in right format
    sum_df_raw <- sum_df_raw %>%
      mutate(n_gen0 = as.numeric(n_gen0),
             n_gen1 = as.numeric(n_gen1),
             n_gen2 = as.numeric(n_gen2),
             # height_change_gen_0 = as.numeric(height_change_gen_0),
             # height_change_gen_1 = as.numeric(height_change_gen_1),
             # height_change_gen_2 = as.numeric(height_change_gen_2),
             # p_val_gen_0 = as.numeric(p_val_gen_0),
             # p_val_gen_1 = as.numeric(p_val_gen_1),
             # p_val_gen_2 = as.numeric(p_val_gen_2),
             p_val_gen_int = as.numeric(p_val_gen_int),
             cv_cor_gen_0 = as.numeric(cv_cor_gen_0),
             cv_cor_gen_1 = as.numeric(cv_cor_gen_1),
             cv_cor_gen_2 = as.numeric(cv_cor_gen_2))
    
    str(sum_df_raw)
    
    
  # For backup
  sum_df <- sum_df_raw
  
  # Filter out models that did not converge
  table(sum_df_raw$converged)
  
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
  
  
# Convert to long format and set filters ---------------------------------------------------

  ## Convert to from wide to long format
  n_gen_long <- sum_df %>%
    dplyr::select(snp, n_gen0, n_gen1, n_gen2, acc_pred, n) %>%
    gather(key = "genotype", value = "n_gen", -snp, -acc_pred, -n) %>%
    mutate(genotype = gsub("n_gen", "", genotype))
  
  head(n_gen_long)
  
  cv_cor_long <- sum_df %>%
    dplyr::select(snp, cv_cor_gen_0, cv_cor_gen_1, cv_cor_gen_2) %>%
    gather(key = "genotype", value = "cv_cor", -snp) %>%
    mutate(genotype = gsub("cv_cor_gen_", "", genotype))
  
  head(cv_cor_long)
  

## Join all long dataframes together
  sum_df_long <- left_join(n_gen_long, cv_cor_long)
  
  head(sum_df_long)
  dim(sum_df_long)
  str(sum_df_long)
  
 
 # Filter out based on number of gene copies
    hist(sum_df_long$n_gen, breaks = 50)

    gen_copies_thresh <- 1 # Used to be 100

   sum_df_long <- sum_df_long %>%
      dplyr::filter(n_gen >= gen_copies_thresh)
    
  ## Filter out heterozygotes?
   table(sum_df_long$genotype)
   # sum_df_long <- sum_df_long %>%
  #    dplyr::filter(genotype != "1")
  
  ## Remove those that predictions weren't made on accession 1
   # table(sum_df_long$acc_pred)
   #  sum_df_long <- sum_df_long %>%
   #    dplyr::filter(acc_pred == 1)
    
  # Subset prediction dataframe to just predictions of SNPs that pass filters
    pred_df_sub <- pred_df_raw %>%
      mutate(snp_gen = paste0(snp, genotype)) %>%
      # Snp-genotype must be in sum_df_long dataframe
      filter(snp_gen %in% paste0(sum_df_long$snp, sum_df_long$genotype)) 
    
    dim(pred_df_sub)
    
    
    ## Calculate height change in warmer temperatures for each snp and genotype
    base0 <- mean(pred_df_sub$pred[pred_df_sub$tmax_sum_dif_unscaled == 0])

    pred_df_long =  pred_df_sub %>%
      dplyr::group_by(snp, genotype) %>%
      do(data.frame(height_change_warmer = mean((.$pred[.$tmax_sum_dif_unscaled > 0]  -
                                               .$pred[.$tmax_sum_dif_unscaled == 0] ) /
                                               abs(.$pred[.$tmax_sum_dif_unscaled == 0]) ) * 100,
                    height_change_warmer_base0 = mean((.$pred[.$tmax_sum_dif_unscaled > 0]
                                                       - base0 ) /  abs(base0)) * 100,
                    height_change_warmer_absolute = mean(.$pred[.$tmax_sum_dif_unscaled > 0])))
    
    head(pred_df_long)
    dim(pred_df_long)
    
    summary(pred_df_long$height_change_warmer_absolute)
    summary(pred_df_long$height_change_warmer)
    summary(pred_df_long$height_change_warmer_base0)
    
    
    # plot(pred_df_long$height_change_warmer,
    #      pred_df_long$height_change_warmer_base0, pch = "." )
    # abline(a = 0, b = 1)
    # cor.test(pred_df_long$height_change_warmer,pred_df_long$height_change_warmer_base0)
    # 
    # plot(pred_df_long$height_change_warmer_absolute,
    #      pred_df_long$height_change_warmer_base0, pch = "." )
    # abline(a = 0, b = 1)
    # cor.test(pred_df_long$height_change_warmer_absolute,
    #          pred_df_long$height_change_warmer_base0 )
    # 
    # head(pred_df_long)
    # dim(pred_df_long)
    
  
 # Join prediction dataframe to summary dataframe
    sum_df_long_w_preds <- left_join(sum_df_long,
                                     pred_df_long)
    
  
# Correct for multiple testing

  library(fdrtool)

  # F values of interaction term
  hist(sum_df$f_val_gen_int, breaks = 50)

  # fdr_fvals = fdrtool(c(sum_df$f_val_gen_int),
  #                statistic = "normal", plot = FALSE)
  
  fdr_fvals = fdrtool(c(sum_df$p_val_gen_int),
                      statistic = "pvalue", plot = TRUE)
  # 
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
  
  par(mfrow = c(2,2))
  
  # P values
  hist(fdr_fvals$pval, breaks = 40, col = "steelblue2",
       xlab = "p-value", main = "Uncorrected p-values", las = 1)
  qq(fdr_fvals$pval, main = "Uncorrected p-values", las = 1)
  
  # Q values
  hist(fdr_fvals$qval, breaks = 40, col = "steelblue2",
       xlab = "q-value", main = "q-values", las = 1)
  qq(fdr_fvals$qval, main = "q-values", las = 1)
  sort(fdr_fvals$qval, decreasing = T)
  
  # Plot figure
  ## Need to hand-edit the axis labesl for q values to say Q instead of P
      # dev.copy(pdf, paste0("./figs_tables/Figure S6 - histogram and qqplots of p values ", Sys.Date(), ".pdf"),
      #          width = 7, height = 7, useDingbats = FALSE)
      # dev.off()

  par(mfrow = c(1,1))


  # Join back to main dataframe
  sum_df$q_val <- fdr_fvals$qval
  sum_df_long_w_preds <- dplyr::left_join(sum_df_long_w_preds,
                                          dplyr::select(sum_df, snp, q_val),
                                          by = "snp")
  summary(sum_df_long_w_preds$q_val)
  
  
  
#    ## Plot interactions at particular genotype
#    # pred_snp = dplyr::left_join(dplyr::select(man_df, SNP), pred_df_raw, by = c("SNP" = "snp"))
#    # pred_snp <- pred_snp %>%
#    #   rename(snp = SNP)
#    # table(pred_snp$snp)
#    # 
#    # snp_filter = "snp_chr8_7837444" # For example of strong interaction
#    # snp_filter = "snp_chr2_10328592"
#    # pred_snp %>%
#    #   dplyr::filter(snp == snp_filter)  %>%
#    # ggplot(., aes(x = tmax_sum_dif_unscaled, y = pred, col = factor(genotype))) + 
#    #   geom_vline(aes(xintercept = 0), lty = 2, size = .25) +
#    #   geom_vline(xintercept = 4.8, size = .25) +
#    #   geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96 * se, fill = factor(genotype)), alpha = .15, color = NA,
#    #               show.legend = F) +
#    #   geom_line(alpha = 1,
#    #             show.legend = T, size = 1) + 
#    #   ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
#    #   xlab("Tmax transfer distance") +
#    #   theme_bw(15) + 
#    #   scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) + 
#    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#    #         panel.grid.minor = element_blank(), 
#    #         axis.line = element_line(colour = "black"))   + 
#    #   labs(color = "Genotype")
#    # 
#    # sum_df_long_w_preds[sum_df_long_w_preds$snp == snp_filter,]
#    # table(dat_snp[, snp_filter])
#    
   
  
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
                    y = -log10(q_val),
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
     scale_y_continuous(breaks = seq(0, 18, 2), limits = c(0, 17)) + 
     ylab("-log10(q-value)") + 
     xlab("Chromosome") + 
     theme_bw(8) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
          # legend.position = "none") +
             )+
     # axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
     NULL
   
   # # Save as file - PNG
     # ggsave(paste0("./figs_tables/fig2/Figure 2 - manhattan plot interaction term wo labels_",
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
                size = 1) +
     
     # Horizontal line
     geom_hline(yintercept = mean(man_df$height_change_warmer_absolute, na.rm = TRUE), 
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

   #   # Save as PDF for legend
   #   ggsave(paste0("./figs_tables/fig2/Figure 2 - manhattan plot growth predictions LEGEND",
   #                 Sys.Date(), ".pdf"),
   #          width = 10, height = 2)
   # 
   # 

  

# Calculate polygenic risk score ------------------------------------------

  # Calculate polygenic risk score for each maternal genotype
  gen_dat_moms <- dplyr::select(gen_dat_clim, id, accession, snp_col_names)
  dim(gen_dat_moms)
  
  # Choose SNPs
   n_snps_values <- c(100, 500, seq(1000, 12000, by = 1000))

  gene_copy_thresholds <- c(1, 25, 50, 100, 250)
  
  # Initialize dataframe for output
  prs_df_summary <- tibble(n_snps = n_snps_values,
                               r2 = NA)
  
 # Order snps based on lowest Q values
   snp_list_ordered <- sum_df$snp[order(sum_df$q_val)]

  # Initialize
    gen_dat_moms$prs <- NULL
    gen_dat_moms$prs_count <- NULL
  
  x = 1

  for(snp_name in snp_list_ordered){

    if(x %% 500 == 0){
      cat("Working on snp number: ", x, "... \n")
    }

    # Subset to genotypes
    sub <- gen_dat_moms[ , snp_name]

    preds_sub <- dplyr::filter(sum_df_long_w_preds, snp == snp_name) %>%
                dplyr::select(genotype, height_change_warmer_absolute, n_gen, q_val)
    
    ## Filter out genotypes below a count value
  #  preds_sub <- dplyr::filter(preds_sub, n_gen > gene_copy_thresh)
    
    counted <- preds_sub$genotype
    
    # Add in median value for missing genotypes
    if(!"0" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "0",
                                  height_change_warmer_absolute = 0,
                                  n_gen = NA,
                                  q_val = NA))
    }
    if(!"1" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "1",
                                               height_change_warmer_absolute = 0,
                                               n_gen = NA,
                                               q_val = NA))
    }
    if(!"2" %in% preds_sub$genotype){
      preds_sub <- rbind(preds_sub, data.frame(genotype = "2",
                                               height_change_warmer_absolute = 0,
                                               n_gen = NA,
                                               q_val = NA))
    }
    
     prs <-  case_when(
        sub == "0" ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "0"],  
        sub == "1"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "1"],
        sub == "2"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "2"],
        TRUE ~ 0
      )
   
   prs_count <- ifelse(as.data.frame(sub)[, 1] %in% counted, yes = 1, no = 0)
   
   if(x == 1){
     gen_dat_moms$prs <- prs
     gen_dat_moms$prs_count <- prs_count
   } else {
     gen_dat_moms$prs <- gen_dat_moms$prs + prs
     gen_dat_moms$prs_count <- gen_dat_moms$prs_count + prs_count

   }
   
   if(x %in% n_snps_values){
     cat("\n||--------------------------------||\n")
     cat("Running model for:", x, " SNPS... \n")
     
     # Adjust for missing data
     gen_dat_moms$prs_by_count <- gen_dat_moms$prs / gen_dat_moms$prs_count
     
     # Save prs score as its own variable for later use
     gen_dat_moms[, paste0("prs_", x)] <-  gen_dat_moms$prs_by_count
     
     ### Join to testing dataframe
     dat_snp_testing_unscaled$accession <- as.numeric(as.character(dat_snp_testing_unscaled$accession))
     test2 <- dplyr::left_join(dat_snp_testing_unscaled, 
                               dplyr::select(gen_dat_moms, accession, prs_by_count), by = "accession")
     
     # Scale prs
     test2$prs_scaled <- (test2$prs_by_count - mean(test2$prs_by_count)) / sd(test2$prs_by_count)
     
     #  test2$prs_scaled <- rnorm(n = length(test2$prs_scaled))
     
     # Set up model terms
     fixed_effects_resids <- paste0(paste0("rgr_resids ~ s(prs_scaled, bs=\"cr\") + 
                                           s(tmax_sum_dif, by = prs_scaled, bs=\"cr\")"))
 
     # Gam with interaction effect
     gam_test_resids = bam(formula = formula(fixed_effects_resids),
                           data = test2,
                           discrete = TRUE, 
                           nthreads = 8,
                           method = "fREML",
                           control = list(trace = FALSE))
     
     print(summary(gam_test_resids))
     
  # Plot visualization
     title <-  paste(" # SNPS: ", x)
     visreg(gam_test_resids, xvar = "tmax_sum_dif", by = "prs_scaled", partial = F, overlay = T, main = title)
     
     cat(paste(title, "\n"))
         
         
    ## Save results   
     prs_df_summary$r2[prs_df_summary$n_snps == x] <- summary(gam_test_resids)$r.sq

   } # End model loop
 
   x = x + 1 # Move onto next SNP
   

  } # End snp loop
  
  
  prs_df_summary
  
  prs_df_summary %>%
    arrange(desc(r2))
  
  ## Plot relationship between # of SNPS used and variance explained in testing set
    ggplot(prs_df_summary, aes(x = n_snps, y = r2)) + 
      geom_point(size = 2) + 
      geom_line() + 
      xlab("Number of SNPs in GEBV estimation") +
      #ylab("R2 adjusted in testing set") +
      ylab(bquote('R'^2[adj]~'training test')) +
      scale_x_continuous(breaks = prs_df_summary$n_snps) +
      scale_y_continuous(breaks = seq(0, 0.02, by = 0.0025), limits = c(-0.0015, 0.015)) + 
      theme_bw(7) + 
      geom_hline(yintercept = 0, lty = 2) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.margin = (margin(0, 0, 0, 0, "cm")),
            axis.text.x = element_text(angle = 60, hjust = 1)) +
      NULL

    
    # # Save as file
    # ggsave(paste0("./figs_tables/Figure S3 - testing variance explained by number of snps",
    #                      Sys.Date(), ".png"),
    #          width = 10, height = 4, units = "cm", dpi = 300)
           #useDingbats=FALSE)
    
    
    
  
   # Set top PRS
    gen_dat_moms$prs_best <- gen_dat_moms$prs_11000
    
    
    prs_best_mean <- mean(gen_dat_moms$prs_best)
    prs_best_sd <- sd(gen_dat_moms$prs_best)
    gen_dat_moms$prs_best_unscaled <- gen_dat_moms$prs_best
    gen_dat_moms$prs_best_scaled <- (gen_dat_moms$prs_best - prs_best_mean) / prs_best_sd
  
    
    summary(gen_dat_moms$prs_best_scaled); hist(gen_dat_moms$prs_best_scaled, breaks = 50)
    
    
    
  
# Plot predictions of PRS -------------------------------------------------

    
    ## Calculate number of top and bottom SNPs 
    dat_snp_prs <- dat_snp_all %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
      dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
                    rgr, rgr_resids, 
                    section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
      dplyr::left_join(., dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled)) %>%
      dplyr::mutate(accession = factor(accession))
    
    head(dat_snp_prs)
    dim(dat_snp_prs)
    
    summary(dat_snp_prs$prs_best)
    summary(dat_snp_prs$prs_best_scaled)
   
    dat_snp_prs_unscaled <- dat_snp_prs
    
    
    # Set formula
     fixed_effects <- paste0(paste0("rgr_resids ~ s(prs_best_scaled, bs=\"cr\") + 
                                           s(tmax_sum_dif, by = prs_best_scaled, bs=\"cr\")"))
 
    # Gam with interaction effect
    gam_prs_all = bam(formula = formula(fixed_effects),
                  data = dat_snp_prs,
                  discrete = TRUE, 
                  nthreads = 8,
                  method = "fREML",
               #   family = "tw",
                  control = list(trace = FALSE))
    
    summary(gam_prs_all)
    
    ## Save output
    
   # Save gam summary to file
   #  sink(file = paste0("./figs_tables/Table S2 - prs_all_gam_model_summary_",
   #                     Sys.Date(), ".txt"))
   #  summary(gam_prs_all)
   # # anova(gam_prs_all)
   #  sink()
    
    
    
    # Run models with testing and training set and output model summaries
        # 
        # # ## Training  
        #   dat_snp_training2 <- dat_snp_training %>%
        #     dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
        #     dplyr::left_join(.,  dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled))
        # 
        #   dim(dat_snp_training2)
        # 
        #   # Gam with interaction effect
        #   gam_prs_training = bam(formula = formula(fixed_effects),
        #                     data = dat_snp_training2,
        #                     discrete = TRUE,
        #                     nthreads = 8,
        #                     method = "fREML",
        #                     #   family = "tw",
        #                     control = list(trace = FALSE))
        # 
        #   summary(gam_prs_training)

        #   ## Save output
        #   
        #   # Save gam summary to file
        #    sink(file = paste0("./figs_tables/Table S2 - prs_training_gam_model_summary_",
        #                       Sys.Date(), ".txt"))
        #    summary(gam_prs_training)
        #   # anova(gam_prs_training)
        #    sink()
        #   
        # #### Testing   
        #    dat_snp_testing2 <- dat_snp_testing %>%
        #      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
        #      dplyr::left_join(.,  dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled))
        #    
        #    dim(dat_snp_testing2)
        #    
        #    # Gam with interaction effect
        #    gam_prs_testing = bam(formula = formula(fixed_effects),
        #                           data = dat_snp_testing2,
        #                           discrete = TRUE, 
        #                           nthreads = 8,
        #                           method = "fREML",
        #                           #   family = "tw",
        #                           control = list(trace = FALSE))
        #    
        #    summary(gam_prs_testing)
        #    
        #    ## Save output
        #    
        #    # Save gam summary to file
        #    sink(file = paste0("./figs_tables/Table S2 - prs_testing_gam_model_summary_",
        #                       Sys.Date(), ".txt"))
        #    summary(gam_prs_testing)
        #    # anova(gam_prs_testing)
        #    sink()   
       

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
                            prs_best_scaled = c(-1, 0, 1)) # SD away from average
    
    # Make predictions
    newdata$pred <- predict(gam_prs_all,
                            newdata = newdata,
                            se.fit = TRUE, type = "response")$fit
    
    newdata$se <- predict(gam_prs_all,
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
                      fill = factor(prs_best_scaled)), alpha = .35) +
      geom_line(aes(col = factor(prs_best_scaled))) + 
      scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
      scale_y_continuous(breaks = seq(0, 1.75, by = 0.25)) + 
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
    # ggsave(filename = paste0("./figs_tables/fig2/Figure 2 - prs growth ", Sys.Date(), ".pdf"),
    #        units = "cm", width = 8, height = 6)
    # 
    
    
    
    
    ## Calculating differences in growth rates
    growth_dif =  newdata %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                         means = scaled_var_means_gbs_all,
                                                         sds = scaled_var_sds_gbs_all))

    # Comparing high to mid
    (growth_dif$pred[growth_dif$prs_best_scaled == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
        growth_dif$pred[growth_dif$prs_best_scaled == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
      abs(growth_dif$pred[growth_dif$prs_best_scaled == 0 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100

    # Compare high to low
    (growth_dif$pred[growth_dif$prs_best_scaled == 1 & growth_dif$tmax_sum_dif_unscaled == 4.8] -
        growth_dif$pred[growth_dif$prs_best_scaled == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) /
     abs(growth_dif$pred[growth_dif$prs_best_scaled == -1 & growth_dif$tmax_sum_dif_unscaled == 4.8]) * 100

    
    
    
    
    
    
    
    
#  Simulate selection based on PRS score like in Hayden & SOD -------------

    
  ### Run full model and save out residuals for cluster analysis
    fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\")  + s(accession, bs = \"re\") + s(locality, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC1_clim, bs=\"cr\") + s(PC2_clim, bs=\"cr\")"))
    
    # Gam with interaction effect
    gam_snp_all2 = bam(formula = formula(fixed_effects),
                      data = dat_snp_all,
                      discrete = TRUE, 
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
    
    
    dat_snp_prs2 <- dat_snp_all2 %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
      dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
                    rgr, rgr_resids, 
                    section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
      dplyr::left_join(., dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled)) %>%
      dplyr::mutate(accession = factor(accession)) %>%
      dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all))

   # Set parameters
   set.seed(1)
   sample_size = 50
   n_reps <- 10000
   
   # Indices for top and bottom quartiles
   upper_ind <- which(dat_snp_prs2$prs_best_scaled >= 1 & 
                        dat_snp_prs2$tmax_sum_dif_unscaled > 0)
   length(upper_ind)
   
   lower_ind <- which(dat_snp_prs2$prs_best_scaled <= 
                        quantile(dat_snp_prs2$prs_best_scaled, .25) 
                       & dat_snp_prs2$tmax_sum_dif_unscaled > 0)
   length(lower_ind)
   
   random_ind <- which(dat_snp_prs2$tmax_sum_dif_unscaled > 0)
   length(random_ind)
   
   climate_match <- dat_snp_prs2  %>%
     dplyr::filter(tmax_sum_dif_unscaled <= .5 & tmax_sum_dif_unscaled >= -.5)  
   dim(climate_match)
   
   climate_colder <- dat_snp_prs2  %>%
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
     
     
     output_df$`Upper quartile GEBV`[x] <-     mean(dat_snp_prs2$rgr_resids[upper_sample])
     output_df$`Lower quartile GEBV`[x] <-     mean(dat_snp_prs2$rgr_resids[lower_sample])
     output_df$None[x] <-                      mean(dat_snp_prs2$rgr_resids[avg_sample])
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
  
  
   
# Plotting current distribution of prs  --------------------------

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
  
  
  
  ## Remove cells that are out of climate range
  
  # # # Calculate ranges
  # tmax_sum_range <- range(unlist(raster::extract(tmax_rast, lobata_range_rough)), na.rm = TRUE)
  # tmin_winter_range <- range(unlist(raster::extract(tmin_winter, lobata_range_rough)),
  #                            na.rm = TRUE)
  # aet_range <- range(unlist(raster::extract(aet, lobata_range_rough)), na.rm = TRUE)
  # cwd_range <- range(unlist(raster::extract(cwd, lobata_range_rough)), na.rm = TRUE)
  # elevation_range <- range(unlist(raster::extract(dem, lobata_range_rough)), na.rm = TRUE)
  # bioclim_04_range <- range(unlist(raster::extract(bioclim_04, lobata_range_rough)),
  #                           na.rm = TRUE)
  # bioclim_15_range <- range(unlist(raster::extract(bioclim_15, lobata_range_rough)),
  #                           na.rm = TRUE)
  # bioclim_18_range <- range(unlist(raster::extract(bioclim_18, lobata_range_rough)),
  #                           na.rm = TRUE)
  # bioclim_19_range <- range(unlist(raster::extract(bioclim_19, lobata_range_rough)),
  #                           na.rm = TRUE)
  # 
  # 
  # # # Set percentage threshold for range extension of climate and spatial variables
  # range_thresh <- 0.0
  # 
  # gam_dist_rast_df_nona <- gam_dist_rast_df_nona %>%
  #   dplyr::filter(tmax_sum >= tmax_sum_range[1] - (abs(tmax_sum_range[1]) * range_thresh) &
  #                   tmax_sum <= tmax_sum_range[2] + (abs(tmax_sum_range[1]) * range_thresh)) %>%
  #   dplyr::filter(tmin_winter >= tmin_winter_range[1] - (abs(tmin_winter_range[1]) * range_thresh) &
  #                   tmin_winter <= tmin_winter_range[2] + (abs(tmin_winter_range[1]) * range_thresh)) %>%
  #   dplyr::filter(cwd >= cwd_range[1] - (abs(cwd_range[1]) * range_thresh) &
  #                   cwd <= cwd_range[2] + (abs(cwd_range[1]) * range_thresh)) %>%
  #   dplyr::filter(bioclim_04 >= bioclim_04_range[1] - (abs(bioclim_04_range[1]) * range_thresh) &
  #                   bioclim_04 <= bioclim_04_range[2] + (abs(bioclim_04_range[1]) * range_thresh)) %>%
  #   dplyr::filter(bioclim_15 >= bioclim_15_range[1] - (abs(bioclim_15_range[1]) * range_thresh) &
  #                   bioclim_15 <= bioclim_15_range[2] + (abs(bioclim_15_range[1]) * range_thresh)) %>%
  #   dplyr::filter(bioclim_18 >= bioclim_18_range[1] - (abs(bioclim_18_range[1]) * range_thresh) &
  #                   bioclim_18 <= bioclim_18_range[2] + (abs(bioclim_18_range[1]) * range_thresh)) %>%
  #   dplyr::filter(elevation >= elevation_range[1] - (abs(elevation_range[1]) * range_thresh) &
  #                   elevation <= elevation_range[2] + (abs(elevation_range[1]) * range_thresh)) %>%
  #   dplyr::filter(longitude >= extent(lobata_range_rough)@xmin - (abs(extent(lobata_range_rough)@xmin) * range_thresh) &
  #                   longitude <= extent(lobata_range_rough)@xmax + (abs(extent(lobata_range_rough)@xmax) * range_thresh)) %>%
  #   dplyr::filter(latitude >= extent(lobata_range_rough)@ymin - (abs(extent(lobata_range_rough)@ymin) * range_thresh) &
  #                   latitude <= extent(lobata_range_rough)@ymax + (abs(extent(lobata_range_rough)@ymax) * range_thresh))
  # 
  # dim(gam_dist_rast_df_nona)
  # 

  
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
  
  
  rda_df2 <- left_join(rda_df, dplyr::select(gen_dat_moms, accession, prs_best_scaled) %>% 
                         distinct(accession, .keep_all = T) %>% 
                         mutate(accession = as.numeric(as.character(accession)))) %>%
     filter(!is.na(accession)) %>%
    filter(!is.na(prs_best_scaled))
  
  pairs.panels(rda_df2[, climate_vars_gam_dist])
  
  car::vif(lm(snp_chr8_30946730 ~ . , 
              data = rda_df2[, c("snp_chr8_30946730", climate_vars_gam_dist)]))
  
  
  rda_df2$locality <- factor(rda_df2$locality)

  gam_dist <- gam(prs_best_scaled ~ s(tmax_sum, bs = "cr") +
                    s(tmin_winter, bs = "cr") +
                    s(cwd, bs = "cr") + 
                   #     s(bioclim_04, bs = "cr") + 
                    s(bioclim_15, bs = "cr") + 
                    s(bioclim_18, bs = "cr") +  
                    s(elevation, bs = "cr") + 
                    s(latitude, bs = "cr") + 
                    s(longitude, bs = "cr"),
                #   te(latitude, longitude, bs = "cr"),
                  #  s(locality, bs = "re"),
                  #  method = "fREML",
                  # discrete = TRUE,
                  #  select = TRUE,
                  data = rda_df2)
  
  summary(gam_dist)
  
  #   ## Save output
  #   
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
  names(rast) <- snp_name # Name the raster
  
  
  prs_rast_gam <- mask(rast, lobata_range_extended)
  plot(prs_rast_gam)
  
  hist(prs_rast_gam, breaks = 50)
  

  
  # Choose climate variables
  climate_vars_core <- climate_vars[!grepl(pattern = "random|PC", x = climate_vars)]
  climate_vars_core
  length(climate_vars_core) # Should be 10

  ## Calculate number of top and bottom SNPs for moms
  dat_snp_prs_mom <- gen_dat_clim_all %>%
    dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality), 
                     by = c("id" = "gbs_name")) %>%
    dplyr::select(id, accession, locality, longitude, latitude, climate_vars_core) %>%
    dplyr::filter(!id %in% exclude_based_on_missing) %>%
    dplyr::mutate(accession = as.numeric(as.character(accession))) %>% 
    dplyr::left_join(., dplyr::select(gen_dat_moms, id, prs_best_scaled))
    
    head(dat_snp_prs_mom)
      
    summary(dat_snp_prs_mom$prs_best_scaled)
    hist(dat_snp_prs_mom$prs_best_scaled, breaks = 50)

 
    
  ## Exclude localities not in main dataset
   dat_snp_prs_mom <-  dat_snp_prs_mom %>%
      dplyr::filter(locality %in% dat_snp_prs$locality)
   
   ## Exclude accessions not in main dataset
   dat_snp_prs_mom <-  dat_snp_prs_mom %>%
     dplyr::filter(!is.na(accession))
   
  # Average by locality
  dat_snp_prs_locality_mean <- dat_snp_prs_mom %>%
                          dplyr::group_by(locality) %>% 
                          summarise_all(mean)
  
  
 # Set up spatial dataframe
  dat_snp_prs_locality_sp <- SpatialPointsDataFrame(dat_snp_prs_locality_mean[, c("longitude", "latitude")],
                                      data = dat_snp_prs_locality_mean[, c("locality", "prs_best_scaled", "latitude", "longitude")],
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  
  
# Make spatial predictions across landscape for number of alleles vs growth ------------------
  ## Recreating Wang et al. 2010 graphs
    
    
    ## Calculate number of top and bottom SNPs for moms
    dat_snp_prs_mom <- gen_dat_clim_all %>%
      dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality), 
                       by = c("id" = "gbs_name")) %>%
      dplyr::select(id, accession, locality, longitude, latitude, climate_vars_core) %>%
      dplyr::filter(!id %in% exclude_based_on_missing) %>%
      dplyr::mutate(accession = as.numeric(as.character(accession))) %>% 
      dplyr::left_join(., dplyr::select(gen_dat_moms, id, prs_best_scaled))
    
    head(dat_snp_prs_mom)
    
    summary(dat_snp_prs_mom$prs_best_scaled)
    hist(dat_snp_prs_mom$prs_best_scaled, breaks = 50)
    
    # Average by locality
    dat_snp_prs_locality_mean <- dat_snp_prs_mom %>%
      dplyr::group_by(locality) %>% 
      summarise_all(mean)
 
    # Set up dataframe   
    dat_snp_prs_locality_sp <- SpatialPointsDataFrame(dat_snp_prs_locality_mean[, c("longitude", "latitude")],
                                                          data = dat_snp_prs_locality_mean[, c("locality", "prs_best_scaled", "latitude", "longitude")],
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
   
    
  ## Resample PRS Raster to match with tmax rasters

      # From gam predictions
      prs_rast_mean <- resample(prs_rast_gam, tmax_rast) # From gam predictions
      
      compareRaster(tmax_rast, prs_rast_mean) # Make sure rasters are equal
      
      
      # Set up scenario where we pick the highest PRS score within a XXXX radius of the point

      ## Change to equal area projection and resolution to 1km
        prs_rast_mean2 <- projectRaster(from = prs_rast_mean,
                                        crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m", 
                                        res = 1000) # Resolution of 1km
        prs_rast_mean2
        
        
        buffer = focalWeight(prs_rast_mean2, d =  25000, # 25 km radius
                             type = "circle")
      
        buffer[buffer > 0] <- 1 # Set all non 0 buffer weights to 1
        buffer
        dim(buffer)
        
        # Run focal statistics to find max values within radius
          prs_rast_mean_focal = focal(prs_rast_mean2, w = buffer, fun = max, na.rm = T)
          
          plot( prs_rast_mean_focal)
          
        # Project back to lat long
          prs_rast_mean_focal <- projectRaster(from =  prs_rast_mean_focal, 
                                to = prs_rast_mean)
    
          prs_rast_mean_focal <- mask( prs_rast_mean_focal, lobata_range_extended)
          plot( prs_rast_mean_focal)
        
        # Check to make sure all positive values
          hist( prs_rast_mean_focal -  prs_rast_mean)
          summary(values( prs_rast_mean_focal -  prs_rast_mean))
        
        # Check to make sure max values aren't above
          summary(values( prs_rast_mean))
          summary(values( prs_rast_mean_focal))
          
      # Save raster    
        prs_rast_max <-  prs_rast_mean_focal
  
        
   ## Get null prediction of height with no change in tmax

     future_null <- tmax_rast
     names(future_null) <- "tmax_sum_dif"
     
     # Need to scale tmax when actual transfer distance is 0
     future_null[!is.na(future_null)] <- forward_transform(0, var = "tmax_sum_dif",
                                                           means = scaled_var_means_gbs_all,
                                                           sds = scaled_var_sds_gbs_all)
     
     future_null_stack <- stack(future_null, prs_rast_mean) #
     
     
     names(future_null_stack) <- c("tmax_sum_dif", "prs_best_scaled") # Rename so predict function works
     future_null_stack
     
     future_null_height_prs <- predict(future_null_stack,
                                       model = gam_prs_all,
                                       progress = "text",
                                       type = "response")
     
     future_null_height_prs
     plot(future_null_height_prs)
     summary(future_null_height_prs)
     
     ## Add baseline data
     
     # Make data frame for prediction - one for top and one for bottom
     newdata_rast_const <- data.frame(section_block = "Block1_1",
                                      height_2014 = 0,
                                      accession = "6",
                                      PC1_gen = 0, PC2_gen = 0, 
                                      PC1_clim = 0, PC2_clim = 0,
                                      locality = "FHL")
     
     
     future_null_height_prs <- future_null_height_prs +  predict(future_null_stack,
             model = gam_snp_all,
             const = newdata_rast_const,
             progress = "text",
             type = "response")
     

   # Initialize stack for height changes
       future_stack_height_change_prs <- future_stack_tmax_dif
       future_stack_height_change_prs_selection <- future_stack_tmax_dif
   
  
    # Loop over scenarios
    for(scenario in 1:nlayers(future_stack_tmax_dif)){
      
      cat("Working on scenario:", names(future_stack_tmax_dif)[scenario], "... \n")
      
      ## Use scaled raster
        rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
        rast_temp <- stack(rast_temp, prs_rast_mean)
        
        names(rast_temp) <- c("tmax_sum_dif", "prs_best_scaled") # Rename so predict function works
        rast_temp
        
        # Calculate baseline Rgr
        baseline_temp <- predict(rast_temp,
                                 model = gam_snp_all,
                                 const = newdata_rast_const,
                                 progress = "text",
                                 type = "response")
        
        
        
        rast_temp_height_prs <- predict(rast_temp,
                                        model = gam_prs_all,
                                        progress = "text",
                                        type = "response")
        
        ## Add baseline RGR
        rast_temp_height_prs <- rast_temp_height_prs + baseline_temp
        
      
      ## Use selection raster
        rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
        rast_temp <- stack(rast_temp, prs_rast_max)
        
        names(rast_temp) <- c("tmax_sum_dif", "prs_best_scaled") # Rename so predict function works
        rast_temp
        
        rast_temp_height_prs_selection <- predict(rast_temp,
                                        model = gam_prs_all,
                                        progress = "text",
                                        type = "response")
        # Add baseline RGR
        rast_temp_height_prs_selection <- rast_temp_height_prs_selection + baseline_temp

      # Dif bw top and null 
     rast_temp_height_change_prs <- (rast_temp_height_prs - future_null_height_prs) / abs(future_null_height_prs) * 100
     # rast_temp_height_change_prs <- rast_temp_height_prs - future_null_height_prs 
     #   rast_temp_height_change_prs <- rast_temp_height_prs
      plot(rast_temp_height_change_prs, main = scenario)
      future_stack_height_change_prs[[scenario]] <- rast_temp_height_change_prs
      names(future_stack_height_change_prs)[scenario] <- names(future_stack_tmax_dif)[scenario] 
      
      
      # Selection
   #  rast_temp_height_change_prs_selection <- rast_temp_height_prs_selection - future_null_height_prs
   #   rast_temp_height_change_prs_selection <- rast_temp_height_prs_selection
      rast_temp_height_change_prs_selection <- (rast_temp_height_prs_selection - future_null_height_prs) / abs(future_null_height_prs) * 100
      plot(rast_temp_height_change_prs_selection,
           main = paste0(scenario, " | selection"))
      future_stack_height_change_prs_selection[[scenario]] <- rast_temp_height_change_prs_selection
      names(future_stack_height_change_prs_selection)[scenario] <- names(future_stack_tmax_dif)[scenario] 
      

    } # End scenario loop
     
    
    # Average across scenarios
     future_stack_height_change_prs_avg <- mean(future_stack_height_change_prs, na.rm = TRUE)
     
     future_stack_height_change_prs_selection_avg <- mean(future_stack_height_change_prs_selection,
                                                          na.rm = TRUE)
     
     hist(future_stack_height_change_prs_avg)
     hist(future_stack_height_change_prs_selection_avg)
     
     
     height_change_dif = future_stack_height_change_prs_selection_avg - future_stack_height_change_prs_avg
     plot(height_change_dif)
     hist(height_change_dif)
     
     # Find min across both rasters
     rast_min <- min(c(min(values(future_stack_height_change_prs_avg), na.rm = T), 
                    min(values(future_stack_height_change_prs_selection_avg), na.rm = T)))
     rast_max <- max(c(max(values(future_stack_height_change_prs_avg), na.rm = T), 
                       max(values(future_stack_height_change_prs_selection_avg), na.rm = T)))
     
     
     hill_cropped <- crop(hill, extent(-125, -118, 33, 42.1)) # crop based on california outline
   
     
     
     
     
     
  # Plot of baseline PRS scores  
    # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +

       # Add kriged raster
     #  tm_shape(future_stack_height_change_prs_selection_avg) +
       tm_shape(prs_rast_mean) +
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
       tm_shape(dat_snp_prs_locality_sp) +
       tm_bubbles(size = .025,
                  shape = 21,
                  col = "black",
                  border.alpha = 0,
                  midpoint = NA)
     
     tm

     tmap_save(tm, filename = paste0("./figs_tables/fig3/Figure 3 - prs map ", Sys.Date(), ".pdf"),
               width = 12, height =12, units = "cm", useDingbats=FALSE)
     
    
     
     
      
   ## Plot of predicted growth rates in 2080 based on current PRS distribution
     
   # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +
       
       # Add kriged raster
       tm_shape(future_stack_height_change_prs_avg) +
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
       tm_shape(dat_snp_prs_locality_sp) +
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
     
     
     
     
     ## Plot of predicted growth rates in 2080 based on MAX PRS
     
     # Plot map
     tm <-  tm_shape(hill_cropped)+
       tm_grid(lwd = .025, labels.inside.frame = TRUE) + # Add lat long lines
       tm_xlab("Longitude") + tm_ylab("Latitude") + # Add axis labels
       tm_raster(palette = gray.colors(n= 9), legend.show = FALSE) + # Add hillshade
       tm_shape(cali_outline) +
       tm_borders("black") +
       # tm_fill("white", alpha = .25) +
       
       # Add kriged raster
       tm_shape(future_stack_height_change_prs_selection_avg) +
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
       tm_shape(dat_snp_prs_locality_sp) +
       tm_bubbles(size = .025,
                  shape = 21,
                  col = "black",
                  border.alpha = 0,
                  midpoint = NA)
     tm
     
     
     # tmaptools::palette_explorer()
     
     # Save plot
     tmap_save(tm, filename = paste0("./figs_tables/fig3/Figure 3 - growth PRS selection 2080 map ", Sys.Date(), ".pdf"),
               width = 12, height =12, units = "cm", useDingbats=FALSE)
     
    
    
# RDA  --------------------------------------------------------------------

 library(vegan)
 
 # Count % of missing data for full gbs dataset
 missing_all <- data.frame(gen_dat_clim_all[ "id"],
                           missing = apply(dplyr::select(gen_dat_clim_all, snp_col_names),
                                           MARGIN = 1, function(x) sum(is.na(x)))/length(snp_col_names))
 
 # Get list of samples to exclude because of missing data
 exclude_based_on_missing <- missing_all$id[which(missing_all$missing >= .2)]

 
# Format data 
   rda_df <-  gen_dat_clim %>%
     dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality),
                      by = c("id" = "gbs_name")) %>%
     dplyr::filter(!(id %in% exclude_based_on_missing))
   
   
   rda_df2 <- left_join(rda_df, dplyr::select(gen_dat_moms, accession, prs_best_scaled) %>% 
                          distinct(accession, .keep_all = T) %>% 
                          mutate(accession = as.numeric(as.character(accession)))) %>%
     filter(!is.na(accession)) %>%
     filter(!is.na(prs_best_scaled))
     
   
   head(rda_df2)
   dim(rda_df2)
   
    

   # # Average allele frequency by locality
   # library(multidplyr)
   # rda_df <- rda_df %>%
   #   partition(locality) %>%
   #   summarise_all(mean, na.rm = TRUE) %>%
   #   collect()
   # 
   # dim(rda_df)
   # 
   # # Testing removing laytonville
   # # rda_df <- rda_df %>%
   # #   filter(locality != "LAY")
   # 
   # ## Remove localities with less than ## individuals
   # count_by_locality <- gen_dat_clim_all %>%
   #   dplyr::select(id) %>%
   #   dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality),
   #                    by = c("id" = "gbs_name")) %>%
   #   dplyr::select(id, locality) %>%
   #   dplyr::filter(!id %in% exclude_based_on_missing) %>%
   #   group_by(locality) %>%
   #   tally()
   # 
   # table(count_by_locality$n)
   # 
   # exclude_locality <- count_by_locality$locality[which(count_by_locality$n < 2)]
   # rda_df <- rda_df %>%
   #   dplyr::filter(!locality %in% exclude_locality)
   # 
   # dim(rda_df)
   # 
   # ## Remove snps based on # of unique values
   # n_distinct <- dplyr::select(rda_df, snp_col_names) %>%
   #   summarise_all(n_distinct)
   # n_distinct <- unlist(n_distinct)
   # 
   # exlude_by_unique_count <- names(which(n_distinct <= 5))
   # exlude_by_unique_count
   # rda_df <- rda_df %>%
   #   dplyr::select(-exlude_by_unique_count)


 # Choose climate variables
   climate_vars_core <- climate_vars[!grepl(pattern = "random|PC", x = climate_vars)]
   climate_vars_core
   length(climate_vars_core) # Should be 10
   
   wna_climate_vars <- climate_vars_core
   
   wna_climate_vars_current <- paste0(wna_climate_vars, "_wna_current")
   wna_climate_vars_lgm <- paste0(wna_climate_vars, "_lgm")
   wna_climate_vars_rcp85 <- paste0(wna_climate_vars, "_rcp85")
   
   
 ## Calculate MEMs
   
   # Calculate distance matrix
   library(geosphere)
   
   mom_lat_longs <- as.data.frame(dplyr::select(rda_df, longitude, latitude))
   dist_mat_gf <- matrix(NA, ncol = nrow(rda_df), nrow = nrow(rda_df))
   
   for(i in 1:nrow(dist_mat_gf)){
     
     if(i %% 25 == 0){
       cat("Working on row:", i, " ... \n")
     }
     
     for(j in 1:nrow(dist_mat_gf)){
       
       dist_mat_gf[i, j] <- distHaversine(mom_lat_longs[i, ], 
                                          mom_lat_longs[j, ])/ 1000 # To convert to km
       
     }
   }
   
   head(dist_mat_gf)
   
   gf_mems = vegan::pcnm(dist_mat_gf)  
   gf_mems$values
   
   # Bind eigenvectors
   rda_df <- bind_cols(rda_df, as.data.frame(gf_mems$vectors[, 1:10]))
   
   mems <- c("PCNM1", "PCNM2", "PCNM3", "PCNM4", "PCNM5",
             "PCNM6", "PCNM7", "PCNM8", "PCNM9", "PCNM10")
   
   
   # Make separate climate datasets
   clim_dat_lgm <- rda_df[, c(wna_climate_vars_lgm)]
   colnames(clim_dat_lgm)
   
   clim_dat_current <- rda_df[, c(wna_climate_vars_current)]
   colnames(clim_dat_current)
   
   clim_dat_rcp85 <- rda_df[, c(wna_climate_vars_rcp85)]
   colnames(clim_dat_rcp85)
   
   spatial_vars_rda <- names(dplyr::select(rda_df, PCNM1:PCNM10)) # 10
   
   
   
   # Scale climate and spatial variables
   rda_all_covariates_scaled <- rda_df %>%
     dplyr::select(-contains("snp_chr")) %>%
     mutate_if(is.numeric, scale)
   
   
   # Scale genotype data
   rda_all_scaled <- flashpcaR::scale2(dplyr::select(rda_df, contains("snp_chr")), impute = 2)
   dim(rda_all_scaled)
   
   rda_top_scaled <- flashpcaR::scale2(dplyr::select(rda_df, top_snps_long$snp), impute = 2)
   dim(rda_top_scaled) 
   
   rda_bottom_scaled <- flashpcaR::scale2(dplyr::select(rda_df, bottom_snps_long$snp), impute = 2)
   dim(rda_bottom_scaled)
   
   rda_outlier_scaled <- flashpcaR::scale2(dplyr::select(rda_df, unique(c(bottom_snps_long$snp,
                                                                top_snps_long$snp))), impute = 2)
   dim(rda_outlier_scaled)
   
   rda_nonoutlier_scaled <- flashpcaR::scale2(dplyr::select(rda_df, snp_col_names[!snp_col_names %in%
                                                            unique(c(bottom_snps_long$snp,
                                                                          top_snps_long$snp))]), impute = 2)
   dim(rda_nonoutlier_scaled)
   
  
   climate_vars_rda_current <- paste0(wna_climate_vars, "_wna_current")
   climate_vars_rda_lgm <- paste0(wna_climate_vars, "_lgm")
   
   ## Testing correlations with BCM and climateWNA datasets
   for(var in wna_climate_vars){
    test = cor.test(pull(gen_dat_clim_all, var),
                   pull(gen_dat_clim_all, paste0(var, "_wna_current")))
    cat(var, ": ", test$estimate, " ... \n")
    
   }
   
   
   
   # Select same number of PCNMs as climate variables
   spatial_vars_rda <- names(dplyr::select(rda_df, PCNM1:PCNM10)) # 10
    

   library(colorfulVennPlot)

   ## Function to calculate variance partitioning
    calc_varpart <- function(genetic_data, label, permtest = FALSE){

      # With all data
      rda_out_full <- rda(as.formula(paste("genetic_data ~ ",
                                           paste(c(spatial_vars_rda,
                                                   climate_vars_rda_current,
                                                   climate_vars_rda_lgm), collapse = "+"))),
                          data = rda_all_covariates_scaled,
                          scale = FALSE)
      
      if(permtest == TRUE){
        cat("RDA out full : ")
       print(anova(rda_out_full, permutations = 99))
       "_____________________________ \n"
      }

      # Partial out spatial component
      rda_out_part_spat <- rda(as.formula(paste("genetic_data ~ ",
                                           paste(c(climate_vars_rda_current,
                                                   climate_vars_rda_lgm), collapse = "+"),
                                           "+Condition(",
                                           paste(spatial_vars_rda, collapse = "+"),")")),
                          data = rda_all_covariates_scaled,
                          scale = FALSE)
      
      if(permtest == TRUE){
        cat("rda_out_part_spat : ")
        print(anova(rda_out_part_spat, permutations = 99))
        "_____________________________ \n"
      }

      # Partial out current climate component
      rda_out_part_clim_current <- rda(as.formula(paste("genetic_data ~ ",
                                                paste(c(spatial_vars_rda,
                                                climate_vars_rda_lgm), collapse = "+"),
                                                "+Condition(",
                                        paste(climate_vars_rda_current, collapse = "+"),")")),
                                       data = rda_all_covariates_scaled,
                                       scale = FALSE)
      
      if(permtest == TRUE){
        cat("rda_out_part_clim_current : ")
        print(anova(rda_out_part_clim_current, permutations = 99))
        "_____________________________ \n"
      }

      # Partial out lgm climate component
      rda_out_part_clim_lgm <- rda(as.formula(paste("genetic_data ~ ",
                                        paste(c(spatial_vars_rda,
                                                climate_vars_rda_current), collapse = "+"),
                                        "+Condition(",
                                        paste(climate_vars_rda_lgm, collapse = "+"),")")),
                                   data = rda_all_covariates_scaled,
                                   scale = FALSE)
      
      if(permtest == TRUE){
        cat("rda_out_part_clim_lgm : ")
        print(anova(rda_out_part_clim_lgm, permutations = 99))
        "_____________________________ \n"
      }

      # Pure Spatial - partial out climate
      rda_out_pure_spat <- rda(as.formula(paste("genetic_data ~ ",
                                paste(spatial_vars_rda, collapse = "+"),
                                "+Condition(",
                                paste(c(climate_vars_rda_current, climate_vars_rda_lgm),
                                        collapse = "+"),")")),
                               data = rda_all_covariates_scaled,
                               scale = FALSE)
      
      if(permtest == TRUE){
        cat(" rda_out_pure_spat : ")
        print(anova( rda_out_pure_spat , permutations = 99))
        "_____________________________ \n"
      }
      
      # Pure current climate
      rda_out_pure_clim_current <-  rda_out_pure_spat <- rda(as.formula(paste("genetic_data ~ ",
                                      paste(climate_vars_rda_current, collapse = "+"),
                                      "+Condition(",
                                      paste(c(spatial_vars_rda, climate_vars_rda_lgm),
                                            collapse = "+"),")")),
                                                             data = rda_all_covariates_scaled,
                                                             scale = FALSE)
      
      if(permtest == TRUE){
        cat(" rda_out_pure_clim_current: ")
        print(anova( rda_out_pure_clim_current, permutations = 99))
        "_____________________________ \n"
      }

      # Pure lgm climate
      rda_out_pure_clim_lgm <-  rda_out_pure_spat <- rda(as.formula(paste("genetic_data ~ ",
                                      paste(climate_vars_rda_lgm, collapse = "+"),
                                      "+Condition(",
                                      paste(c(spatial_vars_rda, climate_vars_rda_current),
                                            collapse = "+"),")")),
                                                         data = rda_all_covariates_scaled,
                                                         scale = FALSE)
      
      if(permtest == TRUE){
        cat(" rda_out_pure_clim_lgm : ")
        print(anova( rda_out_pure_clim_lgm, permutations = 99))
        "_____________________________ \n"
      }

      # Based on total variance (included unexplained variance)
   #  overall_iner <- summary(rda_out_full)$tot.chi
      
      # Based on explained variance
      overall_iner <- summary(rda_out_full)$constr.chi

      ## Set up venn diagram
# fit1 <- euler(c("Spatial" = round(summary(rda_out_pure_spat)$constr.chi/overall_iner*100, 1),
# "LGM climate" =   round(summary(rda_out_pure_clim_lgm)$constr.chi/overall_iner*100, 1),
# "Current climate" = round(summary(rda_out_pure_clim_current)$constr.chi/overall_iner*100, 1),
# "Spatial&Current climate" = round(summary(rda_out_part_clim_lgm)$constr.chi/overall_iner*100, 1),
# "Spatial&LGM climate" = round(summary(rda_out_part_clim_current)$constr.chi/overall_iner*100, 1),
# "Current climate&LGM climate" = round(summary(rda_out_part_spat)$constr.chi/overall_iner*100, 1),
# "Spatial&Current climate&LGM climate" = round(summary(rda_out_full)$constr.chi/overall_iner*100, 1)),
# shape = "circle")
#
#       plot(fit1,
#            quantities = TRUE,
#          #  fills = list(fill = c("red", "steelblue", "green")),
#            main = label)
#         #  legend = list(labels = c("Spatial", "Current climate", "LGM climate")))
#
      # Equal area venn diagram
      y <-  c("100" = round(summary(rda_out_pure_spat)$constr.chi/overall_iner*100, 2),
              "001" = round(summary(rda_out_pure_clim_lgm)$constr.chi/overall_iner*100, 2),
              "010" = round(summary(rda_out_pure_clim_current)$constr.chi/overall_iner*100, 2),
              "110" = round(summary(rda_out_part_clim_lgm)$constr.chi/overall_iner*100, 2),
              "101" = round(summary(rda_out_part_clim_current)$constr.chi/overall_iner*100, 2),
              "011" = round(summary(rda_out_part_spat)$constr.chi/overall_iner*100, 2),
              "111" = round(summary(rda_out_full)$constr.chi/overall_iner*100, 2))

     # plot.new()
      graphics.off()
     # plot.new()
      plotVenn3d(y, labels = c("Spatial", "Current climate", "LGM climate"),
                 Title = label,
                 Colors = c("#99cc66","#9999cc","#ff9999",
                            "#ffff99","#99ccff","#cc99cc","#ffffff"),
                 rot = 180)


    ## Just current climate and spatial

      # y <-  c("10" = round(summary(rda_out_pure_spat)$constr.chi/overall_iner*100, 1),
      #         "01" =   round(summary(rda_out_pure_clim_current)$constr.chi/overall_iner*100, 1),
      #         "11" = round(summary(rda_out_full)$constr.chi/overall_iner*100, 1))
      # 
      # plot.new()
      # # plot.new()
      # plotVenn2d(y, labels = c("Spatial", "Current climate"), Colors = "white")
    }

  # Calculate variance partitions
    calc_varpart(genetic_data = rda_top_scaled, label = "Warm advantageous")
    calc_varpart(genetic_data = rda_bottom_scaled, label = "Warm disadvantageous")
    calc_varpart(genetic_data = rda_outlier_scaled, label = "Outlier SNPs")
    
    calc_varpart(genetic_data = rda_nonoutlier_scaled, label = "Non-Outlier SNPs")
    
  
  # For all Snps
     calc_varpart(genetic_data = rda_all_scaled, label = "All SNPs", permtest = TRUE)
    # dev.copy(pdf, paste0("./figs_tables/Figure S3 - venn diagram ", Sys.Date(), ".pdf"),
    #          width = 3, height = 3, useDingbats = FALSE)
    # dev.off()
     
     # TO get overall r2
     rda_out_full <- rda(as.formula(paste("rda_all_scaled ~ ",
                                          paste(c(spatial_vars_rda,
                                                  climate_vars_rda_current,
                                                  climate_vars_rda_lgm), collapse = "+"))),
                         data = rda_all_covariates_scaled,
                         scale = FALSE)
     
     rda_out_full
     
     

  #  RDA on LGM climate
    rda_all_mod <- rda(as.formula(paste("rda_all_scaled ~ ",
                                        paste(c(wna_climate_vars_lgm), collapse = "+"),
                                        "+Condition(",
                                        paste(c(spatial_vars_rda, wna_climate_vars_current), 
                                              collapse = "+"),")")),
                       data = rda_all_covariates_scaled,
                       scale = FALSE)
    
    rda_out_full <- rda(as.formula(paste("genetic_data ~ ",
                                         paste(c(spatial_vars_rda,
                                                 climate_vars_rda_current,
                                                 climate_vars_rda_lgm), collapse = "+"))),
                        data = rda_all_covariates_scaled,
                        scale = FALSE)rda_all_mod
    
    # anova(rda_all_mod, permutations = 99)
    # anova(rda_all_mod, permutations = 99, by = "axis")
    # anova(rda_all_mod, permutations = 99, by = "terms")
    

    
 ## Biplot highlighting outlier SNPs
    
  # Outliers function from Forester et al.
    outliers <- function(loadings, sd_cutoff){
      # find loadings +/-z sd from mean loading
      lims <- mean(loadings) + c(-1, 1) * sd_cutoff * sd(loadings)
      loadings[loadings < lims[1] | loadings > lims[2]] # locus names in these tails
    }


# Identify outliers
 loadings <- summary(rda_all_mod)$species # Extracts loadings for each SNP (named species here)
 head(loadings)
 
 cand_rda1 <- outliers(loadings[, 1], sd_cutoff = 2.5)  # 2.5 used in Forrester et al. 2018
 cand_rda1
 cand_rda2 <- outliers(loadings[, 2], sd_cutoff = 2.5)
 cand_rda2
 outlier_snps_names <- c(names(cand_rda1), names(cand_rda2))
 length(outlier_snps_names)
 
 # Assign arrows to colors
 
 arrow_names <- names(rda_all_mod$terminfo$ordered)
 # Remove PCNM if partialling
 arrow_names <- arrow_names[!grepl(pattern = "PCNM", arrow_names)]
 arrow_cols <- arrow_names
 arrow_cols[grepl("PCNM", arrow_names)] <- "black"
 # temp variables
 arrow_cols[grepl("tave", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("MWMT", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("MCMT", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("TD", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("DD_0", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("DD5", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("DD_18", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("DD18", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("NFFD", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("bFFP", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("eFFP", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("FFP", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("PAS", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("EMT", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("EXT", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("tmax_sum", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("tmax", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("tmin", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("tmin_winter", arrow_names)] <- "#D86F27"
 arrow_cols[grepl("bioclim_04", arrow_names)] <- "#D86F27"
 
# Precipitation variables
 arrow_cols[grepl("MAP", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("cwd", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("bioclim_15", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("bioclim_18", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("bioclim_19", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("MSP", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("AHM", arrow_names)] <- "#65BADA" # double check
 arrow_cols[grepl("SHM", arrow_names)] <- "#65BADA" # double check
 arrow_cols[grepl("PAS", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("Eref", arrow_names)] <- "#65BADA" # double check
 arrow_cols[grepl("CMD", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("RH", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("PPT_sm", arrow_names)] <- "#65BADA"
 arrow_cols[grepl("PPT_wt", arrow_names)] <- "#65BADA"
 
 arrow_cols
 
 
 ## How many RDA outliers also snp outliers?
   sum(outlier_snps_names %in% c(top_snps_long$snp, bottom_snps_long$snp))
   outlier_snps_names[outlier_snps_names %in% c(top_snps_long$snp, bottom_snps_long$snp)]
   
    double_outliers <- outlier_snps_names[outlier_snps_names %in% c(top_snps_long$snp, bottom_snps_long$snp)]
    double_outliers
    
    top_double <- outlier_snps_names[outlier_snps_names %in% c(top_snps_long$snp)]
    top_double
    
    bottom_double <- outlier_snps_names[outlier_snps_names %in% c(bottom_snps_long$snp)]
    bottom_double
    
    # Make a contingency table
    fisher.test(matrix(c(length(double_outliers),
                         length(outlier_snps_names) - length(double_outliers),
                         length(unique(c(top_snps_long$snp, bottom_snps_long$snp))) - length(double_outliers) ,
                         length(rda_all_mod$colsum) -length(double_outliers) -  (length(outlier_snps_names) - length(double_outliers)) - (length(unique(c(top_snps_long$snp, bottom_snps_long$snp))) - length(double_outliers))), ncol = 2, byrow = T))
    
    
    # Top
    fisher.test(matrix(c(length(top_double),
                         length(outlier_snps_names) - length(top_double),
                         length(unique(c(top_snps_long$snp, bottom_snps_long$snp))) - length(top_double) ,
                         length(rda_all_mod$colsum) -length(top_double) -  (length(outlier_snps_names) - length(top_double)) - (length(unique(c(top_snps_long$snp, bottom_snps_long$snp))) - length(top_double))),
                       ncol = 2, byrow = T))
    
    # Bottom
    fisher.test(matrix(c(length(bottom_double),
                         length(outlier_snps_names) - length(bottom_double),
                         length(unique(c(bottom_snps_long$snp, bottom_snps_long$snp))) - length(bottom_double) ,
                         length(rda_all_mod$colsum) -length(bottom_double) -  (length(outlier_snps_names) - length(bottom_double)) - (length(unique(c(bottom_snps_long$snp, bottom_snps_long$snp))) - length(bottom_double))),
                       ncol = 2, byrow = T))
    

 # Labels of outlier SNPs
   snp_labels <- rep("", times = nrow(loadings))
   snp_labels[rownames(loadings) %in% double_outliers] <- rownames(loadings)[rownames(loadings) %in% double_outliers]
   snp_labels <- gsub(x = snp_labels, pattern = "snp_", replacement = "")
   table(snp_labels)
 
   outlier_snps_names[outlier_snps_names %in% c(top_snps_long$snp)] # Both climate outliers in top snps
   outlier_snps_names[outlier_snps_names %in% c(bottom_snps_long$snp)]
 
 ## Figure out what climate variable has the strongest correlation
   
   sort(abs(summary(rda_all_mod)$biplot[, "RDA1"]), decreasing = T)
   sort(abs(summary(rda_all_mod)$biplot[, "RDA2"]), decreasing = T)
   
   rda_outlier_cors <- data.frame(snp = rep(NA, length(double_outliers)),
                                  rda_axis = rep(NA, length(double_outliers)),
                                  rda_cor = rep(NA, length(double_outliers)),
                                  topvar1 = rep(NA, length(double_outliers)),
                                  topvar1_dir = rep(NA, length(double_outliers)))
   x = 1
   
   for(snp in double_outliers){
     
     rda_outlier_cors$snp[x] <- snp
     
     # Check if in RDA1
     if(snp %in% names(cand_rda1)){
       rda_outlier_cors$rda_axis[x] <- "RDA1"
       rda_outlier_cors$rda_cor[x] <- cand_rda1[names(cand_rda1) == snp]
      # Find strongest associations - could speed up by calculate summary outtside function
       rda_outlier_cors$topvar1[x] <- names(sort(abs(summary(rda_all_mod)$biplot[, "RDA1"]), decreasing = T))[1]
       rda_outlier_cors$topvar1_dir[x] <- summary(rda_all_mod)$biplot[, "RDA1"][names(summary(rda_all_mod)$biplot[, "RDA1"]) ==  rda_outlier_cors$topvar1[x]]
     }
     
     # Check if in RDA2
     if(snp %in% names(cand_rda2)){
       rda_outlier_cors$rda_axis[x] <- "RDA2"
       rda_outlier_cors$rda_cor[x] <- cand_rda2[names(cand_rda2) == snp]
       # Find strongest associations - could speed up by calculate summary outtside function
       rda_outlier_cors$topvar1[x] <- names(sort(abs(summary(rda_all_mod)$biplot[, "RDA2"]), decreasing = T))[1]
       rda_outlier_cors$topvar1_dir[x] <- summary(rda_all_mod)$biplot[, "RDA2"][names(summary(rda_all_mod)$biplot[, "RDA2"]) ==  rda_outlier_cors$topvar1[x]]
     }
     
     # Plot relationship
     # plot(c(pull(rda_all_covariates_scaled, rda_outlier_cors$topvar1[x])), 
     #      rda_all_scaled[, snp], xlab = rda_outlier_cors$topvar1[x], ylab = snp, main = snp, pch = 19, col = "dodgerblue")
     # abline(lm(rda_all_scaled[, snp] ~ c(pull(rda_all_covariates_scaled, rda_outlier_cors$topvar1[x]))), lwd = 2) 
     
     x = x +1
   }
   
   rda_outlier_cors
   
# Make biplot
 cex_points = 1.1
 axis_lim <- .5
 plot(rda_all_mod, type = "n", scaling = 3, 
      xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim), las = 1) # Create empty plot 
 
 # All snps
 points(rda_all_mod, display = "species", pch = 21,
        col = "black", bg = "grey95", scaling = 3, cex = cex_points)
 
 # Outlier snps -rda
 points(rda_all_mod, display = "species", pch = 21, # Outliers
        col = "black", bg = "grey50", scaling = 3, cex = cex_points,
        select = rownames(loadings) %in% outlier_snps_names)
 # Bottom snps
 points(rda_all_mod, display = "species", pch = 21,
        col = "black", bg = "#af8dc3", scaling = 3, cex = cex_points,
        select = rownames(loadings) %in% bottom_snps_long$snp)
 
 # Top snps
 points(rda_all_mod, display = "species", pch = 21, cex = cex_points,
        col = "black", bg = "#7fbf7b", scaling = 3,
        select = rownames(loadings) %in% top_snps_long$snp)
  
 # Add double outliers - growth and climate
  points(rda_all_mod, display = "species", pch = 23, cex = cex_points * 1.5,
         col = "black", bg = "#7fbf7b", scaling = 3,
         select = rownames(loadings) %in% top_double)
  points(rda_all_mod, display = "species", pch = 23, cex = cex_points * 1.5,
         col = "black", bg = "#af8dc3", scaling = 3,
         select = rownames(loadings) %in% bottom_double)
  #text(rda_all_mod, display = "species", labels = snp_labels, cex = .75, scaling = 3)
  
 # Environmental vectors
  text(rda_all_mod, scaling = 3, display = "bp",
       col= arrow_cols, lwd = 3,
       cex = 1) 
  
  # Save as pdf
  dev.copy(pdf, paste0("./figs_tables/fig3/Figure 3 - rda biplot ", Sys.Date(), ".pdf"),
           width = 5, height = 5, useDingbats = FALSE)
  dev.off()

  
  ## Calculate sum of loadings for outliers and non outliers
  
  loadings_df <- as.data.frame(loadings)
  loadings_df$snp <- rownames(loadings_df)
  head(loadings_df)
  
  top_temp <- data.frame(snp = top_snps_long$snp,
                         outlier_type = "Warm advantageous")
  bottom_temp <- data.frame(snp = bottom_snps_long$snp,
                         outlier_type = "Warm disadvantageous")
  all_temp <- data.frame(snp = snp_col_names[!snp_col_names %in% c(top_snps_long$snp,
                                                                   bottom_snps_long$snp)],
                         outlier_type = "non-outlier")
  
  rda_loadings <- rbind(top_temp, bottom_temp, all_temp) %>%
                  dplyr::left_join(., loadings_df)
  
  head(rda_loadings)
  dim(rda_loadings)
  
  ggplot(rda_loadings, aes(x = outlier_type, y = abs(RDA1))) + 
    geom_boxplot()
  
  rda_loadings %>%
    dplyr::select(-snp) %>%
    group_by(outlier_type) %>%
    mutate(RDA1 = abs(RDA1),
           RDA2 = abs(RDA2),
           RDA3 = abs(RDA3)) %>%
    summarise_all(median, na.rm = TRUE)

  
  

# Gradient Forest Test - the return! --------------------------------------


  # Growth predictions by GBS mom
  ## Scale number of top and bottom snps - careful about running this twice!!
  dat_snp_prs_mom_unscaled <- dat_snp_prs_mom
  
  
  # Set up gradient forest

  library(gradientForest)

  
  climate_vars_core <- climate_vars[!grepl(pattern = "random|PC", x = climate_vars)]
  
  climate_vars_core
#  wna_climate_vars <- climate_vars_core[!climate_vars_core %in% c("tave", "tmax", "tmin")]
  
  
  wna_climate_vars_current <- paste0(wna_climate_vars, "_wna_current")
  wna_climate_vars_lgm <- paste0(wna_climate_vars, "_lgm")
  wna_climate_vars_rcp85 <- paste0(wna_climate_vars, "_rcp85")
  
  pairs.panels(gen_dat_clim_all[, wna_climate_vars_lgm ])
  
 
  ## Select genotype data
  # Count % of missing data for full gbs dataset
  missing_all <- data.frame(gen_dat_clim_all[ "id"],
                            missing = apply(dplyr::select(gen_dat_clim_all, snp_col_names),
                                            MARGIN = 1, function(x) sum(is.na(x)))/length(snp_col_names))
  
  # Get list of samples to exclude because of missing data
  exclude_based_on_missing <- missing_all$id[which(missing_all$missing > .2)]
  
  ## Set up df
  gen_dat_gf <- gen_dat_clim_all %>%
    dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality), 
                     by = c("id" = "gbs_name")) %>%
    dplyr::select(id, accession, locality, snp_col_names, wna_climate_vars_current,
                  wna_climate_vars_lgm, wna_climate_vars_rcp85, latitude, longitude) %>%
    dplyr::filter(!id %in% exclude_based_on_missing)#  %>%
  #  dplyr::filter(!is.na(accession))
  
  head(gen_dat_gf)
  dim(gen_dat_gf)
  
  # Average allele frequency by locality
  #   library(multidplyr)
  #   gen_dat_gf <- gen_dat_gf %>%
  #     partition(locality) %>%
  #     summarise_all(mean, na.rm = TRUE) %>%
  #     collect()
  # 
  #   dim(gen_dat_gf)
  # 
  # ## Remove localities with less than ## individuals
  #   count_by_locality <- gen_dat_clim_all %>%
  #     dplyr::select(id) %>%
  #     dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality),
  #                      by = c("id" = "gbs_name")) %>%
  #     dplyr::select(id, locality) %>%
  #     dplyr::filter(!id %in% exclude_based_on_missing) %>%
  #     group_by(locality) %>%
  #     tally()
  # 
  #   table(count_by_locality$n)
  # 
  #   exclude_locality <- count_by_locality$locality[which(count_by_locality$n < 2)]
  #   gen_dat_gf <- gen_dat_gf %>%
  #     dplyr::filter(!locality %in% exclude_locality)
  # 
  #   dim(gen_dat_gf)
  # 
  # # Remove snps based on # of unique values
  #   n_distinct <- dplyr::select(gen_dat_gf, snp_col_names) %>%
  #     summarise_all(n_distinct)
  #   n_distinct <- unlist(n_distinct)
  # 
  #   exlude_by_unique_count <- names(which(n_distinct <= 5))
  #   exlude_by_unique_count
  #   gen_dat_gf <- gen_dat_gf %>%
  #   dplyr::select(-exlude_by_unique_count)
  # 
  # 
  ## Calculate MEMs
  
  # Calculate distance matrix
  library(geosphere)
  
  mom_lat_longs <- as.data.frame(dplyr::select(gen_dat_gf, longitude, latitude))
  dist_mat_gf <- matrix(NA, ncol = nrow(gen_dat_gf), nrow = nrow(gen_dat_gf))
  
  for(i in 1:nrow(dist_mat_gf)){
    
    if(i %% 25 == 0){
      cat("Working on row:", i, " ... \n")
    }
    
    for(j in 1:nrow(dist_mat_gf)){
      
      dist_mat_gf[i, j] <- distHaversine(mom_lat_longs[i, ], 
                                      mom_lat_longs[j, ])/ 1000 # To convert to km
      
    }
  }
  
  head(dist_mat_gf)

  gf_mems = vegan::pcnm(dist_mat_gf)  
  gf_mems$values
  
  # Bind eigenvectors
  gen_dat_gf <- bind_cols(gen_dat_gf, as.data.frame(gf_mems$vectors[, 1:10]))
  
  mems <- c("PCNM1", "PCNM2", "PCNM3", "PCNM4", "PCNM5", "PCNM6", "PCNM7")


  # Make separate climate datasets
  clim_dat_lgm <- gen_dat_gf[, c(wna_climate_vars_lgm)]
  colnames(clim_dat_lgm) <- gsub(pattern = "_lgm", "", colnames(clim_dat_lgm))
  colnames(clim_dat_lgm)
  
  clim_dat_current <- gen_dat_gf[, c(wna_climate_vars_current)]
  colnames(clim_dat_current) <- gsub(pattern = "_wna_current", "", 
                                     colnames(clim_dat_current))
  colnames(clim_dat_current)
  
  clim_dat_rcp85 <- gen_dat_gf[, c(wna_climate_vars_rcp85)]
  colnames(clim_dat_rcp85) <- gsub(pattern = "_rcp85", "", 
                                     colnames(clim_dat_rcp85))
  colnames(clim_dat_rcp85)

  pairs.panels(clim_dat_lgm)

  ## Factors
    library(tidyimpute)
    gf_top <- dplyr::select(gen_dat_gf, contains("snp_"))
  #  gf_top <- gf_top[, !colnames(gf_top) %in% unique(c(bottom_snps_long$snp, top_snps_long$snp))]
   # gf_top <- dplyr::select(gen_dat_gf, bottom_snps_long$snp)
  #  gf_top <- dplyr::select(gen_dat_gf, unique(c(bottom_snps_long$snp, top_snps_long$snp)))
   # gf_top_imputed <- flashpcaR::scale2(gf_top, impute = 2)
    gf_top_imputed <- tidyimpute::impute_mean(gf_top) # Mean impute
    summary(gf_top_imputed[, 1:20])

    maxLevel <- floor(log2(0.368*nrow(gf_top)/2)) #account for correlations, see ?gradientForest
    
    # Randomize climate data
    # for(col in 1:ncol(clim_dat_lgm)){
    #   clim_dat_lgm[, col] <-  clim_dat_lgm[sample(nrow(clim_dat_lgm)), col]
    # }
    

   # Run gradient Forest
           gf <- gradientForest(cbind(clim_dat_lgm,
                                         gf_top_imputed),
                             predictor.vars = colnames(clim_dat_lgm),
                            response.vars = colnames(gf_top_imputed),
                            trace = TRUE,
                           ntree = 50,
                        #   nbin = 201,
                           maxLevel = maxLevel, # From Fitzpatrick R script
                           corr.threshold = 0.5)

        gf$result # Rsquared of positive loci
        gf$species.pos.rsq

        cat("\n", mean(gf$result), "\n")
      
        importance(gf)

    plot(gf, plot.type = "O")
    
    
    # Calculate and map "genetic offset" under climate change ----------------------
    # Script assumes:
    # (1) a dataframe of transformed env. variables for CURRENT climate
    # (e.g., predGI5 from above).
    #
    # (2) a dataframe named env_trns_future containing extracted raster data of
    # env. variables for FUTURE a climate scenario, same structure as env_trns
    
    # first transform FUTURE env. variables
    clim_dat_lgm <- as.data.frame(clim_dat_lgm)
    clim_dat_current <- as.data.frame(clim_dat_current)
    clim_dat_rcp85 <- as.data.frame(clim_dat_rcp85)
    
    predGI5 <- predict(gf, as.data.frame(clim_dat_lgm))
    projGI5 <- predict(gf, as.data.frame(clim_dat_current))
    
    # calculate euclidean distance between current and future genetic spaces
    genOffsetGI5 <- sqrt((projGI5[,1]-predGI5[,1])^2+
                         (projGI5[,2]-predGI5[,2])^2+
                         (projGI5[,3]-predGI5[,3])^2+
                         (projGI5[,4]-predGI5[,4])^2+
                         (projGI5[,5]-predGI5[,5])^2+
                         (projGI5[,6]-predGI5[,6])^2+
                         (projGI5[,7]-predGI5[,7])^2)
    
   # genOffsetGI5 <- sqrt((projGI5$tmax-predGI5$tmax)^2)
    
    # Get offset values by climate variable
    offset <- NA
    for(col in 1:ncol(projGI5)){
      offset[col] <- mean(sqrt((projGI5[,col]-predGI5[,col])^2)) 
      names(offset)[col] <- colnames(projGI5)[col]
    }
    
    sort(offset, decreasing = T)

    summary(genOffsetGI5)
    hist(genOffsetGI5, breaks = 50)

    ## Correlate with predicted change in growth rates
    cors <- dplyr::left_join(dplyr::select(gen_dat_gf, id, locality, 
                                           wna_climate_vars_lgm, mems), 
                           dplyr::select(gen_dat_moms, prs_best_scaled, id))# %>%
      # group_by(locality) %>%
      # summarise_all(mean)

    cors$offset <-  genOffsetGI5
    cors$locality <- factor(cors$locality)
    table(cors$locality)
    dim(cors)
    
    
    # cors <- cors %>%
    #   group_by(locality) %>%
    #   summarise_all(mean)
    
    plot(cors$offset, cors$prs_best_scaled)
    cor.test(cors$offset, cors$prs_best_scaled)

    pairs.panels(cors[, c("offset", "PCNM1", "PCNM2", "PCNM3", "PCNM4", "PCNM5")])
    
    g1 <- gam(prs_best_scaled ~ s(offset, bs = "cr"), # +
             #     s(locality, bs = "re"),
              data = cors)
    summary(g1)
    visreg(g1, xvar = "offset", scale = "response", partial = T)
    visreg(g1, xvar = "offset", scale = "response", partial = F)
    
    # LMER version
    library(lmerTest)
    lm1 <- lmer(prs_best_scaled ~ offset + PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + (1 | locality),
                data = cors)
    lm1 <- lm(prs_best_scaled ~ offset + PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6,
                data = cors)
    summary(lm1)
    visreg(lm1, xvar = "offset")
    visreg(lm1, xvar = "offset", partial = F)
    
    
    
    
  ## Testing out plots - http://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf
    
    plot(cumimp(gf, "tmax_sum"), type = "l")
    plot(cumimp(gf, "cwd"), type = "l")
    plot(cumimp(gf, "bioclim_19"), type = "l")
    
    plot(gf, plot.type = "C", show.overall = F) # SNPs
    
    plot(gf, plot.type = "C", common.scale = T, show.species = F) # Overall
    
    plot(gf, plot.type = "P")
    
    
    
    
    
    
    
    

    ## Need to make rasters with each environmental variable
    ## And then extract this data to a dataframe to use for predictions in GF

    library(raster)

    tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
    tmin_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")

    process_raster <- function(rast){
      rast_temp <- aggregate(rast, fact = 10) # Make smaller
      rast_latlon <- projectRaster(rast_temp, crs = CRS("+proj=longlat +datum=WGS84"))
      plot(rast_latlon)
      return(rast_latlon)
    }

    tmax_rast <- process_raster(tmax_rast)
    tmin_rast <- process_raster(tmin_rast)

    # Find lat long points
    latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))

    gf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                          longitude = latlon[, 1],
                          latitude = latlon[, 2],
                          tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                          tmin_winter = raster::extract(tmin_rast, 1:ncell(tmin_rast)))
    gf_rast_df

    gf_rast_nona <-  gf_rast_df %>%
      filter(!is.na(tmax_sum))


    # transform env using gf models, see ?predict.gradientForest
    gf_pred <- predict(gf, gf_rast_nona[, colnames(gf_rast_nona) %in% climate_vars_gf]) # remove cell column before transforming


    # map continuous variation - reference SNPs
    refRGBmap <- pcaToRaster(gf_pred, gf_rast_df, tmax_rast, gf_rast_nona$cell_id)
    plotRGB(refRGBmap)

    # Mapping spatial genetic variation --------------------------------------------
    ###### functions to support mapping #####
    # builds RGB raster from transformed environment
    # snpPreds = dataframe of transformed variables from gf or gdm model
    # rast = a raster mask to which RGB values are to be mapped
    # cellNums = cell IDs to which RGB values should be assigned
    pcaToRaster <- function(gf_pred, gf_pred_wna, rast, mapCells){
      require(raster)

      pca <- prcomp(gf_pred, center=TRUE, scale.=FALSE)

      ##assigns to colors, edit as needed to maximize color contrast, etc.
      a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
      r <- a1+a2; g <- -a2; b <- a3+a2-a1

      ##scales colors
      scalR <- (r-min(r))/(max(r)-min(r))*255
      scalG <- (g-min(g))/(max(g)-min(g))*255
      scalB <- (b-min(b))/(max(b)-min(b))*255

      ## Join back into data frame that includes NAs
      temp <- data.frame(cell_id = mapCells,
                 scalR = scalR,
                 scalG = scalG,
                 scalB = scalB)

      temp_wna <- left_join(gf_pred_wna, temp, by = "cell_id")


      ## Check to make sure lengths match up
      if(nrow(temp_wna) != ncell(rast)){
        stop("Error: length of raster and of PCA predictions do not match")
      }

      ##assigns color to raster
      rast1 <- rast2 <- rast3 <- rast
      values(rast1) <- temp_wna$scalR
      values(rast2) <- temp_wna$scalG
      values(rast3) <- temp_wna$scalB
      ##stacks color rasters
      outRast <- stack(rast1, rast2, rast3)
      return(outRast)
    }

    # Function to map difference between spatial genetic predictions
    # predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
    # predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
    # rast = a raster mask to which Procrustes residuals are to be mapped
    # mapCells = cell IDs to which Procrustes residuals values should be assigned
    RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
      require(vegan)
      PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
      PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
      diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
      residMap <- residuals(diffProcrust)
      rast[mapCells] <- residMap
      return(list(max(residMap), rast))
    }

    # OK, on to mapping. Script assumes:
    # (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
    # and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
    #
    # (2) a raster mask of the study region to which the RGB data will be written

    # transform env using gf models, see ?predict.gradientForest
    predRef <- predict(gfRef, env_trns[,-1]) # remove cell column before transforming
    predGI5 <- predict(gfGI5, env_trns[,-1])

    # map continuous variation - reference SNPs
    refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
    plotRGB(refRGBmap)
    writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

    # map continuous variation - GI5 SNPs
    GI5RGBmap <- pcaToRaster(predGI5, mask, env_trns$cell)
    plotRGB(refRGBmap)
    writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)

    # Difference between maps (GI5 and reference)
    diffGI5 <- RGBdiffMap(predRef, predGI5, rast=mask, mapCells=env_trns$cell)
    plot(diffGI5[[2]])
    writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
    ################################################################################


    # Calculate and map "genetic offset" under climate change ----------------------
    # Script assumes:
    # (1) a dataframe of transformed env. variables for CURRENT climate
    # (e.g., predGI5 from above).
    #
    # (2) a dataframe named env_trns_future containing extracted raster data of
    # env. variables for FUTURE a climate scenario, same structure as env_trns

    # first transform FUTURE env. variables
    projGI5 <- predict(gfGI5, env_trns_future[,-1])

    # calculate euclidean distance between current and future genetic spaces
    genOffsetGI5 <- sqrt((projGI5[,1]-predGI5[,1])^2+(projGI5[,2]-predGI5[,2])^2
                         +(projGI5[,3]-predGI5[,3])^2+(projGI5[,4]-predGI5[,4])^2
                         +(projGI5[,5]-predGI5[,5])^2+(projGI5[,6]-predGI5[,6])^2
                         +(projGI5[,7]-predGI5[,7])^2)

    # assign values to raster - can be tricky if current/future climate
    # rasters are not identical in terms of # cells, extent, etc.
    mask[env_trns_future$cell] <- genOffsetGI5
    plot(mask)



    ################################################################################
    # END SCRIPTS FOR GRADIENT FOREST MODELING
    ################################################################################



  
  
  
 
#### Calculating Fst per snp
 
 library(adegenet)
 library(hierfstat)
 
 fst <- gen_dat_clim_all %>%
   dplyr::select(id, snp_col_names)  %>%
   dplyr::left_join(., dplyr::select(sample_key, gbs_name, locality), 
                    by = c("id" = "gbs_name")) %>%
   dplyr::filter(!(id %in% exclude_based_on_missing))
 
 dim(fst)
 

 fst[fst == 0] <- "1/1" # homozygote reference
 fst[fst == 1] <- "1/2" # homozygote reference
 fst[fst == 2] <- "2/2" # homozygote reference
 head(fst)
 
 fst_genind <- df2genind(dplyr::select(fst, -id, -locality), 
                         sep = "/", ploidy = 2,
                         pop = fst$locality)
 
 fst_hier <- genind2hierfstat(fst_genind)
 
 stats <- basic.stats(fst_hier)
 
 
 
 test = data.frame(locus = rownames(stats$perloc),
            fst = stats$perlo$Fst)
 
 test$outlier_type = "overall"
 test$outlier_type[test$locus %in% top_snps_long$snp] <- "beneficial"
 test$outlier_type[test$locus %in% mid_snps_long$snp] <- "overall"
 test$outlier_type[test$locus %in% bottom_snps_long$snp] <- "detrimental"
 
 test$outlier_type <- factor(test$outlier_type, levels = c("beneficial", "detrimental", "overall"))
 
 table(test$outlier_type)
 
 library(ggridges)
 
 ggplot(test, aes(fst, fill = outlier_type)) + 
   geom_histogram(col = "white", bins = 40) + 
   facet_wrap(~outlier_type, scales = "free_y", ncol = 1) +
   theme_bw(15)
 
 
 ggplot(test, aes(fst, y = outlier_type, fill = outlier_type)) + 
  # geom_density_ridges(col = "white", bins = 40, alpha = 0.75, scale = 2) + 
   stat_density_ridges(quantile_lines = TRUE, quantiles = 2, col = "grey10", alpha = .9) + 
 #  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE) +
   xlim(c(-.1, .25)) +
   scale_fill_manual(values =  c("#7fbf7b", "#af8dc3",  "grey50")) + 
#   facet_wrap(~outlier_type, scales = "free_y", ncol = 1) +
   theme_bw(15)
 
 test %>% 
   group_by(outlier_type) %>%
   summarise(mean = mean(fst),
                 sd = sd(fst),
             se = sd(fst)/sqrt(n()))
 
 test$outlier_type <- as.factor(test$outlier_type)
 
 kruskal.test(fst ~ outlier_type, data = test)
 
 t.test(test$fst[test$outlier_type == "beneficial"],
             test$fst[test$outlier_type == "overall"] )
 
 t.test(test$fst[test$outlier_type == "detrimental"],
             test$fst[test$outlier_type == "overall"] )

 

 
# Modeling spatial distributions of outlier alleles ------------------------
   
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
       rast <- aggregate(rast, fact = 10) # Make smaller by averaging across cells
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
 
   ## Remove cells that are out of climate range
     
     # # Calculate ranges
     # tmax_sum_range <- range(unlist(raster::extract(tmax_rast, lobata_range_rough)), na.rm = TRUE)
     # tmin_winter_range <- range(unlist(raster::extract(tmin_winter, lobata_range_rough)), 
     #                            na.rm = TRUE)
     # aet_range <- range(unlist(raster::extract(aet, lobata_range_rough)), na.rm = TRUE)
     # cwd_range <- range(unlist(raster::extract(cwd, lobata_range_rough)), na.rm = TRUE)
     # elevation_range <- range(unlist(raster::extract(dem, lobata_range_rough)), na.rm = TRUE)
     # bioclim_04_range <- range(unlist(raster::extract(bioclim_04, lobata_range_rough)),
     #                           na.rm = TRUE)
     # bioclim_15_range <- range(unlist(raster::extract(bioclim_15, lobata_range_rough)),
     #                           na.rm = TRUE)
     # bioclim_18_range <- range(unlist(raster::extract(bioclim_18, lobata_range_rough)),
     #                           na.rm = TRUE)
     # bioclim_19_range <- range(unlist(raster::extract(bioclim_19, lobata_range_rough)),
     #                           na.rm = TRUE)
     
     
     # # Set percentage threshold for range extension of climate and spatial variables
     # range_thresh <- 0.05
     # 
     # gam_dist_rast_df_nona <- gam_dist_rast_df_nona %>%
     #   dplyr::filter(tmax_sum >= tmax_sum_range[1] - (abs(tmax_sum_range[1]) * range_thresh) &
     #                   tmax_sum <= tmax_sum_range[2] + (abs(tmax_sum_range[1]) * range_thresh)) %>%
     #   dplyr::filter(tmin_winter >= tmin_winter_range[1] - (abs(tmin_winter_range[1]) * range_thresh) &
     #                 tmin_winter <= tmin_winter_range[2] + (abs(tmin_winter_range[1]) * range_thresh)) %>%
     #   dplyr::filter(cwd >= cwd_range[1] - (abs(cwd_range[1]) * range_thresh) &
     #                   cwd <= cwd_range[2] + (abs(cwd_range[1]) * range_thresh)) %>%
     #   dplyr::filter(bioclim_04 >= bioclim_04_range[1] - (abs(bioclim_04_range[1]) * range_thresh) &
     #                   bioclim_04 <= bioclim_04_range[2] + (abs(bioclim_04_range[1]) * range_thresh)) %>%
     #   dplyr::filter(bioclim_15 >= bioclim_15_range[1] - (abs(bioclim_15_range[1]) * range_thresh) &
     #                   bioclim_15 <= bioclim_15_range[2] + (abs(bioclim_15_range[1]) * range_thresh)) %>%
     #   dplyr::filter(bioclim_18 >= bioclim_18_range[1] - (abs(bioclim_18_range[1]) * range_thresh) &
     #                   bioclim_18 <= bioclim_18_range[2] + (abs(bioclim_18_range[1]) * range_thresh)) %>%
     #   dplyr::filter(elevation >= elevation_range[1] - (abs(elevation_range[1]) * range_thresh) &
     #                   elevation <= elevation_range[2] + (abs(elevation_range[1]) * range_thresh)) %>%
     #   dplyr::filter(longitude >= extent(lobata_range_rough)@xmin - (abs(extent(lobata_range_rough)@xmin) * range_thresh) &
     #                   longitude <= extent(lobata_range_rough)@xmax + (abs(extent(lobata_range_rough)@xmax) * range_thresh)) %>%
     #   dplyr::filter(latitude >= extent(lobata_range_rough)@ymin - (abs(extent(lobata_range_rough)@ymin) * range_thresh) &
     #                   latitude <= extent(lobata_range_rough)@ymax + (abs(extent(lobata_range_rough)@ymax) * range_thresh))  
     # 
     # dim(gam_dist_rast_df_nona)


  # Choose climate variables   
     climate_vars_gam_dist <-  c("tmax_sum", 
                           "tmin_winter", 
                            "cwd",
                      #     "aet",
                        #   "bioclim_04",
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
     
     
   ## Calculate PCA on climate and spatial variables
     
    #  ## Scale climate variables based on GBS moms
    #  x = 1 # used in loop
    #  
    #  scaled_var_means_gam_dist <- NA
    #  scaled_var_sds_gam_dist <- NA
    #  
    #  gen_dat_top_scaled <- gen_dat_top
    # 
    #  for(var in climate_vars_gam_dist){
    #    
    #    # For all seedlings  
    #    scaled_var_means_gam_dist[x] <- mean(gen_dat_top[[var]], na.rm = TRUE)
    #    scaled_var_sds_gam_dist[x] <- sd(gen_dat_top[[var]], na.rm = TRUE)
    #    gen_dat_top_scaled[var] <- (gen_dat_top_scaled[var] - scaled_var_means_gam_dist[x]) / scaled_var_sds_gam_dist[x]
    # 
    #    x = x + 1
    #  }
    #  
    #  # Assign names  
    #  names(scaled_var_means_gam_dist) <- climate_vars_gam_dist
    #  names(scaled_var_sds_gam_dist) <-  climate_vars_gam_dist
    # 
    #  summary(gen_dat_top_scaled[ , climate_vars_gam_dist]) # make sure they are centered and scaled
    #  
    #  ## Look at interaction in latlong
    #  gen_dat_top_scaled$latlong <- gen_dat_top_scaled$latitude * gen_dat_top_scaled$longitude
    #  
    #  
    #  ## PCA on climate variables
    #  pca_clim = prcomp(gen_dat_top_scaled[ , climate_vars_gam_dist], center = FALSE, scale = FALSE) 
    #  
    #  biplot(pca_clim)
    #  summary(pca_clim)
    #  screeplot(pca_clim, bstick = TRUE)
    #  
    #  ## Visualize PCA results
    #    fviz_contrib(pca_clim, choice = "var", axes = 1)
    #    fviz_contrib(pca_clim, choice = "var", axes = 2)
    #    fviz_contrib(pca_clim, choice = "var", axes = 3)
    #    fviz_contrib(pca_clim, choice = "var", axes = 4)
    #    
    #   
    #    fviz_pca_var(pca_clim, axes = c(1, 2), col.var="contrib")
    #    fviz_pca_var(pca_clim, axes = c(1, 4), col.var="contrib")
    #    
    #    
    #    ## Save plot of contributions
    #    
    #    fviz_pca_var(pca_clim, axes = c(1, 2), col.var="contrib")
    #    ggsave(filename = paste0("./figs_tables/Figure S5 - climate PCA 1 ", 
    #                             Sys.Date(), ".png"),
    #           units = "cm",
    #           width = 20, height = 15)
    #    
    #    fviz_pca_var(pca_clim, axes = c(1, 3), col.var="contrib")
    #    
    #    ggsave(filename = paste0("./figs_tables/Figure S5 - climate PCA 2 ",
    #                             Sys.Date(), ".png"),
    #           units = "cm",
    #           width = 20, height = 15)
    #    
    #    
    #    
    #    
    #      ## PC1 - cwd, long, lat, bioclim_18
    #    ## PC2 - bioclim4, bioclim15, tmin winter
    #    ## PC3 - latlong, elevation, tmax_sum,
    #    ## PC4 - latlong, tmax_sum, elevation
    #    
    #    sort(pca_clim$rotation[, 1])
    #    sort(pca_clim$rotation[, 2])
    #    sort(pca_clim$rotation[, 3])
    #    sort(pca_clim$rotation[, 4])
    #    
    #    
    #  pcs <- as.data.frame(pca_clim$x[, 1:7])
    #    
    #  # Add PCs as columns to dataframes
    #  gen_dat_top <- bind_cols(gen_dat_top, pcs)
    #  gen_dat_bottom <- bind_cols(gen_dat_bottom, pcs)
    #  gen_dat_mid <- bind_cols(gen_dat_mid, pcs)
    #    
    #    
    # # Convert raster df to PC space
    #    
    #    ## Forward scale variables first
    #    for(var in climate_vars_gam_dist){
    #      
    #      # For all seedlings  
    #      gam_dist_rast_df_nona[var] <- forward_transform(gam_dist_rast_df_nona[var],
    #                                                 var = var, 
    #                                                 means = scaled_var_means_gam_dist,
    #                                                 sds = scaled_var_sds_gam_dist)
    #    }
    #    
    #   summary(gam_dist_rast_df_nona)
    #   
    #   ## Make lat long interaction
    #   gam_dist_rast_df_nona$latlong <- gam_dist_rast_df_nona$latitude * gam_dist_rast_df_nona$longitude
    #    
    #   # Generate predictions
    #   gam_dist_rast_df_nona = bind_cols(gam_dist_rast_df_nona,
    #                               data.frame(predict(pca_clim, gam_dist_rast_df_nona[ , -1])))
    #   
    #   head(gam_dist_rast_df_nona)
    #   
    #   summary(gam_dist_rast_df_nona)
    #   
    #  
    #   climate_vars_gam_dist <- c("PC1","PC2","PC3")
    #  
    #  
     
     
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
     
     
     rda_df2 <- left_join(rda_df, dplyr::select(gen_dat_moms, accession, prs_best_scaled) %>% 
                            distinct(accession, .keep_all = T) %>% 
                            mutate(accession = as.numeric(as.character(accession)))) %>%
     #  filter(!is.na(accession)) %>%
       filter(!is.na(prs_best_scaled))
     
  
     rda_df2$locality <- factor(rda_df2$locality)

         gam_dist <- gam(prs_best_scaled ~ s(tmax_sum, bs = "cr") +
                             s(tmin_winter, bs = "cr") +
                             s(cwd, bs = "cr") + 
                        #     s(bioclim_04, bs = "cr") + 
                             s(bioclim_15, bs = "cr") + 
                             s(bioclim_18, bs = "cr") +  
                             s(elevation, bs = "cr") + 
                             s(latitude, bs = "cr") + 
                             s(longitude, bs = "cr"), 
                       #      s(latitude, by = longitude, bs = "cr"),
                           #  s(locality, bs = "re"),
                    #  method = "fREML",
                     # discrete = TRUE,
                   #  select = TRUE,
                     data = rda_df2)
         
         summary(gam_dist)
         
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
             names(rast) <- snp_name # Name the raster
             
             
             rast2 <- mask(rast, lobata_range_extended)
             plot(rast2)
             
             hist(rast2, breaks = 50)
             
            
     
  ## Look at deviance explained
     dev_df %>%
       dplyr::select(-snp) %>%
       dplyr::group_by(mode) %>%
       dplyr::summarise_all(mean)
     
     summary(lm(deviance ~ mode, data = dev_df))
     
     anova(lm(deviance ~ mode, data = dev_df))
     
     dev_df %>%
       ggplot(., aes(x = mode, y = deviance, fill = mode)) +
       geom_boxplot() + geom_jitter(width = .1) + theme_bw(15) +
       NULL
     
     
     dev_df %>%
       ggplot(., aes(x = mode, y = rsq_adj, fill = mode)) +
       geom_boxplot() + geom_jitter(width = .1) + theme_bw(15) +
       NULL

   # # Partial plots
   #   pp_df$var <- as.character(pp_df$var)
   #   pp_df$var[pp_df$var == "aet"] <- "Actual Evap."
   #   pp_df$var[pp_df$var == "bioclim_04"] <- "Temp. seasonality"
   #   pp_df$var[pp_df$var == "bioclim_15"] <- "Precip. seasonality"
   #   pp_df$var[pp_df$var == "bioclim_18"] <- "Summer precip."
   #   pp_df$var[pp_df$var == "latitude"] <- "Latitude"
   #   pp_df$var[pp_df$var == "longitude"] <- "Longitude"
   #   pp_df$var[pp_df$var == "tmax_sum"] <- "Tmax"
   #   pp_df$var[pp_df$var == "tmin_winter"] <- "Tmin"
   #   pp_df$var[pp_df$var == "elevation"] <- "Elevation"
   #   pp_df$var[pp_df$var == "random"] <- "Random variable"
   #   pp_df$var <- factor(pp_df$var)
   #   table(pp_df$var)
     
     ## Reverse order of snps
    #  pp_df$mode
    #  pp_df$mode <- factor(pp_df$mode, levels = c("bottom_snps", 
    #                                              "random_snps", 
    #                                              "top_snps"))
    # table(pp_df$mode)
    
 # Looking at rasters
    stack_top
  
  ## Weighted means      
  # top_snps_stack <- weighted.mean(stack_top,
  #                                   w = dev_df$deviance[dev_df$mode == "top_snps"])
  # bottom_snps_stack <- weighted.mean(stack_bottom,
  #                                 w = dev_df$deviance[dev_df$mode == "bottom_snps"])
  
  # Sum
    # top_snps_stack <- sum(stack_top)
    # bottom_snps_stack <- sum(stack_bottom)
  
  # Regular means
    top_snps_stack <- mean(stack_top, na.rm = TRUE)
    bottom_snps_stack <- mean(stack_bottom, na.rm = TRUE)
    
  # Regular means - with only models with deviance explained > 1%
    # top_snps_stack <- mean(stack_top[[which(dev_df$deviance[dev_df$mode == "top_snps"] >= .01)]], na.rm = TRUE)
    # bottom_snps_stack <- mean(stack_bottom[[which(dev_df$deviance[dev_df$mode == "bottom_snps"] >= .01)]], na.rm = TRUE)  
    # 
    
  # Medians  
    top_snps_stack <- calc(stack_top, median, na.rm = TRUE)
    bottom_snps_stack <- calc(stack_bottom, median, na.rm = TRUE)
    
  # Crop to just lobata range  
    top_snps_stack <- mask(top_snps_stack,  lobata_range_extended)
    bottom_snps_stack <- mask(bottom_snps_stack,  lobata_range_extended)
    
    
  ## Figure 4 - spatial maps of outlier genotypes ####  
    
    ## Set raster themes
    # Height theme
    genoTheme <- modifyList(rasterTheme(region = brewer.pal('PRGn',
                                                              n = 11)),                   
                              list(regions=list(alpha=0.8)))

    # Set max pixels
    max_pixels = 250000
    
    ## Plots of changes in height - RCP 85 top
    levelplot(hill_cropped , margin = FALSE,
              maxpixels = max_pixels, 
              par.settings = GrTheme()) + 
      latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey25")) + 
      as.layer(levelplot(top_snps_stack, 
                         under = TRUE, 
                         maxpixels = max_pixels,
                         par.settings = genoTheme)) + 
      latticeExtra::layer(sp.polygons(lobata_range, col = "black", lwd = 0.5)) + 
      latticeExtra::layer(sp.polygons(lobata_range_extended, col = "black", lwd = 1.5))
    
    # Main plot
    dev.copy(png, paste0("./figs_tables/fig3/Figure 3 - beneficial alleles map ", Sys.Date(), ".png"),
             res = 300, units = "cm", width = 16, height = 16)
    dev.off()
    
    # Save a quick PDF version for the legend
    pdf(paste0("./figs_tables/fig3/Figure 3 - beneficial alleles map legend ", Sys.Date(), ".pdf"),
        width = 4, height = 3)
    levelplot(top_snps_stack, margin = FALSE,
              maxpixels = max_pixels,
              par.settings =  genoTheme)
    dev.off()

    
## Plot Detrimental SNPS
    
    levelplot(hill_cropped , margin = FALSE,
              maxpixels = max_pixels, 
              par.settings = GrTheme()) + 
      latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey25")) + 
      as.layer(levelplot(bottom_snps_stack, 
                         under = TRUE, 
                         maxpixels = max_pixels,
                         par.settings = genoTheme)) + 
      latticeExtra::layer(sp.polygons(lobata_range, col = "black", lwd = 0.5)) + 
      latticeExtra::layer(sp.polygons(lobata_range_extended, col = "black", lwd = 1.5))
    
    # Main plot
    dev.copy(png, paste0("./figs_tables/fig3/Figure 3 - detrimental alleles map ", Sys.Date(), ".png"),
             res = 300, units = "cm", width = 16, height = 16)
    dev.off()
    
    # Save a quick PDF version for the legend
    pdf(paste0("./figs_tables/fig3/Figure 3 - detrimental alleles map legend ", Sys.Date(), ".pdf"),
        width = 4, height = 3)
    levelplot(bottom_snps_stack, margin = FALSE,
              maxpixels = max_pixels,
              par.settings =  genoTheme)
    dev.off()
    
# Figure 3 - partial dependence plots ####

  pp_df %>%
    dplyr::filter(mode != "random_snps") %>%
    dplyr::filter(var != "Random variable") %>%
   ggplot(., 
          aes(x = x_val,
          #    y = prob,
           y = prob_centered,
           #   y = prob_scaled,
              group = snp, col = mode, fill = mode), alpha = 0.75) + 
  #  geom_line(lwd = .25, alpha = 0.55) + # Plot individual lines
    ylab("Allele frequency \n (centered)") + xlab("") + 
    theme_bw(8) + 
    # theme(legend.position="none") +
       scale_color_manual(values = c("#af8dc3", "#7fbf7b")) + 
     geom_smooth(aes(group = mode), lwd = .75) +
     facet_wrap(~var, scales = "free_x", ncol = 4) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          strip.background = element_rect(color = "white", fill="white"),
          strip.text.x = element_text(face="bold")) +
     NULL
  
  # # Save file
  ggsave(filename = paste0("./figs_tables/Figure 3 - partial plots ", Sys.Date(), ".pdf"),
         units = "cm",
         width = 14, height = 4, useDingbats = FALSE)


  

  
# Make spatial predictions across landscape for number of alleles vs growth -------------------------------
  
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
  
  ## Creates vector with list of file paths to all .tif raster files
  ## Directory path where future climate scenarios are located
    dir_name <- "./data/gis/climate_data/BCM/future/"
    raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    raster_files <- raster_files[grep("*[[:digit:]].tif$", raster_files)] # Only those with .tif extension
    raster_files <- raster_files[grep("jja_ave", raster_files)] # Only tmax sum files
    
    raster_files
  
  ## Select only 85 scenario
    raster_files <- raster_files[grep("rcp85", raster_files)]
    
    future_stack <- stack()
    x = 1
  
  ## Loop through raster files and add to stack
  for(file in raster_files){
    
    # Load in raster
    raster_temp <- raster(file)
    
    cat("Working on file:", file, "...\n")
    
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
  
  tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
  
  
  # Calculate difference between future and current
  future_stack_tmax_dif <- future_stack - tmax_rast
  names(future_stack_tmax_dif) <- names(future_stack) # Reassign name
  
  # Scale values
  future_stack_tmax_dif_scaled <- (future_stack_tmax_dif - scaled_var_means_gbs_only["tmax_sum_dif"]) /
    scaled_var_sds_gbs_only["tmax_sum_dif"]
  
  
  # Make data frame for prediction - one for top and one for bottom
  newdata_rast_top <- data.frame(section_block = "IFG_1",
                                 height_2014 = 0,
                                 accession = "1",
                                 PC1_gen = 0, PC2_gen = 0, PC3_gen = 0, 
                                 n_top_snps = forward_transform(x = high[1],
                                                                var = 1,
                                                                means = n_top_snps_mean,
                                                                sd = n_top_snps_sd),
                                 n_bottom_snps = forward_transform(x = 0,
                                                                   var = 1,
                                                                   means = n_bottom_snps_mean,
                                                                   sd = n_bottom_snps_sd))
  
  
  newdata_rast_mid <- data.frame(section_block = "IFG_1",
                                 height_2014 = 0,
                                 accession = "1",
                                 PC1_gen = 0, PC2_gen = 0, PC3_gen = 0, 
                                 n_top_snps = forward_transform(x = mid[1],
                                                                var = 1,
                                                                means = n_top_snps_mean,
                                                                sd = n_top_snps_sd),
                                 n_bottom_snps = forward_transform(x = mid[2],
                                                                   var = 1,
                                                                   means = n_bottom_snps_mean,
                                                                   sd = n_bottom_snps_sd))

  newdata_rast_bottom <- data.frame(section_block = "IFG_1",
                                    height_2014 = 0,
                                    accession = "1",
                                    PC1_gen = 0, PC2_gen = 0, PC3_gen = 0, 
                                    n_top_snps = forward_transform(x = 0,
                                                                   var = 1,
                                                                   means = n_top_snps_mean,
                                                                   sd = n_top_snps_sd),
                                    n_bottom_snps = forward_transform(x = low[2],
                                                                      var = 1,
                                                                      means = n_bottom_snps_mean,
                                                                      sd = n_bottom_snps_sd))
  ## Get null prediction of height with no change in tmax
  future_null <- tmax_rast
  names(future_null) <- "tmax_sum_dif"
  
  # Need to scale tmax when actual transfer distance is 0
  future_null[!is.na(future_null)] <- forward_transform(0, var = "tmax_sum_dif",
                                                        means = scaled_var_means_gbs_only,
                                                        sds = scaled_var_sds_gbs_only)
  future_null
  
  # Top snps
  future_null_height_top <- predict(future_null,
                                    model = gam_snp,
                                    const = newdata_rast_top,
                                    progress = "text",
                                    type = "response")
  
  future_null_height_top
  summary(future_null_height_top)
  
  # Mid snps
  future_null_height_mid <- predict(future_null,
                                       model = gam_snp,
                                       const = newdata_rast_mid,
                                       progress = "text",
                                       type = "response")
  
  future_null_height_mid
  summary(future_null_height_mid)

  
  # Bottom snps
  future_null_height_bottom <- predict(future_null,
                                       model = gam_snp,
                                       const = newdata_rast_bottom,
                                       progress = "text",
                                       type = "response")
  
  future_null_height_bottom
  summary(future_null_height_bottom)
  
  
  
  # Initialize stack for height changes
  future_stack_height_change_top <- future_stack_tmax_dif
  future_stack_height_change_bottom <- future_stack_tmax_dif
  future_stack_height_change_mid <- future_stack_tmax_dif
  future_stack_height_change_topbottom <- future_stack_tmax_dif
  
  # Loop over scenarios
  for(scenario in 1:nlayers(future_stack_tmax_dif)){
    
    cat("Working on scenario:", names(future_stack_tmax_dif)[scenario], "... \n")
    
    ## Use scaled raster
    rast_temp <- future_stack_tmax_dif_scaled[[scenario]]
    names(rast_temp) <- "tmax_sum_dif" # Rename so predict function works
    
    rast_temp_height_top <- predict(rast_temp,
                                    model = gam_snp,
                                    const = newdata_rast_top,
                                    progress = "text",
                                    type = "response")
    
    rast_temp_height_mid <- predict(rast_temp,
                                    model = gam_snp,
                                    const = newdata_rast_mid,
                                    progress = "text",
                                    type = "response")
    
    
    rast_temp_height_bottom <- predict(rast_temp,
                                       model = gam_snp,
                                       const = newdata_rast_bottom,
                                       progress = "text",
                                       type = "response")
    
    
    # Dif bw top and null 
    rast_temp_height_change_top <- (rast_temp_height_top - future_null_height_top) / future_null_height_top * 100
    future_stack_height_change_top[[scenario]] <- rast_temp_height_change_top
    names(future_stack_height_change_top)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
    # Dif bw mid and null 
    rast_temp_height_change_mid <- (rast_temp_height_mid - future_null_height_mid) / future_null_height_mid * 100
    future_stack_height_change_mid[[scenario]] <- rast_temp_height_change_mid
    names(future_stack_height_change_mid)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
    
    # Dif bw bottom and null 
    rast_temp_height_change_bottom <- (rast_temp_height_bottom - future_null_height_bottom) / future_null_height_bottom * 100
    future_stack_height_change_bottom[[scenario]] <- rast_temp_height_change_bottom
    names(future_stack_height_change_bottom)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
    
    # Dif bw top and bottom
    # rast_temp_height_change_topbottom <- (rast_temp_height_top - rast_temp_height_bottom) / rast_temp_height_bottom * 100
    # future_stack_height_change_topbottom[[scenario]] <- rast_temp_height_change_topbottom
    # names(future_stack_height_change_topbottom)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
  } # End scenario loop
  
  
  # Plot summaries across scenarios  
  # histogram(future_stack_height_change_topbottom)
  # bwplot(future_stack_height_difference)
  
  
  ## Average across emissions scenarios
  future_stack_height_change_top_avg <- mean(future_stack_height_change_top, na.rm = TRUE)
  future_stack_height_change_mid_avg <- mean(future_stack_height_change_mid, na.rm = TRUE)
  future_stack_height_change_bottom_avg <- mean(future_stack_height_change_bottom, na.rm = TRUE)
  #future_stack_height_change_topbottom_avg <- mean(future_stack_height_change_topbottom, na.rm = TRUE)
  
  ## Project to latlong
  future_stack_height_change_top_avg <- projectRaster(future_stack_height_change_top_avg, 
                                                      crs = CRS("+proj=longlat +datum=WGS84"))
  future_stack_height_change_mid_avg <- projectRaster(future_stack_height_change_mid_avg, 
                                                      crs = CRS("+proj=longlat +datum=WGS84"))
  future_stack_height_change_bottom_avg <- projectRaster(future_stack_height_change_bottom_avg, 
                                                         crs = CRS("+proj=longlat +datum=WGS84"))
  # future_stack_height_change_topbottom_avg <- projectRaster(future_stack_height_change_topbottom_avg, 
  #                                                           crs = CRS("+proj=longlat +datum=WGS84"))
  
  ## Crop them to extended valley oak range
  future_stack_height_change_top_avg <- mask(future_stack_height_change_top_avg, 
                                             lobata_range_extended)
  future_stack_height_change_mid_avg <- mask(future_stack_height_change_mid_avg, 
                                             lobata_range_extended)
  future_stack_height_change_bottom_avg <- mask(future_stack_height_change_bottom_avg, 
                                                lobata_range_extended)
  # future_stack_height_change_topbottom_avg <- mask(future_stack_height_change_topbottom_avg, 
  #                                                  cali_outline)
  
  
  quantile(future_stack_height_change_top_avg, probs = c(0.005, 0.995))
  quantile(future_stack_height_change_mid_avg, probs = c(0.005, 0.995))
  quantile(future_stack_height_change_bottom_avg, probs = c(0.005, 0.995))
 # quantile(future_stack_height_change_topbottom_avg, probs = c(0.01, 0.99))
  
  
## Figure 2 - spatial predictions ####
  
  ## Set raster themes
    ## hill shade theme: gray colors with semitransparency
    # hsTheme <- modifyList(GrTheme(),                   
    #                       list(axis.line = list(col = "transparent"))
    
    # Height theme
    heightTheme <- modifyList(rasterTheme(region = brewer.pal('RdYlBu',
                                                              n = 11)),                   
                              list(regions=list(alpha=0.8)))
    
    # Set max pixels
    max_pixels = 250000
    
  ## Plots of changes in height - RCP 85 top
    levelplot(hill_cropped , margin = FALSE,
              maxpixels = max_pixels, 
              par.settings = GrTheme()) + 
     latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey25")) + 
     as.layer(levelplot(future_stack_height_change_top_avg, 
                       under = TRUE, 
                       maxpixels = max_pixels,
                       par.settings = heightTheme)) + 
      latticeExtra::layer(sp.polygons(lobata_range, col = "black", lwd = 0.5)) + 
      latticeExtra::layer(sp.polygons(lobata_range_extended, col = "black", lwd = 1.5))
  
    # Main plot
    dev.copy(png, paste0("./figs_tables/Figure 2 - top genotypes height change map ", Sys.Date(), ".png"),
             res = 300, units = "cm", width = 16, height = 16)
    dev.off()
    
    # Save a quick PDF version for the legend
      pdf(paste0("./figs_tables/Figure 2 - top genotypes height change legend ", Sys.Date(), ".pdf"),
          width = 4, height = 3)
      levelplot(future_stack_height_change_top_avg, margin = FALSE,
                maxpixels = max_pixels,
                par.settings =  heightTheme)
      dev.off()
    
  
## Plots of changes in height - RCP 85 bottom
  levelplot(hill_cropped , margin = FALSE,
            maxpixels = max_pixels, 
            par.settings = GrTheme()) + 
    latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey25")) + 
    as.layer(levelplot(future_stack_height_change_bottom_avg, 
                       under = TRUE, 
                       maxpixels = max_pixels,
                       par.settings = heightTheme)) + 
    latticeExtra::layer(sp.polygons(lobata_range, col = "black", lwd = 0.5)) + 
    latticeExtra::layer(sp.polygons(lobata_range_extended, col = "black", lwd = 1.5))
  
      # Main plot
      dev.copy(png, paste0("./figs_tables/Figure 2 - bottom genotypes height change map ", Sys.Date(), ".png"),
               res = 300, units = "cm", width = 16, height = 16)
      dev.off()
      
      # Save a quick PDF version for the legend
      pdf(paste0("./figs_tables/Figure 2 - bottom genotypes height change legend ", Sys.Date(), ".pdf"),
          width = 4, height = 3)
      levelplot(future_stack_height_change_bottom_avg, margin = FALSE,
                maxpixels = max_pixels,
                par.settings =  heightTheme)
      dev.off()
  
  
  ## Compare histograms of heights
  height_change_vals <- rbind(data.frame(type = "top_snps",
                                           height_change = na.omit(values(future_stack_height_change_top_avg))),
                              data.frame(type = "bottom_snps", 
                                         height_change = na.omit(values(future_stack_height_change_bottom_avg))),
                              data.frame(type = "mid_snps", 
                                         height_change = na.omit(values(future_stack_height_change_mid_avg))))

  ggplot(height_change_vals, aes(height_change)) + 
    geom_histogram(aes(fill = type), col = "white", bins = 50, size = .01, position = "identity", alpha  = 0.8) + 
    # geom_density(aes(fill = type, col = type), alpha = 0.5) +
    xlab("% Change in relative growth rate \n in warmer temperatures") + ylab("Frequency") +
    # scale_x_continuous(breaks = c(seq(-10, 5, 2.5)),
    #                    limits = c(-10, 5)) +
    scale_fill_manual(values = c("#7fbf7b", "#af8dc3", "grey50")) + # Flipped from normal order
    scale_color_manual(values = c("#7fbf7b", "#af8dc3", "grey50")) + # Flipped from normal order
    theme_bw(8) + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    NULL
  
  ggsave(filename = paste0("./figs_tables/Figure 2 - histogram of height predictions ", Sys.Date(), ".pdf"),
         units = "cm", width = 6, height = 4)
     

# Gradient forest test ----------------------------------------------------
#   
#   library(gradientForest)
#   
#   climate_vars_gf <- c("tmax_sum", 
#                       # "tmax_sum_lgm",
#                       # "tmin_winter", 
#                      # "bioclim_04",
#                        #"tmin_winter_lgm", 
#                        #  "DD5",
#                      #  "DD5_lgm",
#                        "random", 
#                       "latitude",
#                      # "mem1",
#                     #  "mem2",
#                     #  "mem3",
#                     #  "mem4",
#                       "longitude", 
#                       "elevation")
#   
#   pairs.panels(gen_dat_clim[, climate_vars_gf])
#   
#   
#   ## RANDOM SUBSAMPLE OF SNPS
#   # gf_sample <- sample(snp_col_names, length(top_snps$snp))
#   
#   
#   gf_scaled <- flashpcaR::scale2(gen_dat_clim[, gf_sample], impute = 2)
#   
#   
#   gen_dat_clim_randomized <- gen_dat_clim[, climate_vars_gf]
#   
# 
#   ## Factors
#   gf_factor <- gen_dat_top[ ,top_snps$snp]
#   gf_factor <- gf_factor %>%
#                  replace(is.na(.), 0) %>%
#                 mutate_all(as.factor) 
#   
#   ## Count number of positives
#   pos_thresh <- 20
#   gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 1)) < pos_thresh)]
#   gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 0)) < pos_thresh)]
#   
#   dim(gf_factor)
#   summary(gf_factor)
#               
#   
#   x = 1
#   r2_ran <- NULL
#   num_pos <- NULL
#   
#   
#   maxLevel <- log2(0.368*nrow(gen_dat_clim_randomized)/2) #account for correlations, see ?gradientForest 
#   
#   
#   while(x < 10){
#     
#     cat("working on:", x, "...\n")
#     
#   #  gen_dat_clim_randomized <- gen_dat_clim_randomized[sample(1:nrow(gen_dat_clim_randomized)),]
#     
#       gf <- gradientForest(cbind(gen_dat_clim_randomized,
#                                        gf_factor ),
#                            predictor.vars = climate_vars_gf, 
#                           response.vars = colnames(gf_factor),
#                          #response.vars = gsub("snp_", "", sum_df$snp),
#                           trace = TRUE,
#                          ntree = 100, mtry = 3,
#                          maxLevel = maxLevel, # From Fitzpatrick R script
#                          corr.threshold = 0.5)
#   
#       gf$result # Rsquared of positive loci
#       gf$species.pos.rsq
#       
#       cat("\n", mean(gf$result), "\n")
#       
#       importance(gf)
#       
#       # IF NULL
#       if(length(gf) == 0){
#         if(x == 1) {next}
#         num_pos[x] <- 0
#         r2_ran[x] <- 0
#         importance_df <- rbind(importance_df, 0)
#         x = x+1
#         next
#       }
#       
#       num_pos[x] <- gf$species.pos.rsq
#       r2_ran[x] <- mean(gf$result)
#       if(x == 1) {importance_df <- importance(gf)}
#       importance_df <- rbind(importance_df, importance(gf))
#       x = x +1
#   }
#   
#   r2_ran
#   mean(r2_ran)
#   
#   num_pos
#   mean(num_pos)
#   
#   importance_df
#   sort(colMeans(importance_df), decreasing = TRUE)
#   
#   plot(gf, plot.type = "O")
#   
#   
#   ## Need to make rasters with each environmental variable
#   ## And then extract this data to a dataframe to use for predictions in GF
#   
#   library(raster)
#   
#   tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
#   tmin_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")
#   
#   process_raster <- function(rast){
#     rast_temp <- aggregate(rast, fact = 10) # Make smaller
#     rast_latlon <- projectRaster(rast_temp, crs = CRS("+proj=longlat +datum=WGS84"))
#     plot(rast_latlon)
#     return(rast_latlon)
#   }
#   
#   tmax_rast <- process_raster(tmax_rast)
#   tmin_rast <- process_raster(tmin_rast)
#   
#   # Find lat long points
#   latlon <- xyFromCell(tmax_rast, 1:ncell(tmax_rast))
#   
#   gf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
#                         longitude = latlon[, 1],
#                         latitude = latlon[, 2],
#                         tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
#                         tmin_winter = raster::extract(tmin_rast, 1:ncell(tmin_rast)))
#   gf_rast_df
#   
#   gf_rast_nona <-  gf_rast_df %>%
#     filter(!is.na(tmax_sum))
# 
#   
#   # transform env using gf models, see ?predict.gradientForest
#   gf_pred <- predict(gf, gf_rast_nona[, colnames(gf_rast_nona) %in% climate_vars_gf]) # remove cell column before transforming
#   
#   
#   # map continuous variation - reference SNPs
#   refRGBmap <- pcaToRaster(gf_pred, gf_rast_df, tmax_rast, gf_rast_nona$cell_id)
#   plotRGB(refRGBmap)
# 
#   # Mapping spatial genetic variation --------------------------------------------
#   ###### functions to support mapping #####
#   # builds RGB raster from transformed environment
#   # snpPreds = dataframe of transformed variables from gf or gdm model
#   # rast = a raster mask to which RGB values are to be mapped
#   # cellNums = cell IDs to which RGB values should be assigned
#   pcaToRaster <- function(gf_pred, gf_pred_wna, rast, mapCells){
#     require(raster)
#     
#     pca <- prcomp(gf_pred, center=TRUE, scale.=FALSE)
#     
#     ##assigns to colors, edit as needed to maximize color contrast, etc.
#     a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
#     r <- a1+a2; g <- -a2; b <- a3+a2-a1
#     
#     ##scales colors
#     scalR <- (r-min(r))/(max(r)-min(r))*255
#     scalG <- (g-min(g))/(max(g)-min(g))*255
#     scalB <- (b-min(b))/(max(b)-min(b))*255
#     
#     ## Join back into data frame that includes NAs
#     temp <- data.frame(cell_id = mapCells,
#                scalR = scalR,
#                scalG = scalG, 
#                scalB = scalB)
#     
#     temp_wna <- left_join(gf_pred_wna, temp, by = "cell_id")
#     
#     
#     ## Check to make sure lengths match up
#     if(nrow(temp_wna) != ncell(rast)){
#       stop("Error: length of raster and of PCA predictions do not match")
#     }
#     
#     ##assigns color to raster
#     rast1 <- rast2 <- rast3 <- rast
#     values(rast1) <- temp_wna$scalR
#     values(rast2) <- temp_wna$scalG
#     values(rast3) <- temp_wna$scalB
#     ##stacks color rasters
#     outRast <- stack(rast1, rast2, rast3)
#     return(outRast)
#   }
#   
#   # Function to map difference between spatial genetic predictions
#   # predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
#   # predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
#   # rast = a raster mask to which Procrustes residuals are to be mapped
#   # mapCells = cell IDs to which Procrustes residuals values should be assigned
#   RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
#     require(vegan)
#     PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
#     PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
#     diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
#     residMap <- residuals(diffProcrust)
#     rast[mapCells] <- residMap
#     return(list(max(residMap), rast))
#   }
#   
#   # OK, on to mapping. Script assumes:
#   # (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
#   # and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#   #
#   # (2) a raster mask of the study region to which the RGB data will be written
#   
#   # transform env using gf models, see ?predict.gradientForest
#   predRef <- predict(gfRef, env_trns[,-1]) # remove cell column before transforming
#   predGI5 <- predict(gfGI5, env_trns[,-1])
#   
#   # map continuous variation - reference SNPs
#   refRGBmap <- pcaToRaster(predRef, mask, env_trns$cell)
#   plotRGB(refRGBmap)
#   writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)
#   
#   # map continuous variation - GI5 SNPs
#   GI5RGBmap <- pcaToRaster(predGI5, mask, env_trns$cell)
#   plotRGB(refRGBmap)
#   writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)
#   
#   # Difference between maps (GI5 and reference) 
#   diffGI5 <- RGBdiffMap(predRef, predGI5, rast=mask, mapCells=env_trns$cell)
#   plot(diffGI5[[2]])
#   writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)
#   ################################################################################
#   
#   
#   # Calculate and map "genetic offset" under climate change ----------------------
#   # Script assumes:
#   # (1) a dataframe of transformed env. variables for CURRENT climate 
#   # (e.g., predGI5 from above).
#   #
#   # (2) a dataframe named env_trns_future containing extracted raster data of 
#   # env. variables for FUTURE a climate scenario, same structure as env_trns
#   
#   # first transform FUTURE env. variables
#   projGI5 <- predict(gfGI5, env_trns_future[,-1])
#   
#   # calculate euclidean distance between current and future genetic spaces  
#   genOffsetGI5 <- sqrt((projGI5[,1]-predGI5[,1])^2+(projGI5[,2]-predGI5[,2])^2
#                        +(projGI5[,3]-predGI5[,3])^2+(projGI5[,4]-predGI5[,4])^2
#                        +(projGI5[,5]-predGI5[,5])^2+(projGI5[,6]-predGI5[,6])^2
#                        +(projGI5[,7]-predGI5[,7])^2)
#   
#   # assign values to raster - can be tricky if current/future climate
#   # rasters are not identical in terms of # cells, extent, etc.
#   mask[env_trns_future$cell] <- genOffsetGI5
#   plot(mask)
#   
#   
#   
#   ################################################################################
#   # END SCRIPTS FOR GRADIENT FOREST MODELING
#   ################################################################################
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#     plot(gf, plot.type = "O")
#     
#     most_important <- names(importance(gf))
#     
#     pred <- expand.grid(latitude = seq(31, 45, by = .1))
#     pred$pred <- predict(gf, newdata = pred)
#     
#     pred
#     
#     plot(pred$latitude, pred$pred$latitude)
#   
#     
#     plot(gf, plot.type = "S", imp.vars = most_important,
#           leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
#           cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
#           0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
#     
#     plot(gf, plot.type = "C", imp.vars = most_important,
#      show.overall = F, legend = T, leg.posn = "topleft",
#      leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
#      cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,
#      0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, + 0.3, 0, 0)))
#     
#     plot(gf, plot.type = "C", imp.vars = most_important,
#          show.species = F, common.scale = T, cex.axis = 0.6,
#         cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
#        + 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,  0.3, 0, 0)))
#     
#     plot(gf, plot.type = "P", show.names = T, horizontal = F,
#          cex.axis = 1, cex.labels = 0.7, line = 2.5)
# 
#     
# 
# 
# # RDA on top snps ---------------------------------------------------------
#     
#     rda_scaled <- flashpcaR::scale2(gen_dat_clim[, top_snps$snp], impute = 2)
#     
#     # rda_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
#     
#     rda_full <- vegan::rda(rda_scaled ~ latitude + 
#                              longitude + elevation +
#                              tmax_sum + tmin_winter,
#                            data = gen_dat_clim)
#     
#     rda_full
#     summary(rda_full)
#     
#     anova(rda_full)
#     
#     anova(rda_full, by="axis", permutations = 100)
#     
#     anova(rda_full, 
#           by = "margin",
#           permutations = 100)
#     
#     round(summary(rda_full)$biplot[, c(1,2)], 2)
#     
#     plot(rda_full)
#     
#     
#     loadings <- summary(rda_full)$species # Extracts loadings for each SNP (named species here)
#     head(loadings)
#     
#     # Plot histogram of the loadings
#     hist(loadings[, 1], xlab = "RDA1", breaks = 30, col = "skyblue2", main = "RDA1")
#     
#     # Outliers function from Forester et al.
#     outliers <- function(loadings, sd_cutoff){
#       # find loadings +/-z sd from mean loading
#       lims <- mean(loadings) + c(-1, 1) * sd_cutoff * sd(loadings)
#       loadings[loadings < lims[1] | loadings > lims[2]] # locus names in these tails
#     }
#     
#     cand_rda1 <- outliers(loadings[, 1], sd_cutoff = 2) 
#     cand_rda1
#     
#     cand_rda2 <- outliers(loadings[, 2], sd_cutoff = 2)
#     cand_rda2
#     
#     outlier_snps_names <- c(names(cand_rda1), names(cand_rda2))
#     
#     # Extract row indices of outlier snps
#     sig_snps <- which(rownames(loadings) %in% outlier_snps_names) 
#     
#     sig_snps <- which(rownames(loadings) %in% top_snps$snp)
#     
#     # Color of outlier SNPs
#     snp_cols <- rep("grey90",  times = nrow(loadings)) # Set default color
#     snp_cols[sig_snps] <- "red"
#     
#     # Labels of outlier SNPs
#     snp_labels <- rep("", times = nrow(loadings))
#     snp_labels[sig_snps] <- outlier_snps_names
#     
#     # Make biplot
#     plot(rda_full, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1,1)) # Create empty plot 
#     points(rda_full, display = "species", pch = 21,
#            col = "black", bg = snp_cols, scaling = 3)
#     
#     points(rda_full, display = "species", pch = 21, # Overplot sig snps
#            col = "black", select = sig_snps, bg = "red", scaling = 3)
#     # text(rda_full, display = "species", labels = snp_labels, cex = .75, scaling = 3)
#     text(rda_full, scaling = 3, display = "bp", col="black") # Environmental vectors  
#     
#     
#     
#     
#     # Plot onto map
#     
#     
#     # Get map of California - already loaded in ggplot2 package
#     ca_map <- dplyr::filter(ggplot2::map_data("state"), region == "california")
#     
#     
#     # Get loadings for each individual
#     ind_loadings <- summary(rda_full)$sites
#     head(ind_loadings)
#     
#     
#     # Plot values for first RDA axis
#     ggplot(ca_map, aes(x = long, y = lat)) +
#       geom_polygon(color = "black", fill = "grey80") +
#       geom_point(data = gen_dat_clim, aes(x = longitude, y = latitude, 
#                                           fill = ind_loadings[, 1]), # Color by first RDA axis
#                  color = "black", pch = 21, size = 5) +
#       scale_fill_gradientn(colours = terrain.colors(10), 
#                            name = "RDA1") +
#       theme_bw() + coord_equal(ratio = 1)
#     
#     
#     
#     ## Calculate correlation between prob of finding beneficial allele and future climate change
#     
#     # Read in future tmax scenarios
#     tmax_CCSM4 <- raster("./data/gis/climate_data/BCM/future/CCSM4_rcp85/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915.tif")
#     tmax_CCSM4 <- projectRaster(tmax_CCSM4, crs = CRS("+proj=longlat +datum=WGS84"))
#     tmax_CCSM4 <- resample(tmax_CCSM4, tmax_rast)
#     compareRaster(tmax_CCSM4, tmax_rast)
#     levelplot(tmax_CCSM4, margin = FALSE)
#     
#     
#     tmax_CNRM <- raster("./data/gis/climate_data/BCM/future/CNRM_rcp85/tmx2070_2099jja_ave_CNRM_rcp85_1529358976/tmx2070_2099jja_ave_CNRM_rcp85_1529358976.tif")
#     tmax_CNRM <- projectRaster(tmax_CNRM, crs = CRS("+proj=longlat +datum=WGS84"))
#     tmax_CNRM <- resample(tmax_CNRM, tmax_rast)
#     compareRaster(tmax_CNRM, tmax_rast)
#     levelplot(tmax_CNRM, margin = FALSE)
#     
#     tmax_IPSL <- raster("./data/gis/climate_data/BCM/future/IPSL_rcp85/tmx2070_2099jja_ave_IPSL_rcp85_1529358959/tmx2070_2099jja_ave_IPSL_rcp85_1529358959.tif")
#     tmax_IPSL <- projectRaster(tmax_IPSL, crs = CRS("+proj=longlat +datum=WGS84"))
#     tmax_IPSL <- resample(tmax_IPSL, tmax_rast)
#     compareRaster(tmax_IPSL, tmax_rast)
#     levelplot(tmax_IPSL, margin = FALSE)
#     
#     tmax_MIROC <- raster("./data/gis/climate_data/BCM/future/MIROC_rcp85/tmx2070_2099jja_ave_MIROC_rcp85_1529357403/tmx2070_2099jja_ave_MIROC_rcp85_1529357403.tif")
#     tmax_MIROC <- projectRaster(tmax_MIROC, crs = CRS("+proj=longlat +datum=WGS84"))
#     tmax_MIROC <- resample(tmax_MIROC, tmax_rast)
#     compareRaster(tmax_MIROC, tmax_rast)
#     levelplot(tmax_MIROC, margin = FALSE)
#     
#     # Stack rasters
#     future_stack <- stack(tmax_CCSM4) 
#     tmax_CNRM,
#     tmax_IPSL, 
#     tmax_MIROC)
# 
# levelplot(future_stack, margin = FALSE)
# 
# # Calculate future differences in temperature
# future_stack_dif <- future_stack - tmax_rast 
# 
# future_stack_dif <- mean(future_stack_dif)
# levelplot(future_stack_dif, margin = FALSE)
# 
# 
# # Correlations in future temp differences and probability of finding allele
# data.frame(tmax_sum_dif = values(future_stack_dif),
#            prob_bene = values(top_snps_stack)) %>%
#   ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
#   geom_smooth(method = "lm") + theme_bw(15) +
#   geom_smooth(col = "forestgreen") 
# 
# cor(values(future_stack_dif), values(bene_stack),
#     use = "na.or.complete")
# 
# 
# # Just lobata range
# future_stack_dif_range <- mask(future_stack_dif, lobata_range_rough)
# 
# plot(future_stack_dif_range, bene_stack_range,
#      las = 1, ylab = "Probability of finding beneficial allele",
#      xlab = "Tmax_sum_dif")
# 
# cor(values(future_stack_dif_range), values(bene_stack_range),
#     use = "na.or.complete")
# 
# data.frame(tmax_sum_dif = values(future_stack_dif_range),
#            prob_bene = values(bene_stack_range)) %>%
#   ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
#   geom_smooth(method = "lm") + geom_smooth(col = "forestgreen") + theme_bw(15)
# 
# 
# 
# rasterCo  
# 
# library(spatialEco)
# 
# test <-  rasterCorrelation(future_stack_dif, bene_stack, s = 7)
# test
# levelplot(test, margin = FALSE)


    
    
    
    






# Simulate data to test out interaction effects in gams -------------------

  dat <- expand.grid(cat = c(1,2,3),
                     x = seq(-3, 3, by = .1),
                     dummy1 = seq(-3, 3, by = .1))
#  dat$cat <- factor(dat$cat)
  
  y_intercept <- 100
    
    
  dat$y <- -.5 + -.05*dat$x + -.15*dat$x^2 + rnorm(nrow(dat), sd = .5) + y_intercept
  
  # dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  
  # Make cat 1 upside down
  dat$y[dat$cat == 1] <- -.5 + -.05*dat$x[dat$cat == 1] + .05*dat$x[dat$cat == 1]^2 + rnorm(nrow(dat[dat$cat == 1,]), sd = .5) + y_intercept

  # Make cat 3 a straight line
 #  dat$y[dat$cat == 2] <- rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  # dat$y[dat$cat == 3] <- rnorm(nrow(dat[dat$cat == 3,]), sd = .5)
  # dat$y[dat$cat == 3] <- rnorm(nrow(dat[dat$cat == 3,]), sd = .05) + y_intercept - 2.5
  # dat$y[dat$cat == 5] <- rnorm(nrow(dat[dat$cat == 5,]), sd = .5) + y_intercept
  # 
  ggplot(dat, aes(x, y, col = factor(cat), fill = factor(cat))) + geom_point() + 
    theme_bw() + geom_smooth()
  
  ## Quadratic model
  # lm1 <- lm(y ~ x + I(x^2) + cat + x * cat + I(x^2)*cat, data = dat)
  # summary(lm1)
  # visreg(lm1, xvar = "x", by = "cat")
 
  # Run model
  gam1 <- bam(y ~ s(x, bs = "cr") + 
                s(cat, bs = "cr", k = 2) + 
                s(x, by = cat, bs = "cr") + 
                s(dummy1, bs = "cr") + s(x, by = dummy1, bs = "cr"), 
              method = "fREML", discrete = TRUE,
              data = dat)
  
  summary(gam1)
  
  anova(gam1)
 # anova.gam(gam1, test = "Cp")
  
  visreg(gam1, xvar = "x", by = "cat")
  visreg(gam1)
  
  gam1$coefficients
  
  ## Testing out tensor products
  gam1 <- bam(y ~ te(x,  k = 2) + te(cat, k = 2) + ti(x, cat, k = 2), 
              method = "fREML", discrete = TRUE,
              data = dat)
  
  
  
  ## Using s smoother will give a > 0.05 p value when the line is flat
  ## I think m = 1 sets comparison level so that it compares smoothers to each other (https://cran.r-project.org/web/packages/itsadug/vignettes/test.html)
  
  # Discrete is not reliable - p values changea lot among simulations
  # m = 1 doesn't seem to do much with 4 categories
  
  # P-values are usually reliable if the smoothing parameters are known, or the model is unpenalized. If smoothing parameters have been estimated then the p-values are typically somewhat too low under the null. This occurs because the uncertainty associated with the smoothing parameters is neglected in the calculations of the distributions under the null, which tends to lead to underdispersion in these distributions, and in turn to p-value estimates that are too low. (In simulations where the null is correct, I have seen p-values that are as low as half of what they should be.) Note however that tests can have low power if the estimated rank of the test statistic is much higher than the EDF, so that p-values can also be too low in some cases.
  
  # If it is important to have p-values that are as accurate as possible, then, at least in the single model case, it is probably advisable to perform tests using unpenalized smooths (i.e. s(...,fx=TRUE)) with the basis dimension, k, left at what would have been used with penalization. Such tests are not as powerful, of course, but the p-values are more accurate under the null. Whether or not extra accuracy is required will usually depend on whether or not hypothesis testing is a key objective of the analysis.
  
  
  
  
  
  

