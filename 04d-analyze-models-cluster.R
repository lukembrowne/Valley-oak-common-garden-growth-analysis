
# Run model setup ---------------------------------------------------------

  source("./04b-setup-models-cluster.R")

# Read in cluster run output ----------------------------------------------

  # Run with predictions
  path_to_summaries <- "./output/run_3700632_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_devdif_ld/model_summaries/"
  path_to_predictions <- "./output/run_3700632_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_devdif_ld/model_predictions/"
  
  ## Randomized phenotypes
    path_to_summaries <- "./output/run_3573817_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_rgrrand_1/model_summaries/" # 0 sig snps
    path_to_summaries <- "./output/run_3576839_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_rgrrand_2/model_summaries/" # 0 sig snps
    path_to_summaries <- "./output/run_3576840_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_rgrrand_3/model_summaries/" # 3 outlier snps
    path_to_summaries <-  "./output/run_3576842_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_rgrrand_4/model_summaries/" # 1 outlier snp
    path_to_summaries <- "./output/run_3576843_tmax_sum_dif_rgr_fREML_discrete_v3_tw_genpc_rgrrand_5/model_summaries/" # 5 outlier snps
  
  
  # Read in files
  sum_df_raw <- plyr::ldply(list.files(path_to_summaries, full = TRUE), read_csv)
  
  dim(sum_df_raw)
  sum_df_raw
  
  summary(sum_df_raw)
  
  
  ## Missing snps
  # Find if there indexes if needed to re-run manually or on the cluster
  which(!snp_col_names %in% sum_df_raw$snp)
  snp_col_names[which(!snp_col_names %in% sum_df_raw$snp)]
  
  # For backup
  sum_df <- sum_df_raw
  

  # Filter out models that did not converge
  table(sum_df_raw$converged)
  
  
  ## Read in model predictions
    pred_df_raw <- plyr::ldply(list.files(path_to_predictions, full = TRUE), 
                               read_csv, col_names = c("snp", "tmax_sum_dif_unscaled", 
                                                       "genotype", "pred", "se"))
    pred_df_raw$genotype <- as.character(pred_df_raw$genotype)
    
    head(pred_df_raw)
    dim(pred_df_raw)
    summary(pred_df_raw)
    
  ## Calculate height change in warmer temperatures for each snp and genotype
  
    base0 <- median(pred_df_raw$pred[pred_df_raw$tmax_sum_dif_unscaled == 0])
    
   pred_df_long =  pred_df_raw %>%
      dplyr::group_by(snp, genotype) %>%
      do(data.frame(height_change_warmer = mean((.$pred[.$tmax_sum_dif_unscaled > 0]  - .$pred[.$tmax_sum_dif_unscaled == 0] ) /  .$pred[.$tmax_sum_dif_unscaled == 0] ) * 100,
         height_change_warmer_base0 = mean((.$pred[.$tmax_sum_dif_unscaled > 0]  - base0 ) /  base0) * 100))
   
   head(pred_df_long)
   dim(pred_df_long)
   
   summary(pred_df_long$height_change_warmer)
    
   
 
 
 
# Sort and arrange data ---------------------------------------------------

  
  
  ## Convert to from wide to long format
  n_gen_long <- sum_df %>%
    dplyr::select(snp, n_gen0, n_gen1, n_gen2, acc_pred, n) %>%
    gather(key = "genotype", value = "n_gen", -snp, -acc_pred, -n) %>%
    mutate(genotype = gsub("n_gen", "", genotype))
  
  head(n_gen_long)
  
  p_val_long <- sum_df %>%
      dplyr::select(snp, p_val_gen_0, p_val_gen_1, p_val_gen_2) %>%
      gather(key = "genotype", value = "p_val", -snp) %>%
      mutate(genotype = gsub("p_val_gen_", "", genotype))
  
    head(p_val_long)
  
  height_change_long <- sum_df %>%
    dplyr::select(snp, height_change_gen_0, height_change_gen_1, height_change_gen_2) %>%
    gather(key = "genotype", value = "height_change", -snp) %>%
    mutate(genotype = gsub("height_change_gen_", "", genotype))
  
  head(height_change_long)
  
  height_4.8_long <- sum_df %>%
    dplyr::select(snp, height_4.8_gen_0, height_4.8_gen_1, height_4.8_gen_2) %>%
    gather(key = "genotype", value = "height_4.8", -snp) %>%
    mutate(genotype = gsub("height_4.8_gen_", "", genotype))
  
  head(height_4.8_long)
  
  height_0_long <- sum_df %>%
    dplyr::select(snp, height_0_gen_0, height_0_gen_1, height_0_gen_2) %>%
    gather(key = "genotype", value = "height_0", -snp) %>%
    mutate(genotype = gsub("height_0_gen_", "", genotype))
  
  head(height_0_long)
  
  
  
  
## Join all long dataframes together
  sum_df_long <- left_join(n_gen_long, p_val_long) %>%
    left_join(. , height_change_long) %>%
    left_join(., height_4.8_long) %>%
    left_join(., height_0_long) %>%
    left_join(., pred_df_long)
  
  head(sum_df_long)
  dim(sum_df_long)
  

  # Filter out based on number of gene copies
    gen_copies_thresh <- 200
    
    sum_df_long <- sum_df_long %>%
      dplyr::filter(n_gen >= gen_copies_thresh)
    
  ## Filter out heterozygotes?
    sum_df_long <- sum_df_long %>%
      dplyr::filter(genotype != "1")
  
  ## Remove those that predictions weren't made on accession 1
    sum_df_long <- sum_df_long %>%
      dplyr::filter(acc_pred == 1)
  
  
  
# Calculate height change based on average at 0 transfer
  
  sum_df_long <- sum_df_long %>%
    mutate(height_change_0base = (height_4.8 - mean(height_0)) / mean(height_0) * 100)
  
  # Compare height measurements
   # plot(sum_df_long$height_change_0base, sum_df_long$height_change, pch = ".")
   # plot(sum_df_long$height_change_0base, sum_df_long$height_change_warmer, pch = ".") 
   # abline(0,1); abline(h=0, lty = 2); abline(v=0, lty=2)
  

# Correct p and q values
  # Vignette - https://bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
  
  # Calculate q values and adjusted p values
  qvals <- qvalue(p = sum_df_long$p_val,
                  fdr.level = 0.05)
  
  pvals_adj <- p.adjust(sum_df_long$p_val, method = "fdr")
  
  summary(qvals)
  summary(pvals_adj)
  hist(pvals_adj, breaks = 40)

  qvals$pi0 # overall proportion of true null hypotheses
  table(qvals$significant)
  
  # hist(qvals)
  
  # Assign back into DF
  sum_df_long <- sum_df_long %>%
    mutate(q_val = qvals$qvalues) %>%
    mutate(p_val_adj = pvals_adj)
  
  head(sum_df_long)
  
  
  
## Identify top and bottom XXX% of SNPs
  
  percent_thresh <- 0.005
  
 top_snps_long <-  sum_df_long %>%
    dplyr::filter(p_val_adj <= 0.01) %>%
    dplyr::arrange(desc(height_change_warmer_base0)) %>%
 #  dplyr::arrange(desc(height_change)) %>%
    head(round(nrow(sum_df) * percent_thresh))
 
 top_snps_long$top_snp = TRUE # Make a column for T / F of top snp
 
 summary(top_snps_long)
 
 dim(top_snps_long)
 head(top_snps_long)
 
 bottom_snps_long <- sum_df_long %>%
   dplyr::filter(p_val_adj <= 0.01) %>%
   dplyr::arrange(height_change_warmer_base0) %>%
 #  dplyr::arrange(height_change) %>%
   head(round(nrow(sum_df) * percent_thresh))
 
 bottom_snps_long$bottom_snp = TRUE
 
 head(bottom_snps_long)
 
 # Also make a comparable dataframe of 'random snps'
 
   set.seed(129) # To make sure same SNPs are chosen everytime - reproducible
   
   mid_snps_long <- sum_df_long %>%
     dplyr::filter(p_val_adj <= 0.01) %>%
     dplyr::filter(!snp %in% c(top_snps_long$snp, bottom_snps_long$snp)) %>%
     # Make sure number of genotype copies is comparable
     dplyr::filter(n_gen <= quantile(c(top_snps_long$n_gen, bottom_snps_long$n_gen),
                                     0.9)) %>% 
     dplyr::filter(n_gen >= gen_copies_thresh) %>%
     sample_n(round(nrow(sum_df) * percent_thresh))
   
   
   ## Any Snps sampled 2x?
   sum(table(mid_snps_long$snp) > 1)
   
     head(mid_snps_long)
     summary(mid_snps_long$height_change_0base)

  
  summary(top_snps_long$n_gen)
  summary(mid_snps_long$n_gen)
  summary(bottom_snps_long$n_gen)
  
  ## Join back to main DF
   
   sum_df_long <- left_join(sum_df_long, dplyr::select(top_snps_long, snp, genotype, top_snp)) %>%
                  left_join(., dplyr::select(bottom_snps_long, snp, genotype, bottom_snp)) %>%
                  mutate(top_snp = replace_na(top_snp, FALSE),
                         bottom_snp = replace_na(bottom_snp, FALSE))
   
   head(sum_df_long)
   
   table(sum_df_long$top_snp)
   
   
## Figure 2 -- Plot overlay of outlier SNPS ---------------------------------------
   
   # Join outlier dataframes with prediction dataframes
   pred_df_outlier <- dplyr::left_join(top_snps_long, pred_df_raw) %>%
         dplyr::mutate(type = "top_snp") %>%
      bind_rows( dplyr::left_join(bottom_snps_long, pred_df_raw) %>%
         dplyr::mutate(type = "bottom_snp")) %>%
       bind_rows( dplyr::left_join(mid_snps_long, pred_df_raw) %>%
                    dplyr::mutate(type = "mid_snp"))%>%
         dplyr::mutate(snp_gen = paste0(snp,genotype))
     #    dplyr::filter(type == "bottom_snp") %>%
         
     
   ggplot(pred_df_outlier, aes(x = tmax_sum_dif_unscaled, y = pred,
                             group = factor(snp_gen),
                              col = factor(type))) +
     geom_vline(aes(xintercept = 0), lty = 2, size = .25) +
     geom_vline(xintercept = 4.8, size = .25) +
     geom_line(alpha = 0.1,
               show.legend = FALSE, size = .10) +
   #  geom_point(alpha = 0.35, size = .25) +
     geom_smooth(aes(x = tmax_sum_dif_unscaled,
                     y = pred, group = factor(type)),
                 se = FALSE, size = .5) +
     
     labs(col = "Genotype") +
     ylim(.225, .325) +
     ylab("Relative growth rate") +
     xlab("Tmax transfer distance") +
     theme_bw(8) + 
     scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5),
                        limits = c(0, 7.5)) +
    # facet_wrap(~type) +
   #  facet_wrap(~ type + snp) +
     scale_color_manual(values = c("#FF6633", "grey50", "#6699CC")) + 
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           axis.line = element_line(colour = "black"),
           legend.position = "none")   
   
   # Save to file
   ggsave(filename = paste0("./figs_tables/Figure 2 - outlier responses ", 
                            Sys.Date(), ".pdf"),
          units = "cm",
          height = 4, width = 6,
          useDingbats = FALSE )
 
  summary(top_snps_long$height_change_warmer_base0)
  summary(bottom_snps_long$height_change_warmer_base0)
  summary(mid_snps_long$height_change_warmer_base0)
  
  summary(top_snps_long$height_change_0base)
  summary(bottom_snps_long$height_change_0base)
  summary(mid_snps_long$height_change_0base)
  
  (mean(top_snps_long$height_4.8) - mean(bottom_snps_long$height_4.8)) / mean(bottom_snps_long$height_4.8) * 100
   
   

 
 
 
# Manhattan plots ---------------------------------------------------------
 
 # library(CMplot)
 
 man_df1 <- data.frame(SNP = snp_col_names, 
                       snp_pos)
 
 man_df1$index <- 1:nrow(man_df1) ## For x axis plotting
 
 man_df <- left_join(man_df1, sum_df_long, by = c("SNP" = "snp"))

 head(man_df)
 
 # Set alpha for p values
 man_df$p_val_adj_alpha <- ifelse(man_df$p_val_adj < 0.05, 1, .35)
 
 # Set up shapes
 man_df <- man_df %>%
   mutate(pch = ifelse(genotype == 0 | genotype == 2, 16, 17)) # Circles for homozygotes, trianges for hets
 
 man_df <- man_df %>%
   mutate(pch = ifelse((genotype == 0 & (top_snp | bottom_snp)) | 
                         (genotype == 2 & ( top_snp | bottom_snp)), 21, pch)) # Outlined shapes
 
 man_df <- man_df %>%
   mutate(pch = ifelse((genotype == 1 & (top_snp | bottom_snp)), 24, pch)) # Outlined shapes
 
 table(man_df$pch)
 
 ## Clean up chromosome names - clear scaffold names
 man_df <- man_df %>%
   mutate(chrom = ifelse(grepl("chr", chrom), chrom, ""))
 
 
 ## Figure 3 - Manhattan plot ####
 
 # Setup

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
                  y = height_change_warmer_base0,
                #  y = height_change,
                  col = height_change_warmer_base0,
                #  col = height_change,
                  label = SNP)) +

   # Vertical lines
   geom_vline(xintercept = chrom_breaks$index, col = "grey50", alpha = 0.8, size = .2) + 
   
   # Horizontal line
   geom_hline(yintercept = 0, lty = 2, size = .25) +
   
   # Add points
   geom_point(pch = man_df$pch, 
             # col = "black",
              alpha = man_df$p_val_adj_alpha, 
              size = .65) +
   
   # Overplot top and bottom SNPS
   geom_point(data = man_df[man_df$top_snp | man_df$bottom_snp, ],
              aes(x = index, 
                  y = height_change_warmer_base0,
                  bg = height_change_warmer_base0),
                  # y = height_change,
                  # bg = height_change),
                  col = "black",
                  pch = man_df[man_df$top_snp | man_df$bottom_snp, ]$pch,
                  size = .65 + .35) +
   
   
   
   # # Add snp names
   #   geom_text_repel(
   #     aes(index, y = height_change_gen_0),
   #     data = subset(man_df, height_change_gen_0 > 0 & p_val_adj_gen_0 < 0.05),
   #     col = "black"
   #   ) +
   #   geom_text_repel(
   #     aes(index, y = height_change_gen_1),
   #     data = subset(man_df, height_change_gen_1 > 0 & p_val_adj_gen_1 < 0.05),
   #     col = "black"
 #   ) +
 #   geom_text_repel(
 #     aes(index, y = height_change_gen_2),
 #     data = subset(man_df, height_change_gen_2 > 0 & p_val_adj_gen_2 < 0.05),
 #     col = "black"
 #   ) +
 
 # Scale colors
 # scale_colour_viridis_c(option = "D", direction = -1, name = "% Change in RGR \n w/ 4.8°C increase") + 
 # scale_fill_viridis_c(option = "D", direction = -1, 
 #                      name = "% Change in RGR \n w/ 4.8°C increase",
 #                      guide = FALSE) + 
 
 scale_colour_distiller(type = "seq", palette = "RdYlBu", 
                        direction = 1, 
                        aesthetics = c("color", "fill"),
                        name = "% Change in RGR \n with 4.8°C increase") + 
   
   # Theme options
   scale_x_continuous(expand = c(0.01, 0.01), 
                      breaks = x_axis_breaks$index,
                      labels = x_axis_breaks$chrom) + 
   scale_y_continuous(breaks = seq(-30, 20, 5)) + 
   ylab("") + 
   xlab("Chromosome") + 
   theme_bw(8) + 
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        # legend.position = "none") +
           )+
   # axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
   NULL
 
 # Save as file - PNG
 ggsave(paste0("./figs_tables/Figure 2 - manhattan plot interaction term wo labels_", 
                      Sys.Date(), ".png"),
          width = 16, height = 4, units = "cm", dpi = 300)
 
 # Save as PDF for legend
 ggsave(paste0("./figs_tables/Figure 2 - manhattan plot interaction term LEGEND", 
               Sys.Date(), ".pdf"),
        width = 10, height = 2)
 
 
 
 
 
 

# Convert outlier SNPs genotypes to precense absence --------------------------

 ## For top snps
 gen_dat_top <- gen_dat_clim_all
 
  for(snp in top_snps_long$snp){
    
    # Subset out good genotype
    ben_gen <- as.numeric(top_snps_long$genotype[top_snps_long$snp == snp])
    
    if(length(ben_gen) > 1) {
      stop("length of ben gen shouldn't be more than 1!!")
    }
    
    # Replace genotype
    gen_dat_top[, snp] <- apply(gen_dat_clim_all[, snp], 1, function(x) {
                if(is.na(x)){ NA
                } else if(x == ben_gen){ 1
                } else {  0
                }
              }) # End apply
 
  }
 
 ## Remove not top snps
  dim(gen_dat_top)

  
  
  ## For bottom snps
  gen_dat_bottom <- gen_dat_clim_all
  
  for(snp in bottom_snps_long$snp){
    
    # Subset out good genotype
    ben_gen <- as.numeric(bottom_snps_long$genotype[bottom_snps_long$snp == snp])
    
    if(length(ben_gen) > 1) {
      stop("length of ben gen shouldn't be more than 1!!")
    }
    
    # Replace genotype
    gen_dat_bottom[, snp] <- apply(gen_dat_clim_all[, snp], 1, function(x) {
      if(is.na(x)){ NA
      } else if(x == ben_gen){ 1
      } else {  0
      }
    }) # End apply
    
  }
  
  dim(gen_dat_bottom)
  
## For middle snps "random snps"
  gen_dat_mid <- gen_dat_clim_all
  
  for(snp in mid_snps_long$snp){
    
    # Subset out good genotype
    ben_gen <- as.numeric(mid_snps_long$genotype[mid_snps_long$snp == snp])
    
    if(length(ben_gen) > 1) {
      stop("length of ben gen shouldn't be more than 1!!")
    }
    
    # Replace genotype
    gen_dat_mid[, snp] <- apply(gen_dat_clim_all[, snp], 1, function(x) {
      if(is.na(x)){ NA
      } else if(x == ben_gen){ 1
      } else {  0
      }
    }) # End apply
    
  }
  
  dim(gen_dat_mid)
  
  
  
  
## Compare deviances ########################################################
  
  
  # Make a dataframe that has deviance differences  
  deviance_df <- sum_df_raw %>%
    dplyr::select(snp, dev_dif, dev_dif_random) %>%
    tidyr::gather("dev_type", "dev_dif", -snp) %>%
    dplyr::mutate(has_outlier_snp = 
                    ifelse(snp %in% c(top_snps_long$snp,
                                      bottom_snps_long$snp), "Outlier genotypes" , 
                           "Non-outlier genotypes"))
  
  deviance_df$has_outlier_snp[deviance_df$dev_type == "dev_dif_random"] <- "Random variable"
  
  
  ## Figure S4 - differences in deviances ####
  ggplot(deviance_df, aes(x = factor(has_outlier_snp), y = dev_dif, 
                          fill = factor(has_outlier_snp)) )+ 
    geom_boxplot(outlier.shape = NA, alpha = 1) + 
    #   geom_jitter(aes(col = factor(has_outlier_snp)), alpha = .5, size = .25, pch = 16) + 
    xlab("") + ylab("Difference in deviance") +
    scale_fill_discrete(name = "") +
    theme_bw(10) +
    NULL
  
  ggsave(paste0("./figs_tables/Figure S4 - deviance explained ", Sys.Date(), ".pdf"),
         units = "cm", height = 8, width = 15)
  
  
  ## Statistical tests
  deviance_df %>%
    dplyr::filter(has_outlier_snp %in% c("Outlier genotypes", "Non-outlier genotypes")) %>%
    wilcox.test(dev_dif ~ has_outlier_snp, data = .)
  
  deviance_df %>%
    dplyr::filter(has_outlier_snp %in% c("Outlier genotypes", "Random variable")) %>%
    wilcox.test(dev_dif ~ has_outlier_snp, data = .)
  
  
  
  
  

  
  
# Random forest on just beneficial alleles ------------------------
   
   library(randomForest)
   #library(AUC)
   library(raster)
   library(pdp)
  
  ## Count number of positives - where allele is found
  pos_thresh <- 1
   

  ## For TOP SNPS
     rf_top <- gen_dat_top[ , snp_col_names]
     dim(rf_top)
     
     ## Count positives
     pos <-  apply(rf_top[ ,top_snps_long$snp], 2, function(x) sum(x == 1, na.rm = TRUE) )
     summary(pos)
     hist(pos, breaks =20)
     
     
     # Remove SNPs below threshold
       # For presences
       ind1 <- which(apply(rf_top, 2, function(x) sum(x == 1, na.rm = TRUE))
                     < pos_thresh)
       if(length(ind1) > 0) rf_top <- rf_top[, -ind1] # Need to check if indexes were found
       
       # For absences
       ind2 <- which(apply(rf_top, 2, function(x) sum(x == 0, na.rm = TRUE)) < pos_thresh)
       if(length(ind2 ) > 0) rf_top <- rf_top[, -ind2]
       
       dim(rf_top)
       
       
       ## How many top SNPs pass this threshold?
       sum(top_snps_long$snp %in% colnames(rf_top))
       
       
   ## For BOTTOM SNPS
       rf_bottom <- gen_dat_bottom[ , snp_col_names]
       dim(rf_bottom)
       
       
       ## Count positives
       pos <-  apply(rf_top[ , bottom_snps_long$snp], 2, function(x) sum(x == 1, na.rm = TRUE) )
       summary(pos)
       hist(pos, breaks =20)
       
       
     # Remove SNPs below threshold
       # For presences
       ind1 <- which(apply(rf_bottom, 2, function(x) sum(x == 1, na.rm = TRUE))
                     < pos_thresh)
       if(length(ind1) > 0) rf_bottom <- rf_bottom[, -ind1] # Need to check if indexes were found
       
       # For absences
       ind2 <- which(apply(rf_bottom, 2, function(x) sum(x == 0, na.rm = TRUE)) < pos_thresh)
       if(length(ind2 ) > 0) rf_bottom <- rf_bottom[, -ind2]
       
       dim(rf_bottom)
       
       
       ## How many bottom SNPs pass this threshold?
       sum(bottom_snps_long$snp %in% colnames(rf_bottom))     
       
       
    ## For mid snps
       rf_mid<- gen_dat_mid[ , snp_col_names]
       dim(rf_mid)
       
       

   # Shape outlines
     cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
     
     lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
     
     lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                         proj4string = CRS("+proj=longlat +datum=WGS84"))
     
    # Elevation map
     dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
     
     # Climate rasters
     tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
     
     tmin_winter <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmn1951_1980djf_ave_HST_1513102910/tmn1951_1980djf_ave_HST_1513102910.tif")
     
     aet <- raster("./data/gis/climate_data/BCM/historical/1951-1980/aet1951_1980_ave_HST_1513103180/aet1951_1980_ave_HST_1513103180.tif")
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
     
     rf_rast_df <- data.frame(cell_id = 1:ncell(tmax_rast),
                              longitude = latlon[, 1],
                              latitude = latlon[, 2],
                              tmax_sum = raster::extract(tmax_rast, 1:ncell(tmax_rast)),
                              tmin_winter = raster::extract(tmin_winter, 1:ncell(tmin_winter)),
                              aet = raster::extract(aet, 1:ncell(aet)),
                              elevation = raster::extract(dem, 1:ncell(dem)),
                              bioclim_04 = raster::extract(bioclim_04, 1:ncell(bioclim_04)),
                              bioclim_15 = raster::extract(bioclim_15, 1:ncell(bioclim_15)),
                              bioclim_18 = raster::extract(bioclim_18, 1:ncell(bioclim_18)),
                              bioclim_19 = raster::extract(bioclim_19, 1:ncell(bioclim_19)),
                              random = rnorm(1:ncell(dem)))
    # rf_rast_df
     
     dim(rf_rast_df)
     
     # Remove NAs from dataframe or random forest will complain
     rf_rast_df_nona <-  rf_rast_df %>%
       dplyr::filter(!is.na(tmax_sum)) %>%
       dplyr::filter(!is.na(elevation)) %>%
       dplyr::filter(!is.na(aet)) %>%
       dplyr::filter(!is.na(bioclim_04))
     
     summary(rf_rast_df_nona)
     dim(rf_rast_df_nona)
   

  # Choose climate variables   
     climate_vars_rf <-  c("tmax_sum", 
                           "tmin_winter", 
                           # "cwd",
                           "aet",
                           "bioclim_04",
                           "bioclim_15",
                           "bioclim_18",
                       #    "bioclim_19", # R = 0.79 with aet; R = 0.84 with bioclim_18
                           "elevation",
                           "latitude", 
                           "longitude", # R = -0.84 with latitude
                            "random")
     
     pairs.panels(gen_dat_top[, climate_vars_rf])
     
     
# Initialize variables
     
     stack_top <- stack() # For Top SNPs
     stack_bottom <- stack() # For Top SNPs
     stack_random_snps <- stack() # For Random SNPs
     stack_random_env <- stack() # For randomizing environmental variables

     start_flag = FALSE
     pp_df <- NULL # Initialize
     
     reps = 1
     
     # Sample same number of SNPs as top snps in random forest
     n_snps_to_sample <- sum(top_snps_long$snp %in% colnames(rf_top)) 
 
     
## Begin Random forest loop         
  for(mode in c("top_snps", "bottom_snps", "random_snps")){   # , random_env
    
    # Start rep loop
    for(rep in 1:reps){  
    
    # Choose SNPS to loop over  - same SNPs for either top snps or random_env
      if(mode == "top_snps" | mode == "random_env") {
        snps <- top_snps_long$snp[top_snps_long$snp %in% colnames(rf_top)]
      }
      
      if(mode == "bottom_snps" | mode == "random_env") {
        snps <- bottom_snps_long$snp[bottom_snps_long$snp %in% colnames(rf_bottom)]
      }
      
      
      # If random_snps, choose same number of top snps, but snps that aren't top or bottom snps
      # But only do it on the first rep so that we remeasure the same 'random snps' each rep
      if(mode == "random_snps") {
        snps <- mid_snps_long$snp[mid_snps_long$snp %in% colnames(rf_mid)]
      }

      # Loop over SNPs
       for(snp in 1:length(snps)){
         
         snp_name <- snps[snp]
         
         cat("Working on mode:", mode, "rep: ", rep, " snp:", snp, "  ", snp_name, "...\n")
         
         
         # Get climate and snp data
         if(mode == "top_snps"){
           X = gen_dat_top[, climate_vars_rf] # Save climate variables
           y = pull(rf_top, snp_name) # Pull out snp values
         }
         
         if(mode == "bottom_snps"){
           X = gen_dat_bottom[, climate_vars_rf] # Save climate variables
           y = pull(rf_bottom, snp_name) # Pull out snp values
         }
         
         if(mode == "random_snps"){
           X = gen_dat_mid[, climate_vars_rf] # Save climate variables
           y = pull(rf_mid, snp_name) # Pull out snp values
         }

         ## IF mode == random_env permute CLIMATE VARIABLES 
         if(mode == "random_env") {
           for(col in 1:ncol(X)){
             X[, col] <- X[sample(1:nrow(gen_dat_top)), col]
           }
           
         }

        # Count number of precense / absence
         n1 <- sum(y == 1, na.rm = TRUE)
         n0 <- sum(y == 0, na.rm = TRUE)
         
         not_NAs <- !is.na(y) # True if not NA
         
         # Convert to factor
         y <- factor(y)
        
         # Error check 
         if(length(levels(y)) != 2){
           stop("Length of levels of y not equal to 2- something is wrong")
         }
         
        # Set up dataframe for random forest
         dat_rf <- cbind(y = y[not_NAs], X[not_NAs, ])
 
         # Run random forest
         rf <- randomForest::randomForest(y ~.,
                            data = dat_rf,
                            importance = TRUE,
                            #   nodesize = 10,
                               sampsize = c(min(n0, n1), min(n0, n1)),
                           # sampsize = c(10, 10),
                            ntree = 1000)

         # plot(test)
         # plotImpHistory(test)
         
         # print(importance(rf))
          # cat("OOB error rate: ", mean(rf$err.rate[, 1]), "\n")
          # cat("Class1 error rate: ", rf$confusion[2,3], "\n")
         
         print(rf)
         
         
         # ROC plot
         # # The closer the curve follows the left-hand border and then the top border of the ROC space, the more accurate the test. The closer the curve comes to the 45-degree diagonal of the ROC space, the less accurate the test.
         # roc0 <- roc(rf$votes[,1],factor(1 * (rf$y==0)))
         # roc1 <- roc(rf$votes[,2],factor(1 * (rf$y==1)))
         # 
         # plot(roc0,
         #      main= snp_name, col = "blue", lwd = 2)
         # plot(roc1, add = TRUE, col = "green", lwd = 2)
         # legend("topright", legend = c(0, 1), col = c("blue", "green"), lty = 1, lwd = 2)
         # 
         # 
         
         #  plot(rf)
         
         #  varImpPlot(rf)
         
         # Save importance and OOB error rates
         if(start_flag == FALSE){
           
          # # Importance 
           # importance_boruta <- data.frame(mode = mode, rep = rep, snp = snp_name,
           #                                 climate_var = rownames(attStats(boruta_rf)),
           #                                 meanImp = attStats(boruta_rf)$meanImp,
           #                                 decision = attStats(boruta_rf)$decision)
           
           
           importance_rf <- data.frame(mode = mode, rep = rep, snp = snp_name,
                                           climate_var =  rownames(randomForest::importance(rf, scale = FALSE)),
                                           randomForest::importance(rf, scale = FALSE))
           
           
            # Error rates
             error_df <- data.frame(mode = mode,
                                      rep = rep, 
                                      snp = snp_name, 
                                      tail(rf$err.rate, 1)) # Last tree error rates
           
         } else { 

           
          # importance_boruta <- rbind(importance_boruta, data.frame(mode = mode, rep = rep, snp = snp_name,
          #                                                          climate_var = rownames(attStats(boruta_rf)),
          #                                                          meanImp = attStats(boruta_rf)$meanImp,
          #                                                          decision = attStats(boruta_rf)$decision))
          
          importance_rf <- rbind(importance_rf, data.frame(mode = mode, rep = rep, 
                                                           snp = snp_name,
                                    climate_var =  rownames(randomForest::importance(rf, 
                                                                            scale = FALSE)),
                                               randomForest::importance(rf, scale = FALSE)))
          
            error_df <- rbind(error_df, data.frame(mode = mode, 
                                                   rep = rep, 
                                                   snp = snp_name, 
                                                   tail(rf$err.rate, 1)))
         } # End else for start flag
  
         # Plot partial plot
          for(climate_var in climate_vars_rf){
           
             pp_out <- pdp::partial(rf, pred.var = climate_var, 
                               plot = FALSE, which.class = 2L,
                               prob = TRUE, smooth = TRUE)
            
             pp_df_temp <- data.frame(snp = snp_name,
                                      mode = mode,
                                      rep = rep,
                                      var = climate_var,
                                      prob = pp_out[, 2],
                                      prob_centered = pp_out[, 2] - mean(pp_out[, 2]),
                                      prob_scaled = scale(pp_out[, 2]),
                                      x_val = pp_out[, 1],
                                      stringsAsFactors = FALSE)

  
           if(!exists("pp_df")){ # If pp_df doesn't exist yet
             pp_df = pp_df_temp
           } else {
             pp_df <- bind_rows(pp_df, pp_df_temp)
           }
           
          } # End climate vars loop

         
         #   
         # }
         
         
         ## Predict across a raster and add to stack
  
             # Second column is positives
        
             rf_rast_df_nona$pred = predict(rf, rf_rast_df_nona, type = "prob")[, 2]
             

             rf_rast_df_temp <- left_join(rf_rast_df, rf_rast_df_nona[, c("cell_id", "pred")])

             rast <- tmax_rast
             values(rast) <- rf_rast_df_temp$pred
             names(rast) <- snp_name # Name the raster

             #  levelplot(rast, margin = FALSE)
          if(mode == "top_snps") stack_top <- stack(stack_top, rast)
          if(mode == "bottom_snps") stack_bottom <- stack(stack_bottom, rast)
          if(mode == "random_snps")   stack_random_snps <- stack(stack_random_snps, rast)
          if(mode == "random_env")   stack_random_env <- stack(stack_random_env, rast)
        
           start_flag <- TRUE # Flip flag to true
  
       } # End loop over SNPs
    
     } # End rep
    
  } # End mode loop ####################
     
     
     
     
     
   ## Importance
     
     importance_rf_raw <- importance_rf
     pp_df_raw <- pp_df
     
     head(importance_rf)
     dim(importance_rf)
     
     ## Change variable names in importance and partial plots 
       table(importance_rf$climate_var)
       
       importance_rf$climate_var <- as.character(importance_rf$climate_var)
       
       importance_rf$climate_var[importance_rf$climate_var == "aet"] <- "Actual Evap."
       importance_rf$climate_var[importance_rf$climate_var == "bioclim_04"] <- "Temp. seasonality"
       importance_rf$climate_var[importance_rf$climate_var == "bioclim_15"] <- "Precip. seasonality"
       importance_rf$climate_var[importance_rf$climate_var == "bioclim_18"] <- "Summer precip."
       importance_rf$climate_var[importance_rf$climate_var == "latitude"] <- "Latitude"
       importance_rf$climate_var[importance_rf$climate_var == "longitude"] <- "Longitude"
       importance_rf$climate_var[importance_rf$climate_var == "tmax_sum"] <- "Tmax"
       importance_rf$climate_var[importance_rf$climate_var == "tmin_winter"] <- "Tmin"
       importance_rf$climate_var[importance_rf$climate_var == "elevation"] <- "Elevation"
       importance_rf$climate_var[importance_rf$climate_var == "random"] <- "Random"
       importance_rf$climate_var <- factor(importance_rf$climate_var)
       
       table(importance_rf$climate_var)
       
       ## Reverse order of snps
         importance_rf$mode
         importance_rf$mode <- factor(importance_rf$mode, levels = c("bottom_snps", 
                                                                     "random_snps",
                                                                     "top_snps"))
          table(importance_rf$mode)
       
       
       # Partial plots
         pp_df$var <- as.character(pp_df$var)
         pp_df$var[pp_df$var == "aet"] <- "Actual Evap."
         pp_df$var[pp_df$var == "bioclim_04"] <- "Temp. seasonality"
         pp_df$var[pp_df$var == "bioclim_15"] <- "Precip. seasonality"
         pp_df$var[pp_df$var == "bioclim_18"] <- "Summer precip."
         pp_df$var[pp_df$var == "latitude"] <- "Latitude"
         pp_df$var[pp_df$var == "longitude"] <- "Longitude"
         pp_df$var[pp_df$var == "tmax_sum"] <- "Tmax"
         pp_df$var[pp_df$var == "tmin_winter"] <- "Tmin"
         pp_df$var[pp_df$var == "elevation"] <- "Elevation"
         pp_df$var[pp_df$var == "random"] <- "Random"
         pp_df$var <- factor(pp_df$var)
         table(pp_df$var)
         
         ## Reverse order of snps
         pp_df$mode
         pp_df$mode <- factor(pp_df$mode, levels = c("bottom_snps", 
                                                     "random_snps", 
                                                     "top_snps"))
        table(pp_df$mode)
       
    ## rf importance 
     # Find order of most important variables for TOP_SNPS
     importance_rf_df_top <- importance_rf %>%
       group_by(mode, snp, climate_var ) %>%
       summarise_all(mean) %>% # Average across reps
       ungroup() %>%
       group_by(mode, climate_var) %>%
       dplyr::select(-snp, -rep) %>%
       group_by(mode, climate_var) %>%
       summarise_all(median) %>%
       arrange(mode, MeanDecreaseAccuracy) %>%
       filter(mode == "top_snps") %>%
       ungroup() %>%
       dplyr::select(-MeanDecreaseAccuracy, -MeanDecreaseGini, -X0, -X1, -mode, -snp) %>%
       dplyr::mutate(order = row_number()) 
     
     # Make a dataframe with mean importance values
     
     lwr <- function(x) {quantile(x, 0.05)} # To find lower and upper quantiles
     upr <- function(x) {quantile(x, 0.95)}
     se <- function(x) {sd(x)/sqrt(length(x))}
    
     importance_rf_means <- 
                 importance_rf %>%
                 group_by(mode, snp, climate_var) %>%  ## Average across reps by snp
                 summarise_all(mean) %>%
                  ungroup() %>%
                 dplyr::select(-snp) %>%
                 group_by(mode, climate_var) %>%
                summarise_all(funs(mean, median, lwr, upr, sd, se))  %>% # Calculate means and upr and lower quantiles
                 left_join(., importance_rf_df_top, by = c("climate_var"))
    
     
     
  ## Figure 4 - rf variable importance    ####
     # Make ggplot
     importance_rf %>%
     group_by(mode, snp, climate_var ) %>%
       summarise_all(funs(mean, median, sd)) %>% # Average across reps
       dplyr::filter(mode != "random_snps") %>%
       dplyr::mutate(climate_var = factor(climate_var,
                                          levels = importance_rf_df_top$climate_var)) %>%
     ggplot(., aes(climate_var, MeanDecreaseAccuracy_mean, fill = mode)) + 
     #  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
      
       geom_hline(yintercept = 0, lty = 2, size = .25) +
       geom_point(aes(col = mode),
                  position = position_jitterdodge(),
                  size = .55, alpha = .75, pch = 16) +
       geom_boxplot(outlier.shape = NA, notch = TRUE, alpha = 0.45, size = .1) + 
       # geom_point(data = importance_rf_means, 
       #            aes(climate_var, MeanDecreaseAccuracy_mean, fill = mode),
       #            position = position_dodge2(width = .65),
       #            size = 3,
       #            pch = 21) + 
      
       # geom_errorbar(aes(ymin = MeanDecreaseAccuracy_mean - MeanDecreaseAccuracy_sd,
       #                   ymax = MeanDecreaseAccuracy_mean + MeanDecreaseAccuracy_sd),
       #                    position = position_dodge(width = .5),
       #               width = 0.1, alpha = 0.5) +
       theme_bw(10) + 
     #  ggtitle("Variable importance") + 
       #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
       ylab("Mean decrease in accuracy") + xlab("") +
       scale_x_discrete(limits = rev(levels(mode))) +
       scale_fill_manual(values = c("#FF6633","#6699CC")) + 
       scale_color_manual(values = c("#FF6633","#6699CC")) + 
      coord_flip() + 
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             legend.position = "none") + 
       NULL
     
     # Save file
     ggsave(filename = paste0("./figs_tables/Figure 4 - variable importance ", 
                              Sys.Date(), ".pdf"),
          units = "cm",
          height = 8, width = 8,
            
            useDingbats = FALSE )

     
 # Looking at rasters
    stack_top
    
    ## Average across reps within each SNP
    # Using stack apply - give indices of each raster that's a repition of a SNP - should go in order from 1-length(top_snps) and then repeat XX number of reps
    stack_top_avg = stackApply(stack_top,
                               indices = rep(1:length(top_snps_long$snp), 
                                                        times = rep), 
                               mean, na.rm = TRUE)
    
    
    stack_bottom_avg = stackApply(stack_bottom, 
                               indices = rep(1:length(bottom_snps_long$snp), 
                                             times = rep), 
                               mean, na.rm = TRUE)
    
    
     
  ## Take mean of rasters, weighted by class1 error - higher weights to lower error
   #  top_snps_stack <- weighted.mean(stack_top,
   #                                  w = 1 - error_df$OOB[error_df$mode == "top_snps"])
   # bottom_snps_stack <- weighted.mean(stack_bottom,
   #                                  w = 1 - error_df$OOB[error_df$mode == "bottom_snps"])
   #  random_snps_stack <- weighted.mean(stack_random_snps,
   #                                     w = 1 - error_df$OOB[error_df$mode == "random_snps"])
   #  random_env_stack <- weighted.mean(stack_random_env,
   #                                    w = 1 - error_df$OOB[error_df$mode == "random_env"])

    top_snps_stack <- mean(stack_top_avg, na.rm = TRUE)
    bottom_snps_stack <- mean(stack_bottom_avg, na.rm = TRUE)
   # random_snps_stack <- mean(stack_random_snps, na.rm = TRUE)
    
    
  ## Figure 4 - spatial maps of outlier genotypes ####  
    
    # Set breaks for color palette
    color_breaks <-  seq(from = min(c(values(top_snps_stack), values(bottom_snps_stack)),
                                    na.rm = TRUE),
                         to = max(c(values(top_snps_stack), values(bottom_snps_stack)),
                             na.rm = TRUE), 
                         length.out = 9)
    
  ## Plot beneficial raster
     levelplot(top_snps_stack, margin = FALSE,
           #    main = "Probability of finding beneficial genotype",
               maxpixels = 1e6,
               par.settings =  rasterTheme(region = brewer.pal('PRGn', n = 9)),
               at = color_breaks) +
       latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
       latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.5)) + 
       latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 1.5))
     
     # Save plot
     dev.copy(png, paste0("./figs_tables/Figure 4 - beneficial alleles map ", Sys.Date(), ".png"),
              res = 300, units = "cm", width = 16, height = 16)
     dev.off()
     
     
     
   # Save a quick PDF version for the legend
       pdf(paste0("./figs_tables/Figure 4 - beneficial alleles map LEGEND ", Sys.Date(), ".pdf"),
          width = 4, height = 3)
         levelplot(top_snps_stack, margin = FALSE,
                 #  main = "Probability of finding beneficial genotype",
                   maxpixels = 1e3,
                   par.settings =  rasterTheme(region = brewer.pal('PRGn', n = 9)),
                   at = color_breaks)
       dev.off()
     
     
     
     ## Maybe make a PDF version that will only have the scale
     
     
     ## Plot raster
     levelplot(bottom_snps_stack, margin = FALSE,
            #   main = "Probability of finding detrimental genotype",
               maxpixels = 1e6,
               par.settings =  rasterTheme(region = brewer.pal('PRGn', n = 9)),
               at = color_breaks) +
       latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
       latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.5)) + 
       latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 1.5))
     
   # Save plot
   dev.copy(png, paste0("./figs_tables/Figure 4 - detrimental alleles map ", Sys.Date(), ".png"),
            res = 300, units = "cm", width = 16, height = 16)
     dev.off()
   
     
  # # Get values by region
  #  # Not 100% sure if extract is doing exactly what we want   
  #   prob_by_region <- data.frame(region = lobata_range$REGION,
  #                                prob   = unlist(lapply(raster::extract(top_snps_stack,
  #                                                         y = lobata_range), mean)))
  #   
  # 
  #   ggplot(prob_by_region, aes(region, prob)) + 
  #     geom_bar(stat = "identity", col = "black", fill = "steelblue2") + theme_bw(15)
  #    
     
  
  
 # Figure 4 - partial dependence plots ####

  pp_df %>%
    dplyr::filter(mode != "random_snps") %>%
    dplyr::mutate(var = factor(var, # Reorder based on importance
                               levels = rev(importance_rf_df_top$climate_var))) %>%
    group_by(mode, snp, var, x_val) %>% # Average across reps
    summarise_all(mean) %>%
   ggplot(., 
          aes(x = x_val,
              y = prob,
           #   y = prob_centered,
           #   y = prob_scaled,
              group = snp, col = mode, fill = mode), alpha = 0.75) + 
  #   geom_line(lwd = .25, alpha = 0.15) + 
    ylab("Probablity of occurrence") + xlab("") + 
    theme_bw(8) + 
    # theme(legend.position="none") +
       scale_color_manual(values = c("#FF6633", "#6699CC")) + 
     geom_smooth(aes(group = mode), lwd = .75) +
     facet_wrap(~var, scales = "free_x") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          strip.background = element_rect(color = "white", fill="white"),
          strip.text.x = element_text(face="bold")) +
     NULL
  
  # Save file
  ggsave(filename = paste0("./figs_tables/Figure 4 - partial plots ", Sys.Date(), ".pdf"),
         units = "cm",
         width = 10, height = 8, useDingbats = FALSE)
   

  
  
  
  
  
  
  
  
  
  # Figure 2 - growth vs copies of alleles ##################
  
  ## Calculate number of top and bottom SNPs 
  dat_snp_count <- dat_snp_unscaled %>%
    dplyr::select(accession, tmax_sum_dif, rgr, section_block, PC1_gen, PC2_gen, PC3_gen, height_2014) %>%
    dplyr::mutate(accession = as.numeric(as.character(accession))) %>% # Change accession to numeric
    dplyr::left_join(., dplyr::select(gen_dat_top, accession, top_snps_long$snp)) %>%
    mutate(n_top_snps = rowSums(dplyr::select(., top_snps_long$snp), na.rm = TRUE)) %>%
    mutate(snp_chr7_9348369 = NULL) %>% # Remove this snp because it's in both top and bottom snps to not messu p counts
    dplyr::left_join(., dplyr::select(gen_dat_bottom, accession, bottom_snps_long$snp)) %>%
    mutate(n_bottom_snps = rowSums(dplyr::select(., bottom_snps_long$snp), na.rm = TRUE)) %>%
    mutate(accession = factor(as.character(accession)))
  
  dim(dat_snp_count)
  
  summary(dat_snp_count$n_top_snps)
  summary(dat_snp_count$n_bottom_snps)
  
  plot(dat_snp_count$n_bottom_snps, dat_snp_count$n_top_snps)
  cor.test(dat_snp_count$n_bottom_snps, dat_snp_count$n_top_snps)
  
  
  ## Scale number of top and bottom snps 
  n_top_snps_mean <- mean(dat_snp_count$n_top_snps)
  n_top_snps_sd <- sd(dat_snp_count$n_top_snps)
  dat_snp_count$n_top_snps_unscaled <- dat_snp_count$n_top_snps
  dat_snp_count$n_top_snps <- (dat_snp_count$n_top_snps - n_top_snps_mean) / n_top_snps_sd
  
  n_bottom_snps_mean <- mean(dat_snp_count$n_bottom_snps)
  n_bottom_snps_sd <- sd(dat_snp_count$n_bottom_snps)
  dat_snp_count$n_bottom_snps_unscaled <- dat_snp_count$n_bottom_snps
  dat_snp_count$n_bottom_snps <- (dat_snp_count$n_bottom_snps - n_bottom_snps_mean) / n_bottom_snps_sd
  
  # Set formula
  fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + s(tmax_sum_dif, bs=\"cr\") + s(n_top_snps, bs=\"cr\") + s(tmax_sum_dif, by = n_top_snps, bs=\"cr\") + s(n_bottom_snps, bs=\"cr\") + s(tmax_sum_dif, by = n_bottom_snps, bs=\"cr\") + s(accession, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC3_gen, bs=\"cr\") + s(tmax_sum_dif, by = PC1_gen, bs=\"cr\") + s(tmax_sum_dif, by = PC2_gen, bs=\"cr\") + s(tmax_sum_dif, by = PC3_gen, bs=\"cr\")"))
  
  
  # Gam with interaction effect
  gam_snp = bam(formula = formula(fixed_effects),
                data = dat_snp_count,
                discrete = TRUE, 
                nthreads = 8,
                method = "fREML",
                #   method = "ML",
                family = "tw",
                control = list(trace = FALSE))
  
  summary(gam_snp)
  
  
  ## Plot predictions
  
  newdata <-  expand.grid(section_block = "IFG_1",
                          height_2014 = 0,
                          accession = "1",
                          tmax_sum_dif = c(seq(min(dat_snp_unscaled$tmax_sum_dif),
                                               max(dat_snp_unscaled$tmax_sum_dif),
                                               length.out = 50)),
                          PC1_gen = 0, PC2_gen = 0, PC3_gen = 0, 
                          n_top_snps_unscaled = c(0,20),
                          n_bottom_snps_unscaled = c(0, 20))
  
  
  newdata$n_top_snps <- forward_transform(x = newdata$n_top_snps_unscaled,
                                          var = 1,
                                          means = n_top_snps_mean,
                                          sd = n_top_snps_sd)
  
  newdata$n_bottom_snps <- forward_transform(x = newdata$n_bottom_snps_unscaled,
                                             var = 1,
                                             means = n_bottom_snps_mean,
                                             sd = n_bottom_snps_sd)
  
  newdata$pred <- predict(gam_snp,
                          newdata = newdata,
                          se.fit = TRUE, type = "response")$fit
  
  newdata$se <- predict(gam_snp,
                        newdata = newdata,
                        se.fit = TRUE, type = "response")$se
  
  
  newdata %>%
    # dplyr::select(tmax_sum_dif, n_top_snps, n_bottom_snps, pred) %>%
    dplyr::mutate(type_count = paste0(n_top_snps_unscaled, ":", n_bottom_snps_unscaled)) %>%
    dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
                                                         means = scaled_var_means_gbs_only, sds = scaled_var_sds_gbs_only)) %>%
    dplyr::filter(type_count %in% c("20:0", "0:20")) %>%
    ggplot(., aes(x = tmax_sum_dif_unscaled, y = pred)) +
    geom_vline(aes(xintercept = 0), lty = 2, size = .25) +
    geom_vline(xintercept = 4.8, size = .25) +
    
    geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96 * se, fill = factor(type_count)), alpha = .25) +
    geom_line(aes(col = factor(type_count))) + 
    theme_bw(10) +
    
    scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
    # ylab(expression(Relative~growth~rate~(cm~cm^-1~yr^-1))) +
    ylab("Relative growth rate") +
    xlab("Tmax transfer distance") +
    ylim(c(.20, .4)) +
    scale_color_manual(values = c("#FF6633", "#6699CC")) +
    scale_fill_manual(values = c("#FF6633", "#6699CC")) +
    theme_bw(8) + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    NULL
  
  ggsave(filename = paste0("./figs_tables/Figure 2 - outlier counts ", Sys.Date(), ".pdf"),
         units = "cm", width = 6, height = 4)
  
  
  
  
  
  # Make spatial predictions across landscape -------------------------------
  
  # Load libraries ----------------------------------------------------------
  
  library(raster)
  library(rasterVis)
  
  
  # Loading in rasters ------------------------------------------------------
  
  ## Read in elevation dem and create a hillshade for mapping
  dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
  
  slope = raster::terrain(dem, opt='slope')
  aspect = raster::terrain(dem, opt='aspect')
  hill = hillShade(slope, aspect, 40, 270)
  
  ## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined_grouped.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
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
  
  ## Loop through raster files and add climate values to dat_mom dataframe
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
                                 n_top_snps = forward_transform(x = 20,
                                                                var = 1,
                                                                means = n_top_snps_mean,
                                                                sd = n_top_snps_sd),
                                 n_bottom_snps = forward_transform(x = 0,
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
                                    n_bottom_snps = forward_transform(x = 20,
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
    
    
    rast_temp_height_bottom <- predict(rast_temp,
                                       model = gam_snp,
                                       const = newdata_rast_bottom,
                                       progress = "text",
                                       type = "response")
    
    
    # Dif bw top and null 
    rast_temp_height_change_top <- (rast_temp_height_top - future_null_height_top) / future_null_height_top * 100
    future_stack_height_change_top[[scenario]] <- rast_temp_height_change_top
    names(future_stack_height_change_top)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
    
    # Dif bw bottom and null 
    rast_temp_height_change_bottom <- (rast_temp_height_bottom - future_null_height_bottom) / future_null_height_bottom * 100
    future_stack_height_change_bottom[[scenario]] <- rast_temp_height_change_bottom
    names(future_stack_height_change_bottom)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
    
    # Dif bw top and bottom
    rast_temp_height_change_topbottom <- (rast_temp_height_top - rast_temp_height_bottom) / rast_temp_height_bottom * 100
    future_stack_height_change_topbottom[[scenario]] <- rast_temp_height_change_topbottom
    names(future_stack_height_change_topbottom)[scenario] <- names(future_stack_tmax_dif)[scenario] # Reassign name
    
  } # End scenario loop
  
  
  # Plot summaries across scenarios  
  # histogram(future_stack_height_change_topbottom)
  # bwplot(future_stack_height_difference)
  
  
  ## Average across emissions scenarios
  future_stack_height_change_top_avg <- mean(future_stack_height_change_top, na.rm = TRUE)
  future_stack_height_change_bottom_avg <- mean(future_stack_height_change_bottom, na.rm = TRUE)
  future_stack_height_change_topbottom_avg <- mean(future_stack_height_change_topbottom, na.rm = TRUE)
  
  ## Project to latlong
  future_stack_height_change_top_avg <- projectRaster(future_stack_height_change_top_avg, 
                                                      crs = CRS("+proj=longlat +datum=WGS84"))
  future_stack_height_change_bottom_avg <- projectRaster(future_stack_height_change_bottom_avg, 
                                                         crs = CRS("+proj=longlat +datum=WGS84"))
  future_stack_height_change_topbottom_avg <- projectRaster(future_stack_height_change_topbottom_avg, 
                                                            crs = CRS("+proj=longlat +datum=WGS84"))
  
  ## Crop them
  future_stack_height_change_top_avg <- mask(future_stack_height_change_top_avg, 
                                             cali_outline)
  future_stack_height_change_bottom_avg <- mask(future_stack_height_change_bottom_avg, 
                                                cali_outline)
  future_stack_height_change_topbottom_avg <- mask(future_stack_height_change_topbottom_avg, 
                                                   cali_outline)
  
  
  
  summary(future_stack_height_change_topbottom)
  quantile(future_stack_height_change_top_avg, probs = c(0.005, 0.995))
  quantile(future_stack_height_change_bottom_avg, probs = c(0.005, 0.995))
  quantile(future_stack_height_change_topbottom_avg, probs = c(0.01, 0.99))
  
  
  ## Figure 2 - spatial predictions ####
  
  ## Plots of changes in height - RCP 85
  levelplot(future_stack_height_change_top_avg,
            margin = FALSE,
            maxpixels = 1000000,
            par.settings = rasterTheme(region = brewer.pal('RdYlBu', n = 9)),
            at = seq(-3, 2, length.out = 18),
            main = "Top genotypes") + 
    latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
    latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.5)) + 
    latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 1.5))
  
  # Main plot
  dev.copy(png, paste0("./figs_tables/Figure 2 - top genotypes map ", Sys.Date(), ".png"),
           res = 300, units = "cm", width = 16, height = 16)
  dev.off()
  
  # For legend
  dev.copy(pdf, paste0("./figs_tables/Figure 2 - top genotypes legend ", Sys.Date(), ".pdf"),
           width = 5, height = 3)
  dev.off()
  
  
  
  ## Plots of changes in height - RCP 85
  levelplot(future_stack_height_change_bottom_avg,
            margin = FALSE,
            maxpixels = 1000000,
            par.settings = rasterTheme(region = brewer.pal('RdYlBu', n = 9)),
            at = seq(-4, -9, length.out = 18),
            main = "Bottom genotypes") + 
    latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
    latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.5)) + 
    latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 1.5))
  
  # Main plot
  dev.copy(png, paste0("./figs_tables/Figure 2 - bottom genotypes map ", Sys.Date(), ".png"),
           res = 300, units = "cm", width = 16, height = 16)
  dev.off()
  
  # For legend
  dev.copy(pdf, paste0("./figs_tables/Figure 2 - bottom genotypes legend ", Sys.Date(), ".pdf"),
           width = 5, height = 3)
  dev.off()
  
  
  ## Compare histograms of heights
  
  
  height_change_vals <-   rbind(data.frame(type = "top_snps",
                                           height_change = na.omit(values(future_stack_height_change_top_avg))),
                                data.frame(type = "bottom_snps", 
                                           height_change = na.omit(values(future_stack_height_change_bottom_avg))))
  
  ggplot(height_change_vals, aes(height_change)) + 
    geom_histogram(aes(fill = type), col = "white", bins = 50, size = .01) + 
    # geom_density(aes(fill = type, col = type), alpha = 0.5) +
    xlab("% change in relative growth rate") + ylab("Frequency") +
    scale_x_continuous(breaks = c(seq(-10, 5, 2.5)),
                       limits = c(-10, 5)) +
    scale_fill_manual(values = c("#6699CC", "#FF6633")) + # Flipped from normal order
    scale_color_manual(values = c("#6699CC", "#FF6633")) + # Flipped from normal order
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
  
  library(gradientForest)
  
  climate_vars_gf <- c("tmax_sum", 
                      # "tmax_sum_lgm",
                      # "tmin_winter", 
                     # "bioclim_04",
                       #"tmin_winter_lgm", 
                       #  "DD5",
                     #  "DD5_lgm",
                       "random", 
                      "latitude",
                     # "mem1",
                    #  "mem2",
                    #  "mem3",
                    #  "mem4",
                      "longitude", 
                      "elevation")
  
  pairs.panels(gen_dat_clim[, climate_vars_gf])
  
  
  ## RANDOM SUBSAMPLE OF SNPS
  # gf_sample <- sample(snp_col_names, length(top_snps$snp))
  
  
  gf_scaled <- flashpcaR::scale2(gen_dat_clim[, gf_sample], impute = 2)
  
  
  gen_dat_clim_randomized <- gen_dat_clim[, climate_vars_gf]
  

  ## Factors
  gf_factor <- gen_dat_top[ ,top_snps$snp]
  gf_factor <- gf_factor %>%
                 replace(is.na(.), 0) %>%
                mutate_all(as.factor) 
  
  ## Count number of positives
  pos_thresh <- 20
  gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 1)) < pos_thresh)]
  gf_factor <- gf_factor[, -which(apply(gf_factor, 2, function(x) sum(x == 0)) < pos_thresh)]
  
  dim(gf_factor)
  summary(gf_factor)
              
  
  x = 1
  r2_ran <- NULL
  num_pos <- NULL
  
  
  maxLevel <- log2(0.368*nrow(gen_dat_clim_randomized)/2) #account for correlations, see ?gradientForest 
  
  
  while(x < 10){
    
    cat("working on:", x, "...\n")
    
  #  gen_dat_clim_randomized <- gen_dat_clim_randomized[sample(1:nrow(gen_dat_clim_randomized)),]
    
      gf <- gradientForest(cbind(gen_dat_clim_randomized,
                                       gf_factor ),
                           predictor.vars = climate_vars_gf, 
                          response.vars = colnames(gf_factor),
                         #response.vars = gsub("snp_", "", sum_df$snp),
                          trace = TRUE,
                         ntree = 100, mtry = 3,
                         maxLevel = maxLevel, # From Fitzpatrick R script
                         corr.threshold = 0.5)
  
      gf$result # Rsquared of positive loci
      gf$species.pos.rsq
      
      cat("\n", mean(gf$result), "\n")
      
      importance(gf)
      
      # IF NULL
      if(length(gf) == 0){
        if(x == 1) {next}
        num_pos[x] <- 0
        r2_ran[x] <- 0
        importance_df <- rbind(importance_df, 0)
        x = x+1
        next
      }
      
      num_pos[x] <- gf$species.pos.rsq
      r2_ran[x] <- mean(gf$result)
      if(x == 1) {importance_df <- importance(gf)}
      importance_df <- rbind(importance_df, importance(gf))
      x = x +1
  }
  
  r2_ran
  mean(r2_ran)
  
  num_pos
  mean(num_pos)
  
  importance_df
  sort(colMeans(importance_df), decreasing = TRUE)
  
  plot(gf, plot.type = "O")
  
  
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
  
  
  
  
  
  
  
  
  
  
    plot(gf, plot.type = "O")
    
    most_important <- names(importance(gf))
    
    pred <- expand.grid(latitude = seq(31, 45, by = .1))
    pred$pred <- predict(gf, newdata = pred)
    
    pred
    
    plot(pred$latitude, pred$pred$latitude)
  
    
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

    


# RDA on top snps ---------------------------------------------------------
    
    rda_scaled <- flashpcaR::scale2(gen_dat_clim[, top_snps$snp], impute = 2)
    
    # rda_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names], impute = 2)
    
    rda_full <- vegan::rda(rda_scaled ~ latitude + 
                             longitude + elevation +
                             tmax_sum + tmin_winter,
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
    
    
    loadings <- summary(rda_full)$species # Extracts loadings for each SNP (named species here)
    head(loadings)
    
    # Plot histogram of the loadings
    hist(loadings[, 1], xlab = "RDA1", breaks = 30, col = "skyblue2", main = "RDA1")
    
    # Outliers function from Forester et al.
    outliers <- function(loadings, sd_cutoff){
      # find loadings +/-z sd from mean loading
      lims <- mean(loadings) + c(-1, 1) * sd_cutoff * sd(loadings)
      loadings[loadings < lims[1] | loadings > lims[2]] # locus names in these tails
    }
    
    cand_rda1 <- outliers(loadings[, 1], sd_cutoff = 2) 
    cand_rda1
    
    cand_rda2 <- outliers(loadings[, 2], sd_cutoff = 2)
    cand_rda2
    
    outlier_snps_names <- c(names(cand_rda1), names(cand_rda2))
    
    # Extract row indices of outlier snps
    sig_snps <- which(rownames(loadings) %in% outlier_snps_names) 
    
    sig_snps <- which(rownames(loadings) %in% top_snps$snp)
    
    # Color of outlier SNPs
    snp_cols <- rep("grey90",  times = nrow(loadings)) # Set default color
    snp_cols[sig_snps] <- "red"
    
    # Labels of outlier SNPs
    snp_labels <- rep("", times = nrow(loadings))
    snp_labels[sig_snps] <- outlier_snps_names
    
    # Make biplot
    plot(rda_full, type = "n", scaling = 3, xlim = c(-1, 1), ylim = c(-1,1)) # Create empty plot 
    points(rda_full, display = "species", pch = 21,
           col = "black", bg = snp_cols, scaling = 3)
    
    points(rda_full, display = "species", pch = 21, # Overplot sig snps
           col = "black", select = sig_snps, bg = "red", scaling = 3)
    # text(rda_full, display = "species", labels = snp_labels, cex = .75, scaling = 3)
    text(rda_full, scaling = 3, display = "bp", col="black") # Environmental vectors  
    
    
    
    
    # Plot onto map
    
    
    # Get map of California - already loaded in ggplot2 package
    ca_map <- dplyr::filter(ggplot2::map_data("state"), region == "california")
    
    
    # Get loadings for each individual
    ind_loadings <- summary(rda_full)$sites
    head(ind_loadings)
    
    
    # Plot values for first RDA axis
    ggplot(ca_map, aes(x = long, y = lat)) +
      geom_polygon(color = "black", fill = "grey80") +
      geom_point(data = gen_dat_clim, aes(x = longitude, y = latitude, 
                                          fill = ind_loadings[, 1]), # Color by first RDA axis
                 color = "black", pch = 21, size = 5) +
      scale_fill_gradientn(colours = terrain.colors(10), 
                           name = "RDA1") +
      theme_bw() + coord_equal(ratio = 1)
    
    
    
    ## Calculate correlation between prob of finding beneficial allele and future climate change
    
    # Read in future tmax scenarios
    tmax_CCSM4 <- raster("./data/gis/climate_data/BCM/future/CCSM4_rcp85/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915/tmx2070_2099jja_ave_CCSM4_rcp85_1529358915.tif")
    tmax_CCSM4 <- projectRaster(tmax_CCSM4, crs = CRS("+proj=longlat +datum=WGS84"))
    tmax_CCSM4 <- resample(tmax_CCSM4, tmax_rast)
    compareRaster(tmax_CCSM4, tmax_rast)
    levelplot(tmax_CCSM4, margin = FALSE)
    
    
    tmax_CNRM <- raster("./data/gis/climate_data/BCM/future/CNRM_rcp85/tmx2070_2099jja_ave_CNRM_rcp85_1529358976/tmx2070_2099jja_ave_CNRM_rcp85_1529358976.tif")
    tmax_CNRM <- projectRaster(tmax_CNRM, crs = CRS("+proj=longlat +datum=WGS84"))
    tmax_CNRM <- resample(tmax_CNRM, tmax_rast)
    compareRaster(tmax_CNRM, tmax_rast)
    levelplot(tmax_CNRM, margin = FALSE)
    
    tmax_IPSL <- raster("./data/gis/climate_data/BCM/future/IPSL_rcp85/tmx2070_2099jja_ave_IPSL_rcp85_1529358959/tmx2070_2099jja_ave_IPSL_rcp85_1529358959.tif")
    tmax_IPSL <- projectRaster(tmax_IPSL, crs = CRS("+proj=longlat +datum=WGS84"))
    tmax_IPSL <- resample(tmax_IPSL, tmax_rast)
    compareRaster(tmax_IPSL, tmax_rast)
    levelplot(tmax_IPSL, margin = FALSE)
    
    tmax_MIROC <- raster("./data/gis/climate_data/BCM/future/MIROC_rcp85/tmx2070_2099jja_ave_MIROC_rcp85_1529357403/tmx2070_2099jja_ave_MIROC_rcp85_1529357403.tif")
    tmax_MIROC <- projectRaster(tmax_MIROC, crs = CRS("+proj=longlat +datum=WGS84"))
    tmax_MIROC <- resample(tmax_MIROC, tmax_rast)
    compareRaster(tmax_MIROC, tmax_rast)
    levelplot(tmax_MIROC, margin = FALSE)
    
    # Stack rasters
    future_stack <- stack(tmax_CCSM4) 
    tmax_CNRM,
    tmax_IPSL, 
    tmax_MIROC)

levelplot(future_stack, margin = FALSE)

# Calculate future differences in temperature
future_stack_dif <- future_stack - tmax_rast 

future_stack_dif <- mean(future_stack_dif)
levelplot(future_stack_dif, margin = FALSE)


# Correlations in future temp differences and probability of finding allele
data.frame(tmax_sum_dif = values(future_stack_dif),
           prob_bene = values(top_snps_stack)) %>%
  ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
  geom_smooth(method = "lm") + theme_bw(15) +
  geom_smooth(col = "forestgreen") 

cor(values(future_stack_dif), values(bene_stack),
    use = "na.or.complete")


# Just lobata range
future_stack_dif_range <- mask(future_stack_dif, lobata_range_rough)

plot(future_stack_dif_range, bene_stack_range,
     las = 1, ylab = "Probability of finding beneficial allele",
     xlab = "Tmax_sum_dif")

cor(values(future_stack_dif_range), values(bene_stack_range),
    use = "na.or.complete")

data.frame(tmax_sum_dif = values(future_stack_dif_range),
           prob_bene = values(bene_stack_range)) %>%
  ggplot(., aes(tmax_sum_dif, prob_bene)) + geom_point(size = 0.15, alpha = 0.25) +
  geom_smooth(method = "lm") + geom_smooth(col = "forestgreen") + theme_bw(15)



rasterCo  

library(spatialEco)

test <-  rasterCorrelation(future_stack_dif, bene_stack, s = 7)
test
levelplot(test, margin = FALSE)


    
    
    
    






# Simulate data to test out interaction effects in gams -------------------

  dat <- expand.grid(cat = c(1,2,3),
                     x = seq(-3, 3, by = .1))
  dat$cat <- factor(dat$cat)
  
  y_intercept <- 100
    
    
  dat$y <- -.5 + -.25*dat$x + -.25*dat$x^2 + rnorm(nrow(dat), sd = .5) + y_intercept
  
  # dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  
  # Make cat 1 upside down
 # dat$y[dat$cat == 1] <- -.5 + -.25*dat$x[dat$cat == 1] + .25*dat$x[dat$cat == 1]^2 + rnorm(nrow(dat[dat$cat == 1,]), sd = .5)
  dat$y[dat$cat == 2] <- -.5 + -.25*dat$x[dat$cat == 2] + .25*dat$x[dat$cat == 2]^2 + rnorm(nrow(dat[dat$cat == 2,]), sd = .5) + y_intercept
  
  # Make cat 2 a straight line
 #  dat$y[dat$cat == 2] <- rnorm(nrow(dat[dat$cat == 2,]), sd = .5)
  # dat$y[dat$cat == 3] <- rnorm(nrow(dat[dat$cat == 3,]), sd = .5)
   dat$y[dat$cat == 3] <- rnorm(nrow(dat[dat$cat == 3,]), sd = .05) + y_intercept
  # dat$y[dat$cat == 5] <- rnorm(nrow(dat[dat$cat == 5,]), sd = .5) + y_intercept
  # 
  ggplot(dat, aes(x, y, col = cat, fill = cat)) + geom_point() + 
    theme_bw() + geom_smooth()
  
  
  ## Quadratic model
  lm1 <- lm(y ~ x + I(x^2) + cat + x * cat + I(x^2)*cat, data = dat)
  summary(lm1)
  visreg(lm1, xvar = "x", by = "cat")
  
  
  
  # Run model
  gam1 <- bam(y ~ s(x, bs = "cr") + cat + s(x, by = cat, bs = "cr"), 
              method = "fREML", discrete = TRUE,
              data = dat)
  
  gam2 <- bam(y ~ s(x, bs = "cr") +cat, 
              method = "fREML", discrete = TRUE,
              data = dat)
  
  summary(gam1)
  
  anova(gam2, gam1, test = "F")
  
  anova(gam1)
  anova.gam(gam1, test = "Cp")
  
  ?visreg(gam1, xvar = "x", by = "cat")
  
  gam1$coefficients
  
  ## Using s smoother will give a > 0.05 p value when the line is flat
  ## I think m = 1 sets comparison level so that it compares smoothers to each other (https://cran.r-project.org/web/packages/itsadug/vignettes/test.html)
  
  # Discrete is not reliable - p values changea lot among simulations
  # m = 1 doesn't seem to do much with 4 categories
  
  # P-values are usually reliable if the smoothing parameters are known, or the model is unpenalized. If smoothing parameters have been estimated then the p-values are typically somewhat too low under the null. This occurs because the uncertainty associated with the smoothing parameters is neglected in the calculations of the distributions under the null, which tends to lead to underdispersion in these distributions, and in turn to p-value estimates that are too low. (In simulations where the null is correct, I have seen p-values that are as low as half of what they should be.) Note however that tests can have low power if the estimated rank of the test statistic is much higher than the EDF, so that p-values can also be too low in some cases.
  
  # If it is important to have p-values that are as accurate as possible, then, at least in the single model case, it is probably advisable to perform tests using unpenalized smooths (i.e. s(...,fx=TRUE)) with the basis dimension, k, left at what would have been used with penalization. Such tests are not as powerful, of course, but the p-values are more accurate under the null. Whether or not extra accuracy is required will usually depend on whether or not hypothesis testing is a key objective of the analysis.
  
  
  
  
  
  

