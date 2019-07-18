

# Load libraries ----------------------------------------------------------
  
  source("./01_clean-process-explore-data.R")
  source("./03_adding-climate-data.R")
  
  dim(dat_all_scaled)
  dim(dat_gbs_all_scaled)
  dim(dat_gbs_training_scaled)
  dim(dat_gbs_testing_scaled)


  library(gamm4)
  library(visreg)
  library(rrBLUP)
  library(qvalue) # devtools::install_github("jdstorey/qvalue")
  library(MuMIn)
  library(vegan)
  library(beepr)
  library(patchwork) # devtools::install_github("thomasp85/patchwork")
  library(rasterVis)
  library(pdp)
  library(ggrepel)
  library(caret)
  library(flashpcaR) # devtools::install_github("gabraham/flashpca/flashpcaR")
  
  
  # function to transform from raw to scaled variables
    forward_transform <- function(x, var, means, sds){
      ( x - means[var]) / sds[var]
    }
    
    # Function to back transform variables
    back_transform <- function(x, var, means, sds){
      (x * sds[var]) + means[var]
    }


# Goal is to run an individual gam for each SNP, like a fancy GWAS, and then correct for multiple testing after
    

# Format and save to .Rdata for cluster -----------------------------------

  # Do PCA on kinship matrix  
    
  # Joining climate and genetic data
  # Do an inner join to avoid having GBS moms with climate data but poor genotyping data
 
  # Just moms in common garden  
   gen_dat_clim <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim)
  
  # All GBS samples
  gen_dat_clim_all <- inner_join(climate_gbs_mom, 
                             dplyr::select(gen_dat_all, -accession), 
                             by = c("id" = "gbs_name"))
  dim(gen_dat_clim_all)
  
  
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
  
  
  ## Calculate kinship matrix
  
  # Need to first convert to -1, 0, 1 matrix
    kin_mat <- as.matrix(gen_dat_clim[, snp_col_names])
    kin_mat[kin_mat == 0] <- -1
    kin_mat[kin_mat == 1] <- 0
    kin_mat[kin_mat == 2] <- 1
    
    kin_mat <- A.mat(kin_mat) # Calculate kinship matrix
  
  
  ## Calculate PCA on kinship matrix
    pca_gen = prcomp(kin_mat, center = TRUE, scale = TRUE) ## PCA on kinship matrix
  
   # biplot(pca_gen)
    summary(pca_gen)
  #  screeplot(pca_gen, bstick = TRUE)
    
  ### Figure S3 - Plotting results of PCA for supplementary material

      # par(mfrow = c(2,1))
      # 
      # # Scale values for plotting in biplot - from: https://stats.stackexchange.com/questions/66926/what-are-the-four-axes-on-pca-biplot
      # scores<-pca_gen$x
      # lam <- pca_gen$sdev[1:3]
      # n <- NROW(scores)
      # lam <- lam * sqrt(n)
      # lam <- lam^1
      # yy<-t(t(pca_gen$rotation[, 1:3]) * lam)
      # xx<-t(t(scores[, 1:3])/lam)
      # 
      # # Pc1 vs pc2
      #   plot(xx[, 1], xx[, 2],
      #        las = 1,
      #        xlab = paste0("PC1: ", summary(pca_gen)$importance[2, 1] * 100, "% variance explained"),
      #        ylab = paste0("PC2: ", summary(pca_gen)$importance[2, 2] * 100, "% variance explained"),
      #        pch = 19,
      #        type = "n")
      #   text(xx[, 1], xx[, 2],
      #        labels = gen_dat_clim$id, cex = .5)
      #   mtext(text = "(a)", side = 3, adj = -.1, padj = -1, cex = 1.25)
      # 
      # 
      #   # Pc1 vs pc3
      #   plot(xx[, 1], xx[, 3],
      #        las = 1,
      #        xlab = paste0("PC1: ", summary(pca_gen)$importance[2, 1] * 100, "% variance explained"),
      #        ylab = paste0("PC3: ", summary(pca_gen)$importance[2, 3] * 100, "% variance explained"),
      #        pch = 19,
      #        type = "n")
      #   text(xx[, 1], xx[, 3],
      #        labels = gen_dat_clim$id, cex = .5)
      #   mtext(text = "(b)", side = 3, adj = -.1, padj = -1, cex = 1.25)
      # 
      #   dev.copy(png, filename = paste0("./figs_tables/Figure S3 - PCA plot_",
      #                                   Sys.Date(), ".png"),
      #            res = 300, width = 1800, height =2400)
      #   dev.off()

    pcs <- as.data.frame(pca_gen$x[, 1:10])
    colnames(pcs) <- paste0(colnames(pcs), "_gen")
    
    bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

    
 
    
 # Set up dataframe with SNP data, height data, and climate data
    dat_snp_all = left_join(dat_gbs_all_scaled,
                      bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                         pcs),  by = "accession")
    
  # Set up factors  
    dat_snp_all$accession <- factor(dat_snp_all$accession)
    dat_snp_all$section_block <- factor(dat_snp_all$section_block)
    dat_snp_all$section <- factor(dat_snp_all$section)
    dat_snp_all$locality <- factor(dat_snp_all$locality)
    
  # Scale PCA axes
    dat_snp_all$PC1_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC1_gen)) / sd(dat_snp_all$PC1_gen)
    dat_snp_all$PC2_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC2_gen)) / sd(dat_snp_all$PC2_gen)
    dat_snp_all$PC3_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC3_gen)) / sd(dat_snp_all$PC3_gen)
    dat_snp_all$PC4_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC4_gen)) / sd(dat_snp_all$PC4_gen)
    dat_snp_all$PC5_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC5_gen)) / sd(dat_snp_all$PC5_gen)
    
  # Save unscaled version  
    dat_snp_all_unscaled <- dat_snp_all
 
    summary(dat_snp_all[, 1050:1060])
    
  ## Testing set  
    dat_snp_testing = left_join(dat_gbs_testing_scaled,
                        bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                  pcs),  by = "accession")
    
    # Set up factors  
    dat_snp_testing$accession <- factor(dat_snp_testing$accession)
    dat_snp_testing$section_block <- factor(dat_snp_testing$section_block)
    dat_snp_testing$section <- factor(dat_snp_testing$section)
    dat_snp_testing$locality <- factor(dat_snp_testing$locality)
    
    # Scale PCA axes
    dat_snp_testing$PC1_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC1_gen)) / sd(dat_snp_testing$PC1_gen)
    dat_snp_testing$PC2_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC2_gen)) / sd(dat_snp_testing$PC2_gen)
    dat_snp_testing$PC3_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC3_gen)) / sd(dat_snp_testing$PC3_gen)
    dat_snp_testing$PC4_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC4_gen)) / sd(dat_snp_testing$PC4_gen)
    dat_snp_testing$PC5_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC5_gen)) / sd(dat_snp_testing$PC5_gen)
    
    # Save unscaled version  
    dat_snp_testing_unscaled <- dat_snp_testing
    
    summary(dat_snp_testing[, 1050:1060])    
    
 
  ## Training set  
    dat_snp_training = left_join(dat_gbs_training_scaled,
                                bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                          pcs),  by = "accession")
    
    # Set up factors  
    dat_snp_training$accession <- factor(dat_snp_training$accession)
    dat_snp_training$section_block <- factor(dat_snp_training$section_block)
    dat_snp_training$section <- factor(dat_snp_training$section)
    dat_snp_training$locality <- factor(dat_snp_training$locality)
    
    # Scale PCA axes
    dat_snp_training$PC1_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC1_gen)) / sd(dat_snp_training$PC1_gen)
    dat_snp_training$PC2_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC2_gen)) / sd(dat_snp_training$PC2_gen)
    dat_snp_training$PC3_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC3_gen)) / sd(dat_snp_training$PC3_gen)
    dat_snp_training$PC4_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC4_gen)) / sd(dat_snp_training$PC4_gen)
    dat_snp_training$PC5_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC5_gen)) / sd(dat_snp_training$PC5_gen)
    
    # Save unscaled version  
    dat_snp_training_unscaled <- dat_snp_training
    
    summary(dat_snp_training[, 1050:1060])    
    
    
    
    
    
    # Save means and sds for prediction
    scaled_snps_means_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                                   function(x) mean(x, na.rm = TRUE))
    scaled_snps_sds_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                                 function(x) sd(x, na.rm = TRUE))
    
    scaled_snps_means_training <- apply(dat_snp_training[, snp_col_names], MARGIN = 2,
                                        function(x) mean(x, na.rm = TRUE))
    scaled_snps_sds_training <- apply(dat_snp_training[, snp_col_names], MARGIN = 2,
                                      function(x) sd(x, na.rm = TRUE))
    
    scaled_snps_means_testing <- apply(dat_snp_testing[, snp_col_names], MARGIN = 2,
                                       function(x) mean(x, na.rm = TRUE))
    scaled_snps_sds_testing <- apply(dat_snp_testing[, snp_col_names], MARGIN = 2,
                                     function(x) sd(x, na.rm = TRUE))
    
  
    ## Scale genotypes
    dat_snp_all[, snp_col_names] <-   apply(dat_snp_all[, 
                                                snp_col_names], MARGIN = 2,
                                        function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    dat_snp_training[, snp_col_names] <-   apply(dat_snp_training[, 
                                                snp_col_names], MARGIN = 2,
                                        function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    dat_snp_testing[, snp_col_names] <-   apply(dat_snp_testing[, 
                                                                  snp_col_names], MARGIN = 2,
                                                 function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    
  
    
# Set a prediction data frame for training dataset
    

    # Set up each climate variable individually
    pred <-  expand.grid(height_2014 = 0,
                         accession = "1",
                         section_block = "Block1_1",
                         locality = "FHL",
                         PC1_clim =0, PC2_clim = 0,
                         PC1_gen = 0, PC2_gen = 0)
    
    for(var in "tmax_sum_dif"){
      
      var_unscaled <- paste0(var, "_unscaled")
      
      # Set scaled climate var dif, with a 0 and XXXX degree increase as well
      var_df <- expand.grid(accession = "1",
                            var_temp = c(seq(min(dat_all_clim[ ,var]),
                                             max(dat_all_clim[ ,var]),
                                             length.out = 150), 0))
      
      # Rename columns and set unscaled climate_var dif
      var_df <- var_df %>%
        rename(!!var_unscaled := var_temp) %>%
        mutate(!!var := forward_transform(x = get(var_unscaled),
                                          var = var,
                                          means = scaled_var_means_gbs_training,
                                          sds = scaled_var_sds_gbs_training))
      
      # Do a left join if the first variable, or else just bind cols, removing accession column  
      if(nrow(pred) == 1){
        pred <- left_join(pred, var_df)
      } else {
        pred <- bind_cols(pred, var_df[, -1])
      }
    }
    
    head(pred)
    
    pred
    
    
  
## Create folds for 80-20 validation
  set.seed(129)
  
 # folds <- caret::createFolds(factor(dat_snp_unscaled$accession), k = 5)
  folds <- NA

  # table(dat_snp_unscaled$accession)
  # table(dat_snp_unscaled$accession[folds[[1]]])
  # table(dat_snp_unscaled$accession[folds[[2]]])
  # table(dat_snp_unscaled$accession[folds[[3]]])
  # table(dat_snp_unscaled$accession[folds[[4]]])
  # table(dat_snp_unscaled$accession[folds[[5]]])
  # 
  # table(dat_snp_unscaled$section_block)
  # table(dat_snp_unscaled$section_block[folds[[1]]])
  # table(dat_snp_unscaled$section_block[folds[[2]]])
  # table(dat_snp_unscaled$section_block[folds[[3]]])
  # table(dat_snp_unscaled$section_block[folds[[4]]])
  # table(dat_snp_unscaled$section_block[folds[[5]]])
  
  
  ### Run full model and save out residuals for cluster analysis
  
  ## Original used in PNAS submission
  fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + s(tmax_sum_dif, bs=\"cr\") + s(accession, bs = \"re\") + s(locality, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\")  + s(tmax_sum_dif, by = PC1_gen, bs=\"cr\") + s(tmax_sum_dif, by = PC2_gen, bs=\"cr\") + s(PC1_clim, bs=\"cr\") + s(PC2_clim, bs=\"cr\")  + s(tmax_sum_dif, by = PC1_clim, bs=\"cr\") + s(tmax_sum_dif, by = PC2_clim, bs=\"cr\")"))
  
  
  # Testing out removing family effects
  fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\") + s(tmax_sum_dif, bs=\"cr\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\")"))
  
  
  # Gam with interaction effect
  gam_snp_all = bam(formula = formula(fixed_effects),
                    data = dat_snp_all,
                    discrete = TRUE, 
                    nthreads = 8,
                    method = "fREML",
                    family = "tw",
                    control = list(trace = FALSE))
  
  summary(gam_snp_all)
  
  dat_snp_all$rgr_resids <- resid(gam_snp_all)
  hist(dat_snp_all$rgr_resids, breaks = 50)
  
  ## Join residuals to training and testing data
  dat_snp_training <- dplyr::left_join(dat_snp_training, dplyr::select(dat_snp_all, accession_progeny, rgr_resids))
  dat_snp_training_unscaled <- dplyr::left_join(dat_snp_training_unscaled,
                                                dplyr::select(dat_snp_all, accession_progeny, rgr_resids))
  
  dat_snp_testing <- dplyr::left_join(dat_snp_testing, dplyr::select(dat_snp_all, accession_progeny, rgr_resids))
  dat_snp_testing_unscaled <- dplyr::left_join(dat_snp_testing_unscaled,
                                                dplyr::select(dat_snp_all, accession_progeny, rgr_resids))
  

# Save data to file that will be uploaded to cluster  
  save(dat_snp_all, dat_snp_all_unscaled,
       dat_snp_training, dat_snp_training_unscaled,
       dat_snp_testing, dat_snp_testing_unscaled,
       snp_col_names,
       pred,
       scaled_var_means_gbs_all, scaled_var_sds_gbs_all,
       scaled_var_means_gbs_testing, scaled_var_sds_gbs_testing,
       scaled_var_means_gbs_training, scaled_var_sds_gbs_training,
       scaled_snps_means_all, scaled_snps_sds_all,
       scaled_snps_means_training, scaled_snps_sds_training,
       scaled_snps_means_testing, scaled_snps_sds_testing,
       folds,
    file = paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))


 