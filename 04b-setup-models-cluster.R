

# Load libraries ----------------------------------------------------------
  
  source("./01_clean-process-explore-data.R")
  source("./03_adding-climate-data.R")
  
  dim(dat_all_scaled)
  dim(dat_gbs_only_scaled)


  # install_github("jdstorey/qvalue")
  library(gamm4)
  library(visreg)
  library(rrBLUP)
  library(qvalue)
  library(MuMIn)
  library(vegan)
  library(beepr)
  library(patchwork)
  library(rasterVis)
  library(pdp)
  library(ggrepel)
  
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
  
  
  ## Calculate PCA on genetic data
    pca_gen = prcomp(bglr_gen_scaled, center = FALSE, scale = FALSE)
   # biplot(pca_gen)
    summary(pca_gen)
  #  screeplot(pca_gen, bstick = TRUE)
    
    pcs <- as.data.frame(pca_gen$x[, 1:10])
    colnames(pcs) <- paste0(colnames(pcs), "_gen")
    
    bglr_gen_scaled <- bind_cols(bglr_gen_scaled, pcs)

    
  # Set up dataframe with SNP data, height data, and climate data
    dat_snp = left_join(dat_gbs_only_scaled,
                      bind_cols(gen_dat_clim[, c("accession", snp_col_names)],
                                         pcs),
                      by = "accession")
  
   
  # Set up factors  
    dat_snp$accession <- factor(dat_snp$accession)
    dat_snp$section_block <- factor(dat_snp$section_block)
    dat_snp$section <- factor(dat_snp$section)
    
    
  # Scale PCA axes
    dat_snp$PC1_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC1_gen)) / sd(dat_snp$PC1_gen)
    dat_snp$PC2_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC2_gen)) / sd(dat_snp$PC2_gen)
    dat_snp$PC3_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC3_gen)) / sd(dat_snp$PC3_gen)
    dat_snp$PC4_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC4_gen)) / sd(dat_snp$PC4_gen)
    dat_snp$PC5_gen <- (dat_snp$PC1_gen - mean(dat_snp$PC5_gen)) / sd(dat_snp$PC5_gen)
    
    
  # Save unscaled version  
    dat_snp_unscaled <- dat_snp
 
    
  # Scale genotypes  
    scaled_snps_means <- apply(dat_snp[, snp_col_names], MARGIN = 2,
          function(x) mean(x, na.rm = TRUE))
    
    scaled_snps_sds <- apply(dat_snp[, snp_col_names], MARGIN = 2,
                               function(x) sd(x, na.rm = TRUE))
    
    dat_snp[, snp_col_names] <-   apply(dat_snp[, 
                                             snp_col_names], MARGIN = 2,
                    function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    
    dat_snp[1:10, 1050:1060]
    
    # Mean impute
    # dat_snp[, snp_col_names] <- apply(dat_snp[, 
    #                                                 snp_col_names], MARGIN = 2,
    #            function(x) ifelse(is.na(x), yes = 0, no = x))
    
    dat_snp[1:10, 1050:1060]
    
    summary(dat_snp[, 1050:1060])
    
    
    
# Set a prediction data frame


 # Prediction frame with a XXXX degree increase
    
  degree_increase = 4.8
    
  # Set up each climate variable individually
  pred <-  expand.grid(height_2014 = 0,
                      accession = "1",
                      section_block = "IFG_1",
                      PC1_gen = 0, PC2_gen = 0, PC3_gen = 0)
  
  for(var in climate_vars_dif){
    
    var_unscaled <- paste0(var, "_unscaled")
    
    
    
    # Set scaled climate var dif, with a 0 and XXXX degree increase as well
    var_df <- expand.grid(accession = "1",
                          var_temp = c(seq(min(dat_all_clim[ ,var]),
                                           max(dat_all_clim[ ,var]),
                                           length.out = 50), 0, degree_increase))
    
    # Rename columns and set unscaled climate_var dif
    var_df <- var_df %>%
      rename(!!var_unscaled := var_temp) %>%
      mutate(!!var := forward_transform(x = get(var_unscaled),
                                              var = var,
                                              means = scaled_var_means_gbs_only,
                                              sds = scaled_var_sds_gbs_only))
    
   # Do a left join if the first variable, or else just bind cols, removing accession column  
    if(nrow(pred) == 1){
    pred <- left_join(pred, var_df)
    } else {
      pred <- bind_cols(pred, var_df[, -1])
    }
  }
  
  head(pred)

  pred

# Save data to file that will be uploaded to cluster  
  # save(dat_snp, dat_snp_unscaled,
  #      snp_col_names, pred,
  #      scaled_snps_means, scaled_snps_sds,
  #      scaled_var_means_gbs_only, scaled_var_sds_gbs_only,
  #   file = paste0("./output/gam_cluster_", Sys.Date(), ".Rdata"))


 