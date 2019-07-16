## Simulate data
library(tidyverse)
library(rrBLUP)
library(gamm4)
library(visreg)

## See if we have more power to detect patterns if we do it based on accession or based on genotype
## Simulate data based on genotype and tmax transfer distance
## Then compare to only have 2 samples (one in each garden) and looking for associations that way

## See if the way we are calculating Breeding value correlates with observed breeding values

## What is the false positive rate?


## Simulate QTL and rescale the phenotype to be the 'optimum' temp, which is the mean in the normal distribution

# Set parameters

set.seed(129)


n_accession <- 300
n_seedlings_per_accession <- 10
n_sites <- 2
n_blocks <- 10


sim_dat <- expand.grid(
            accession = as.character(1:n_accession),
            block = as.character(1:(n_blocks/2)))
dim(sim_dat)


sim_dat <- rbind(cbind(site = "1", sim_dat),
                 cbind(site = "2", sim_dat))
dim(sim_dat)

sim_dat$section_block <- paste0(sim_dat$site, "_", sim_dat$block)
sim_dat$section_block <- factor(sim_dat$section_block)

sim_dat$ID <- as.character(1:nrow(sim_dat))




# Generate effect sizes for ...

# Block
  block_effect <- rnorm(n = levels(sim_dat$section_block), sd = .25)
  names(block_effect) <- levels(sim_dat$section_block)
  block_effect

# Site
#  site_effect <- rnorm(n = n_sites, sd = .25)
  site_effect <- c(.5, -.5)
  names(site_effect) <- levels(sim_dat$site_effect)
  site_effect

# accession
  accession_effect <- rnorm(n = n_accession, sd = .25)
  names(accession_effect) <- levels(sim_dat$accession)
  accession_effect
  
  
# Simulate QTL - which is the optimum temperature for each genotype, determined by many loci
  
  # SCRIPT PARAMETERS
  r.snp <- 0.9   # Proportion of variance in trait explained by SNPs.
  d <- 0.0001    # Proportion of additive genetic variance due to QTLs.
  n <- n_accession    # Number of samples.
  p <- 2000        # Number of markers (SNPs).
  i <- seq(from = 1, to = p, length.out = 10)  # Indices of the QTLs ("causal variants").
  
  # GENERATE DATA SET
  # Generate the minor allele frequencies. They are uniformly
  # distributed on [0.05,0.5].
  maf <- 0.05 + 0.45 * runif(p)
  
  # Simulate genotype data X from an idealized population according to the 
  # minor allele frequencies generated above.
  X.snp <- (runif(n*p) < maf) +
    (runif(n*p) < maf)
  X.snp <- matrix(as.double(X.snp),n,p,byrow = TRUE)
  X.snp.unscaled <- X.snp
  
  # Center the columns of X.
  rep.row <- function (x, n)
    matrix(x,n,length(x),byrow = TRUE)
  X.snp <- X.snp - rep.row(colMeans(X.snp),n)

  # Generate (small) polygenic additive effects for the SNPs.
  u.snp <- rnorm(p)
 
  # Generate (large) QTL effects for the SNPs.
  beta.snp    <- rep(0,p)
  beta.snp[i] <- rnorm(length(i))

  
  # Adjust the additive effects so that we control for the proportion of
  # additive genetic variance that is due to QTL effects (d) and the
  # total proportion of variance explained (r). That is, we adjust beta
  # and u so that
  #
  #   r = a/(a+1)
  #   d = b/a,
  #
  # where I've defined
  #
  #   a = (u + beta)'*cov(X)*(u + beta),
  #   b = beta'*cov(X)*beta.
  #
  # Note: this code only works if d or r are not exactly 0 or exactly 1.
  # Snps
  st.snp   <- c(r.snp/(1-r.snp) * d/var(X.snp %*% beta.snp))
  beta.snp <- sqrt(st.snp) * beta.snp
  sa.snp   <- max(Re(polyroot(c(c(var(X.snp %*% beta.snp) - r.snp/(1-r.snp)),
                            2*sum((X.snp %*% beta.snp) * (X.snp %*% u.snp))/n,
                            c(var(X.snp %*% u.snp))))))^2
  u.snp    <- sqrt(sa.snp) * u.snp
  
  # Generate the quantitative trait measurements.
  y.total <- c(X.snp %*% (u.snp + beta.snp) + rnorm(n))
  
  # Check that the data are simulated in the way we requested.
  var.snp <- c(var(X.snp %*% (u.snp + beta.snp))/var(y.total))
  cat("Proportion of variance in Y explained by genotypes =  ",var.snp,"\n")
  

  ## Output in GEMMA format
  ### Output SNP information in BIMBAM format
  snp.sim <- as.data.frame(X.snp.unscaled) %>%
    `colnames<-`(paste0("snp_",
                        str_pad(as.character(1:ncol(X.snp.unscaled)),
                                4, pad = "0"))) %>%
    mutate(ID = paste0("ind_",
                       str_pad(as.character(1:nrow(X.snp.unscaled)), 4, pad = "0"))) %>%
    gather(key = "site", value = "value", -ID) %>%
    spread(key = ID, value = value) %>%
    mutate(allele_type1 = "C", allele_type2 = "T") %>%
    select(site, allele_type1, allele_type2, everything())
  
  head(snp.sim[, 1:20])
  
  dim(snp.sim)
  
  # Write to file
  # write_tsv(snp.sim,
  #           "./data_simulation/snp.sim.txt",
  #           col_names = FALSE)
  
  ## Make SNP annotation file
  snp.sim.annotation <- data.frame(snp.name = snp.sim$site,
                                   snp.position = str_split(snp.sim$site,
                                                            "_", simplify = T)[,2],
                                   snp.chrom = "1")
  
  head(snp.sim.annotation)
  
  # write_tsv(x = snp.sim.annotation,
  #           path = "./data_simulation/snp.sim.annotation.txt",
  #           col_names = FALSE)
  
  
  # Use simulated phenotype as 'optimum' temps for growth response
  tmax_opt <- scales::rescale(y.total, to = c(-7.5, 2.5))
  names(tmax_opt) <- levels(sim_dat$accession)
  
  
# Generate tmax differences
  site_temp <- c(35, 25) # Set temperatures of sites
  names(site_temp) <- c("1", "2")
  
  accession_temp <- rnorm(n = n_accession, mean = 30, sd = 2.5)
  names(accession_temp) <- levels(sim_dat$accession)
  hist(accession_temp)
  summary(accession_temp)
  
  # Calculate tmax dif
  for(i in 1:nrow(sim_dat)){
    sim_dat$tmax_dif[i] <- site_temp[sim_dat$site[i]] -  accession_temp[sim_dat$accession[i]] 
  }
  
  hist(sim_dat$tmax_dif)
  
  # Save effect sizes of tmax dif
  
  ## No accession variation based on QTL
  # sim_dat$tmax_dif_effect <- dnorm(sim_dat$tmax_dif,
  #                          mean = -5, sd = 4)
  # # Rescale
  # sim_dat$tmax_dif_effect <- scales::rescale( sim_dat$tmax_dif_effect,
  #                                             to = c(-.25, .25))
  # 
  # plot(sim_dat$tmax_dif, sim_dat$tmax_dif_effect)

  
  ## accession variation based on QTL
  for(i in 1:nrow(sim_dat)){
    sim_dat$tmax_dif_effect[i] <- dnorm(sim_dat$tmax_dif[i], 
                                        mean = tmax_opt[sim_dat$accession[i]],
                                        sd = 4)
    }
  
  # Rescale
  sim_dat$tmax_dif_effect <- scales::rescale(sim_dat$tmax_dif_effect,
                                              to = c(-.25, .25))
  
  plot(sim_dat$tmax_dif, sim_dat$tmax_dif_effect)
  
  
  sim_dat %>%
    filter(accession %in% as.character(1:10)) %>%
  ggplot(. ,aes(x = tmax_dif, y = tmax_dif_effect, col = accession)) + 
    geom_point(show.legend = F) + geom_line()
  
  
# Simulate growth rates
  
for(i in 1:nrow(sim_dat)){
  
  sim_dat$rgr[i] <- 1 +
                   block_effect[sim_dat$section_block[i]] + 
                   site_effect[sim_dat$site[i]] + 
                  accession_effect[sim_dat$accession[i]] + 
                  sim_dat$tmax_dif_effect[i] + 
                  rnorm(n = 1, mean = 0, sd = .25)
  
}  
  
  
  hist(sim_dat$rgr)
  
  
  
 #  
 #  ## Gam version
  gam <- bam(rgr ~ section_block + s(accession, bs = "re") +  s(tmax_dif, bs = "cr"),
             data = sim_dat,
             discrete = TRUE,
             nthreads = 8,
             method = "fREML",
             # accession = "gaussian",
          #   accession = "tw",
             control = list(trace = FALSE))

  summary(gam)

  visreg(gam, partial = F)
   
  

# Choose training and testing sets ----------------------------------------

  # 70/30 sampling
  training_moms <- sample(levels(sim_dat$accession),
                          size = length(levels(sim_dat$accession))*0.70)
  
  testing_moms <- levels(sim_dat$accession)[!levels(sim_dat$accession)
                                                     %in% training_moms]
  
  any(training_moms %in% testing_moms) # Should both be false
  any(testing_moms %in% training_moms)
  
  
  ## Save two datasets for training and testing
  sim_dat_training <- sim_dat %>%
    dplyr::filter(accession %in% training_moms)
  sim_dat_testing <- sim_dat %>%
    dplyr::filter(accession %in% testing_moms)
  
  dim(sim_dat_training)
  dim(sim_dat_testing)
  
  table(sim_dat_testing$section_block)
 # table(sim_dat_testing$locality)
  table(sim_dat_testing$accession)
  
  table(sim_dat_training$section_block)
 # table(sim_dat_training$locality)
  table(sim_dat_training$accession)  
  
  
  
  
  
  

# Set up model for cluster ------------------------------------------------

  # function to transform from raw to scaled variables
  forward_transform <- function(x, var, means, sds){
    ( x - means[var]) / sds[var]
  }
  
  # Function to back transform variables
  back_transform <- function(x, var, means, sds){
    (x * sds[var]) + means[var]
  }
  
  
  
  # Getting genetic data
  colnames(X.snp.unscaled) <- paste0("snp_", 1:ncol(X.snp.unscaled))
  snp_col_names <- colnames(X.snp.unscaled)

  ## Calculate kinship matrix
  
  # Need to first convert to -1, 0, 1 matrix
  kin_mat <- as.matrix(X.snp)
  kin_mat[kin_mat == 0] <- -1
  kin_mat[kin_mat == 1] <- 0
  kin_mat[kin_mat == 2] <- 1
  
  kin_mat <- rrBLUP::A.mat(kin_mat) # Calculate kinship matrix
  
  
  ## Calculate PCA on kinship matrix
  pca_gen = prcomp(kin_mat, center = TRUE, scale = TRUE) ## PCA on kinship matrix
  
  # biplot(pca_gen)
  summary(pca_gen)
  #  screeplot(pca_gen, bstick = TRUE)
  
  ### Figure S3 - Plotting results of PCA for supplementary material
  

  # Scale values for plotting in biplot - from: https://stats.stackexchange.com/questions/66926/what-are-the-four-axes-on-pca-biplot
  scores<-pca_gen$x
  lam <- pca_gen$sdev[1:3]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  lam <- lam^1
  yy<-t(t(pca_gen$rotation[, 1:3]) * lam)
  xx<-t(t(scores[, 1:3])/lam)

  # Pc1 vs pc2
    plot(xx[, 1], xx[, 2],
         las = 1,
         xlab = paste0("PC1: ", summary(pca_gen)$importance[2, 1] * 100, "% variance explained"),
         ylab = paste0("PC2: ", summary(pca_gen)$importance[2, 2] * 100, "% variance explained"),
         pch = 19,
         type = "n")
    text(xx[, 1], xx[, 2],
         labels = 1:nrow(X.snp), cex = .5)

  pcs <- as.data.frame(pca_gen$x[, 1:10])
  colnames(pcs) <- paste0(colnames(pcs), "_gen")
  
  snp_col_names <- colnames(X.snp.unscaled)
  gen_dat <- bind_cols(as.data.frame(X.snp.unscaled), pcs)
  gen_dat$accession <- factor(names(tmax_opt))
  
  
  # Set up dataframe with SNP data, height data, and climate data
  dat_snp_all = left_join(sim_dat,
                gen_dat,  by = "accession")
  
  # Set up factors  
  dat_snp_all$accession <- factor(dat_snp_all$accession)
  dat_snp_all$section_block <- factor(dat_snp_all$section_block)
  dat_snp_all$section <- factor(dat_snp_all$section)
#  dat_snp_all$locality <- factor(dat_snp_all$locality)
  
  # Scale PCA axes
  dat_snp_all$PC1_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC1_gen)) / sd(dat_snp_all$PC1_gen)
  dat_snp_all$PC2_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC2_gen)) / sd(dat_snp_all$PC2_gen)
  dat_snp_all$PC3_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC3_gen)) / sd(dat_snp_all$PC3_gen)
  dat_snp_all$PC4_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC4_gen)) / sd(dat_snp_all$PC4_gen)
  dat_snp_all$PC5_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC5_gen)) / sd(dat_snp_all$PC5_gen)
  
  # Save unscaled version  
  dat_snp_all_unscaled <- dat_snp_all
  
  summary(dat_snp_all[, 500:510])
  
  # ## Testing set  
  dat_snp_testing = left_join(sim_dat_testing,
              gen_dat,  by = "accession")

  # Set up factors
  dat_snp_testing$accession <- factor(dat_snp_testing$accession)
  dat_snp_testing$section_block <- factor(dat_snp_testing$section_block)
  dat_snp_testing$section <- factor(dat_snp_testing$section)
#  dat_snp_testing$locality <- factor(dat_snp_testing$locality)

  # Scale PCA axes
  dat_snp_testing$PC1_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC1_gen)) / sd(dat_snp_testing$PC1_gen)
  dat_snp_testing$PC2_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC2_gen)) / sd(dat_snp_testing$PC2_gen)
  dat_snp_testing$PC3_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC3_gen)) / sd(dat_snp_testing$PC3_gen)
  dat_snp_testing$PC4_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC4_gen)) / sd(dat_snp_testing$PC4_gen)
  dat_snp_testing$PC5_gen <- (dat_snp_testing$PC1_gen - mean(dat_snp_testing$PC5_gen)) / sd(dat_snp_testing$PC5_gen)

  # Save unscaled version
  dat_snp_testing_unscaled <- dat_snp_testing

  summary(dat_snp_testing[, 550:560])


  ## Training set
  dat_snp_training = left_join(sim_dat_training,
                               gen_dat,  by = "accession")

  # Set up factors
  dat_snp_training$accession <- factor(dat_snp_training$accession)
  dat_snp_training$section_block <- factor(dat_snp_training$section_block)
  dat_snp_training$section <- factor(dat_snp_training$section)
  #dat_snp_training$locality <- factor(dat_snp_training$locality)

  # Scale PCA axes
  dat_snp_training$PC1_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC1_gen)) / sd(dat_snp_training$PC1_gen)
  dat_snp_training$PC2_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC2_gen)) / sd(dat_snp_training$PC2_gen)
  dat_snp_training$PC3_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC3_gen)) / sd(dat_snp_training$PC3_gen)
  dat_snp_training$PC4_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC4_gen)) / sd(dat_snp_training$PC4_gen)
  dat_snp_training$PC5_gen <- (dat_snp_training$PC1_gen - mean(dat_snp_training$PC5_gen)) / sd(dat_snp_training$PC5_gen)

  # Save unscaled version
  dat_snp_training_unscaled <- dat_snp_training

  summary(dat_snp_training[, 550:560])

  
  
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
                       section_block = "1_1",
                       locality = "FHL",
                       PC1_clim =0, PC2_clim = 0,
                       PC1_gen = 0, PC2_gen = 0)
  
  pred <- left_join(pred, expand.grid(accession = "1",
                              tmax_dif_unscaled = c(seq(min(sim_dat$tmax_dif),
                                               max(sim_dat$tmax_dif),
                                               length.out = 150), 0)))
  
  pred$tmax_dif <- scale(pred$tmax_dif_unscaled)
  
  head(pred)
  
  pred
  
  folds <- NA
  
  
  ### Run full model and save out residuals for cluster analysis
  fixed_effects <- paste0(paste0("rgr ~ section_block  + s(tmax_dif, bs=\"cr\") + s(accession, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\")  + s(tmax_dif, by = PC1_gen, bs=\"cr\") + s(tmax_dif, by = PC2_gen, bs=\"cr\")"))
  
  
  # Gam with interaction effect
  gam_snp_all = bam(formula = formula(fixed_effects),
                    data = dat_snp_all,
                    discrete = TRUE, 
                    nthreads = 8,
                    method = "fREML",
               #    family = "tw",
                    control = list(trace = FALSE))
  
  summary(gam_snp_all)
  
  dat_snp_all$rgr_resids <- resid(gam_snp_all)
  hist(dat_snp_all$rgr_resids, breaks = 50)
  
  ## Join residuals to training and testing data
  dat_snp_training <- dplyr::left_join(dat_snp_training, 
                                       dplyr::select(dat_snp_all, ID, rgr_resids))
  dat_snp_training_unscaled <- dplyr::left_join(dat_snp_training_unscaled,
                                                dplyr::select(dat_snp_all, ID, rgr_resids))
  
  dat_snp_testing <- dplyr::left_join(dat_snp_testing, dplyr::select(dat_snp_all, ID, rgr_resids))
  dat_snp_testing_unscaled <- dplyr::left_join(dat_snp_testing_unscaled,
                                               dplyr::select(dat_snp_all, ID, rgr_resids))
  
  
  # Save data to file that will be uploaded to cluster  
  save(dat_snp_all, dat_snp_all_unscaled,
       dat_snp_training, dat_snp_training_unscaled,
       dat_snp_testing, dat_snp_testing_unscaled,
       snp_col_names,
       pred,
       sim_dat, sim_dat_testing, sim_dat_training,
       accession_effect,
       beta.snp, u.snp, block_effect, site_effect, tmax_opt,
       X.snp, X.snp.unscaled,
    #   scaled_var_means_gbs_all, scaled_var_sds_gbs_all,
     #  scaled_var_means_gbs_testing, scaled_var_sds_gbs_testing,
     #  scaled_var_means_gbs_training, scaled_var_sds_gbs_training,
       scaled_snps_means_all, scaled_snps_sds_all,
       scaled_snps_means_training, scaled_snps_sds_training,
       scaled_snps_means_testing, scaled_snps_sds_testing,
       folds,
    file = paste0("./data_simulation//gam_cluster_sim_", Sys.Date(), ".Rdata"))
  
  
  
  
  
  

# Analyzing cluster output ------------------------------------------------

  # Path to results
  path <- "run_381428_tmax_sum_dif_sim_test2"
  
  path_to_summaries <- paste0("./data_simulation/", path, "/model_summaries/")
  path_to_predictions <- paste0("./data_simulation/", path, "/model_predictions/")
  
  # Load data
  load(file = paste0("./data_simulation/", path,
                     "/gam_cluster_sim_2019-07-15.Rdata" ))
  
  
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
  base0 <- mean(pred_df_sub$pred[pred_df_sub$tmax_dif_unscaled == 0])
  
  pred_df_long =  pred_df_sub %>%
    dplyr::group_by(snp, genotype) %>%
    do(data.frame(height_change_warmer = mean((.$pred[.$tmax_dif_unscaled > 0]  -
                                                 .$pred[.$tmax_dif_unscaled == 0] ) /
                                                abs(.$pred[.$tmax_dif_unscaled == 0]) ) * 100,
                  height_change_warmer_base0 = mean((.$pred[.$tmax_dif_unscaled > 0]
                                                     - base0 ) /  abs(base0)) * 100,
                  height_change_warmer_absolute = mean(.$pred[.$tmax_dif_unscaled > 0])))
  
  head(pred_df_long)
  dim(pred_df_long)
  
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
  
  fdr_fvals = fdrtool(c(sum_df$f_val_gen_int),
                      statistic = "normal", plot = FALSE)
  
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
  
  # Local false discovery rate
  hist(fdr_fvals$lfdr, breaks = 40)
  sum(fdr_fvals$lfdr < 0.05)
  summary(fdr_fvals$lfdr)
  
  
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
  ## Using iterative filtering
  snp_pos <-data.frame(chrom = "1", pos = 1:ncol(X.snp))
  snp_pos
  

  man_df1 <- data.frame(SNP = snp_col_names, snp_pos)
  
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
  
  man_df$color <- "steelblue2"
  
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
  
  

  # Calculate polygenic risk score for each maternal genotype
  gen_dat_moms <- dplyr::select(gen_dat, accession, snp_col_names)
  dim(gen_dat_moms)
  
  # Choose SNPs
  n_snps_values <- c(1000, seq(2000, 11000, by = 2000))
  
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
    #  dat_snp_testing_unscaled$accession <- as.numeric(as.character(dat_snp_testing_unscaled$accession))
      test2 <- dplyr::left_join(dat_snp_testing_unscaled, 
                                dplyr::select(gen_dat_moms, accession, prs_by_count), by = "accession")
      
      # Scale prs
      test2$prs_scaled <- (test2$prs_by_count - mean(test2$prs_by_count)) / sd(test2$prs_by_count)
      
      #  test2$prs_scaled <- rnorm(n = length(test2$prs_scaled))
      
      # Set up model terms
      fixed_effects_resids <- paste0(paste0("rgr_resids ~ s(prs_scaled, bs=\"cr\") + 
                                            s(tmax_dif, by = prs_scaled, bs=\"cr\")"))
      
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
      visreg(gam_test_resids, xvar = "tmax_dif", by = "prs_scaled", partial = F, overlay = T, main = title)
      
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
  gen_dat_moms$prs_best <- gen_dat_moms$prs_10000
  
  
  prs_best_mean <- mean(gen_dat_moms$prs_best)
  prs_best_sd <- sd(gen_dat_moms$prs_best)
  gen_dat_moms$prs_best_unscaled <- gen_dat_moms$prs_best
  gen_dat_moms$prs_best_scaled <- (gen_dat_moms$prs_best - prs_best_mean) / prs_best_sd
  
  
  summary(gen_dat_moms$prs_best_scaled); hist(gen_dat_moms$prs_best_scaled, breaks = 50)
  
  
  
  
  
  
# Compare estimated to ‘true’ breeding value ------------------------------

  tmax_seq <- seq(-10, 10, by = .1)
  
  
  acc_prs <- data.frame(accession = NA, 
                        prs_true = NA)
  x = 1
  
  for(acc in names(tmax_opt)){
    
   acc_preds <-  dnorm(tmax_seq, 
          mean = tmax_opt[acc],
          sd = 4)
   
   
   acc_prs[x, "accession"] <- acc
   acc_prs[x, "prs_true"] <- mean(acc_preds[tmax_seq > 0])
   
   x = x + 1
   
   ## Plot lines
   if(x == 1){
     plot(tmax_seq, acc_preds, type = "l", col = sample(colors(), 1), las = 1)
   } else {
     lines(tmax_seq, acc_preds, type = "l", col = sample(colors(), 1))
   }
   
  }
  
  acc_prs
  
  summary(acc_prs$prs_true)
  hist(acc_prs$prs, breaks = 50)
  
  
  plot(scale(acc_prs$prs_true), gen_dat_moms$prs_best_scaled, pch = 19)
  abline(a = 0, b = 1, lwd = 2)
  cor.test(acc_prs$prs_true, gen_dat_moms$prs_best)
  

  

# Plot predictions of PRS -------------------------------------------------
  
  
  ## Calculate number of top and bottom SNPs 
  dat_snp_prs <- dat_snp_all %>%
  #  dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
    dplyr::select(accession, site, tmax_dif, 
                  rgr, rgr_resids, 
                  section_block) %>%
                  #PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
    dplyr::left_join(., dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled)) %>%
    dplyr::mutate(accession = factor(accession))
  
  head(dat_snp_prs)
  dim(dat_snp_prs)
  
  summary(dat_snp_prs$prs_best)
  summary(dat_snp_prs$prs_best_scaled)
  
  dat_snp_prs_unscaled <- dat_snp_prs
  
  # Set formula
  fixed_effects <- paste0(paste0("rgr_resids ~ s(prs_best_scaled, bs=\"cr\") + 
                                 s(tmax_dif, by = prs_best_scaled, bs=\"cr\")"))
  
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
  newdata <-  expand.grid(section_block = "1_1",
                          locality = "FHL",
                          height_2014 = 0,
                          accession = "6",
                          tmax_dif = c(seq(min(dat_snp_all_unscaled$tmax_dif),
                                               max(dat_snp_all_unscaled$tmax_dif),
                                               length.out = 150), 
                                           forward_transform(4.8, "tmax_dif",
                                               setNames(mean(dat_snp_all_unscaled$tmax_dif),
                                                        "tmax_dif"),
                                               setNames(sd(dat_snp_all_unscaled$tmax_dif),
                                                        "tmax_dif"))),
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
    # dplyr::mutate(tmax_dif_unscaled = back_transform(tmax_dif, var = "tmax_dif",
    #                                  means = setNames(mean(dat_snp_all_unscaled$tmax_dif),
    #                                                                       "tmax_dif"),
    #                                    sds = setNames(sd(dat_snp_all_unscaled$tmax_dif),
    #                                                                     "tmax_dif"))) %>%
    ggplot(., aes(x = tmax_dif, y = pred)) +
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
  
  
  # ### Run full model and save out residuals for cluster analysis
  # fixed_effects <- paste0(paste0("rgr ~ section_block + s(height_2014, bs=\"cr\")  + s(accession, bs = \"re\") + s(locality, bs = \"re\") + s(PC1_gen, bs=\"cr\") + s(PC2_gen, bs=\"cr\") + s(PC1_clim, bs=\"cr\") + s(PC2_clim, bs=\"cr\")"))
  # 
  # # Gam with interaction effect
  # gam_snp_all2 = bam(formula = formula(fixed_effects),
  #                    data = dat_snp_all,
  #                    discrete = TRUE, 
  #                    nthreads = 8,
  #                    method = "fREML",
  #                    family = "tw",
  #                    control = list(trace = FALSE))
  # 
  # summary(gam_snp_all2)
  # 
  # dat_snp_all2 <- dat_snp_all
  # 
  # dat_snp_all2$rgr_resids <- resid(gam_snp_all2)
  # hist(dat_snp_all2$rgr_resids, breaks = 50)
  # 
  # ## Mean center the residuals
  # dat_snp_all2$rgr_resids <-  dat_snp_all2$rgr_resids - mean(dat_snp_all2$rgr_resids)
  # summary(dat_snp_all2$rgr_resids)
  # 
  # 
  # dat_snp_prs2 <- dat_snp_all2 %>%
  #   dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
  #   dplyr::select(accession, site, locality, tmax_sum, tmax_sum_dif, 
  #                 rgr, rgr_resids, 
  #                 section_block, PC1_gen, PC2_gen, PC1_clim, PC2_clim, height_2014) %>%
  #   dplyr::left_join(., dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled)) %>%
  #   dplyr::mutate(accession = factor(accession)) %>%
  #   dplyr::mutate(tmax_sum_dif_unscaled = back_transform(tmax_sum_dif, var = "tmax_sum_dif",
  #                                                        means = scaled_var_means_gbs_all,
  #                                                        sds = scaled_var_sds_gbs_all))
  # 
  # # Set parameters
  # set.seed(1)
  # sample_size = 50
  # n_reps <- 10000
  # 
  # # Indices for top and bottom quartiles
  # upper_ind <- which(dat_snp_prs2$prs_best_scaled >= 1 & 
  #                      dat_snp_prs2$tmax_sum_dif_unscaled > 0)
  # length(upper_ind)
  # 
  # lower_ind <- which(dat_snp_prs2$prs_best_scaled <= 
  #                      quantile(dat_snp_prs2$prs_best_scaled, .25) 
  #                    & dat_snp_prs2$tmax_sum_dif_unscaled > 0)
  # length(lower_ind)
  # 
  # random_ind <- which(dat_snp_prs2$tmax_sum_dif_unscaled > 0)
  # length(random_ind)
  # 
  # climate_match <- dat_snp_prs2  %>%
  #   dplyr::filter(tmax_sum_dif_unscaled <= .5 & tmax_sum_dif_unscaled >= -.5)  
  # dim(climate_match)
  # 
  # climate_colder <- dat_snp_prs2  %>%
  #   dplyr::filter(tmax_sum_dif_unscaled <= -2)  
  # dim(climate_colder)
  # 
  # # Initialize output
  # output_df <- tibble(rep = 1:n_reps,
  #                     `Upper quartile GEBV` = NA,
  #                     `Lower quartile GEBV` = NA,
  #                     `None` = NA,
  #                     `Climate Match` = NA,
  #                     `Climate Colder` = NA)
  # 
  # 
  # for(x in 1:n_reps){
  #   
  #   if(x %% 500 == 0){
  #     cat("Working on x: ", x, " .... \n")
  #   }
  #   
  #   upper_sample <- sample(upper_ind, size = sample_size)
  #   lower_sample <- sample(lower_ind, size = sample_size)
  #   avg_sample <- sample(random_ind, size = sample_size)
  #   climate_match_sample <- sample(1:nrow(climate_match), size = sample_size)
  #   climate_colder_sample <- sample(1:nrow(climate_colder), size = sample_size)
  #   
  #   
  #   output_df$`Upper quartile GEBV`[x] <-     mean(dat_snp_prs2$rgr_resids[upper_sample])
  #   output_df$`Lower quartile GEBV`[x] <-     mean(dat_snp_prs2$rgr_resids[lower_sample])
  #   output_df$None[x] <-                      mean(dat_snp_prs2$rgr_resids[avg_sample])
  #   output_df$`Climate Match`[x] <-            mean(climate_match$rgr_resids[climate_match_sample])
  #   output_df$`Climate Colder`[x] <-          mean(climate_colder$rgr_resids[climate_colder_sample])
  #   
  # }
  # 
  # 
  # # Plot results 
  # output_df %>%
  #   gather(key = type, value = rgr_resid, -rep) %>%
  #   mutate(type = factor(type, levels = c("None", "Climate Match", "Climate Colder",
  #                                         "Upper quartile GEBV", "Lower quartile GEBV"))) %>%
  #   dplyr::filter(type != "Lower quartile GEBV") %>%
  #   ggplot(., aes(x = type, y = rgr_resid, fill = type)) + 
  #   geom_hline(yintercept = 0, lty = 2)+
  #   # geom_violin(draw_quantiles = c(.5)) +
  #   geom_boxplot(outlier.shape = NA) + 
  #   xlab("Choice criteria") +
  #   ylab("Adjusted Relative Growth Rate") + 
  #   scale_fill_manual(values =  c("#227eb8", "#fdc086", "#ffff99", "#7fbf7b")) + # "#af8dc3")) +
  #   scale_x_discrete(labels = c("Random", "Climate match \n -0.5° to 0.5°\nTmax transfer", 
  #                               "Climate match \n <= -2°\nTmax transfer",
  #                               "GEBV \n >= 1 SD")) +
  #   #"GEBV \n Lower quartile")) + 
  #   scale_y_continuous(breaks = seq(-.25, .25, by = 0.05),
  #                      limits = c(-.1, .15)) +
  #   theme_bw(8) + 
  #   theme(panel.border = element_blank(), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         axis.line = element_line(colour = "black"),
  #         axis.text.x = element_text(size = 8),
  #         axis.title.x = element_text(size = 10),
  #         axis.title.y = element_text(size = 10),
  #         legend.position = "none") +
  #   NULL
  # 
  # # Save output
  # # ggsave(filename = paste0("./figs_tables/fig2/Figure 2 - rgr predictions by choice ",
  # #                          Sys.Date(), ".pdf"),
  # #        units = "cm",
  # #        width = 9, height = 8)
  # 
  # 
  # 
  # ## Comparing different scenarios 
  # growth_by_scenario <-  output_df %>%
  #   gather(key = type, value = rgr_resid, -rep) %>%
  #   mutate(type = factor(type, levels = c("None", "Climate Match", "Climate Colder",
  #                                         "Upper quartile GEBV", "Lower quartile GEBV"))) %>%
  #   dplyr::filter(type != "Lower quartile GEBV") %>%
  #   group_by(type) %>%
  #   summarise_all(mean)
  # 
  # growth_by_scenario
  # 
  # growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "Climate Match"]
  # 
  # growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "Climate Colder"]
  # 
  # growth_by_scenario$rgr_resid[growth_by_scenario$type == "Upper quartile GEBV"] - growth_by_scenario$rgr_resid[growth_by_scenario$type == "None"]
  # 
  # 
  # 
  

  
  

