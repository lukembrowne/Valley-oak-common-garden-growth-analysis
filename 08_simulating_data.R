

# Simulate data -------------------------------------------------------

# This script will generate a simulated dataset of relative growth rates of progeny planted across Tmax differences
# This script is designed to be run through the command line, which will vary the simulation parameters and allow the simulations to be run in batches

# Load libraries
library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv
library(tidyverse)
library(visreg)
library(rrBLUP)


# Read in command line arguments
args<-commandArgs(trailingOnly=TRUE)

print(args)

n_sites <- as.numeric(args[1])
qtl_ve <- as.numeric(args[2])
rep <- as.numeric(args[3])

set.seed(rep) # Set seed based on rep

# Experiment parameters
n_accession <- 300
n_blocks <- 5 # This equals the number of seedlings / accession per site

# QTL parameters
r.snp <- 0.99   # Proportion of variance in trait explained by SNPs.
d <- qtl_ve    # Proportion of additive genetic variance due to QTLs.
n <- n_accession    # Number of samples.
p <- 12000      # Number of markers (SNPs).
i <- seq(from = 1, to = p, length.out = 100)  # Indices of the QTLs ("causal variants").

# Tmax optimum function parameters
tmax_opt_sd <- 4

# Set up main dataframe
sim_dat <- expand.grid(site = as.character(1:n_sites),
            accession = as.character(1:n_accession),
            block = as.character(1:(n_blocks)))
dim(sim_dat)

sim_dat$section_block <- paste0(sim_dat$site, "_", sim_dat$block)
sim_dat$section_block <- factor(sim_dat$section_block)

sim_dat$ID <- as.character(1:nrow(sim_dat))

# table(sim_dat$accession, sim_dat$section_block)
# table(sim_dat$accession, sim_dat$site)
table(sim_dat$site)



# Generate effect sizes for ...
  
  # Block
    block_effect <- rnorm(n = levels(sim_dat$section_block), sd = .25)
    names(block_effect) <- levels(sim_dat$section_block)
    block_effect
  
  # Site
    site_effect <- rnorm(n = n_sites, sd = .25)
    names(site_effect) <- levels(sim_dat$site_effect)
    site_effect
  
  # accession
    accession_effect <- rnorm(n = n_accession, sd = .25)
    names(accession_effect) <- levels(sim_dat$accession)
   # accession_effect
  
  
# Simulate QTL - which is the optimum temperature for each genotype, determined by many loci

    # Small script to simulated data that models the an idealized
    # genome-wide association study, in which all single nucleotide
    # polymorphisms (SNPs) have alternative alleles at high frequencies in
    # the population (>5%), and all SNPs contribute some (small) amount of
    # variation to the quantitative trait. The SNPs are divided into two
    # sets: those that contribute a very small amount of variation, and
    # those that contribute a larger amount of variable (these are the
    # "QTLs"---quantitative trait loci).
    #
    # Adjust the five parameters below (r, d, n, p, i) to control how the
    # data are enerated.
    #
    # Peter Carbonetto
    # University of Chicago
    # January 23, 2017
    #
    
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
  # Note: this code only works if d or r are not exactly 0 or exactly 1.
  # Snps
  st.snp   <- c(r.snp/(1-r.snp) * d/var(X.snp %*% beta.snp))
  beta.snp <- sqrt(st.snp) * beta.snp
  sa.snp   <- max(Re(polyroot(c(c(var(X.snp %*% beta.snp) - r.snp/(1-r.snp)),
                            2*sum((X.snp %*% beta.snp) * (X.snp %*% u.snp))/n,
                            c(var(X.snp %*% u.snp))))))^2
  u.snp    <- sqrt(sa.snp) * u.snp

  # Set polygenic effects to 0 if QTL variation explained is really high
  if(qtl_ve == 0.99){
    u.snp = rep(0, p)
  }

  if(qtl_ve == 0.01){
    beta.snp = rep(0,p)
  }

  # Generate the quantitative trait measurements.
  y.total <- c(X.snp %*% (u.snp + beta.snp) + rnorm(n))
  
  # Check that the data are simulated in the way we requested.
  var.snp <- c(var(X.snp %*% (u.snp + beta.snp))/var(y.total))
  var.poly <- c(var(X.snp %*% (u.snp))/var(y.total))
  var.qtl <- c(var(X.snp %*% (beta.snp))/var(y.total))
  cat("Proportion of variance in Y explained by both polygenic and QTL effects =  ",var.snp,"\n")
  cat("Proportion of variance in Y explained by polygenic effects =  ",var.poly,"\n")
  cat("Proportion of variance in Y explained by QTL effects =  ",var.qtl,"\n")

  # Print out summaries
  summary(beta.snp)
  summary(u.snp)
  

  
# Generate tmax differences
  site_temp <- seq(25, 35, length.out = n_sites)
  
  # Set temp of home site for accession - assuming here no correlation between level of adaptation and home temperature
  accession_temp <- rnorm(n = n_accession, mean = 30, sd = 4)
  names(accession_temp) <- levels(sim_dat$accession)
  hist(accession_temp)
  summary(accession_temp)
  
  # Calculate tmax dif
  for(i in 1:nrow(sim_dat)){
    sim_dat$tmax_dif[i] <- site_temp[sim_dat$site[i]] -  accession_temp[sim_dat$accession[i]] 
  }
  

  # Figure for supplementary section
  hist(sim_dat$tmax_dif, breaks = 50, col = "steelblue2", xlab = "Tmax difference",
       las = 1, main = "Histogram of simulated Tmax differences for 2 planting sites")
  
  
# Use simulated phenotype as 'optimum' tmax_difs for growth response
  # Where growth rates of their progeny should peak
  tmax_opt <- scales::rescale(y.total,
                              to = c(-7.5, 2.5)) # This sets maximum and min of how far off optimum genotypes are
  names(tmax_opt) <- levels(sim_dat$accession)
  
  
  
# Save effect sizes of tmax dif
  
  ## accession variation based on QTL
  # Simulate where progeny are on transfer function for each accession, assuming transfer function comes from a normal distribution
  for(i in 1:nrow(sim_dat)){
    sim_dat$tmax_dif_effect[i] <- dnorm(sim_dat$tmax_dif[i], 
                                  mean = tmax_opt[sim_dat$accession[i]],
                                  sd = tmax_opt_sd)
    }
  
  # Rescale effect sizes
  tmax_dif_effect_mean <- mean(sim_dat$tmax_dif_effect)
  tmax_dif_effect_sd <- sd(sim_dat$tmax_dif_effect)
  sim_dat$tmax_dif_effect_scaled <- ((sim_dat$tmax_dif_effect - tmax_dif_effect_mean) / 
                                      tmax_dif_effect_sd )
  
  # Divide by XXX to reduce overall effect size
  sim_dat$tmax_dif_effect_scaled <- sim_dat$tmax_dif_effect_scaled / 5
  summary(sim_dat$tmax_dif_effect_scaled)
  
  
  
  
  
# Plot transfer function for each accession and calculate true breeding values
  tmax_seq <- seq(min(sim_dat$tmax_dif), max(sim_dat$tmax_dif), by = .1)
  
  # Make dataframe to save 'true' GEBV
  acc_prs <- data.frame(accession = NA, 
                        prs_true = NA) 
  
 # plot_n <- sample(names(tmax_opt), size = 5)
  plot_n <- names(tmax_opt)
  already_plotted = FALSE
  x = 1
  for(acc in names(tmax_opt)){
    
    acc_preds <-  dnorm(tmax_seq, 
                        mean = tmax_opt[acc],
                        sd = tmax_opt_sd)
  
    acc_prs[x, "accession"] <- acc
    acc_prs[x, "prs_true"] <- mean(acc_preds[tmax_seq > 0])

    rand_color <- sample(colors(), 1)
    
    ## Plot lines
    if(acc %in% plot_n && !already_plotted){
      plot(tmax_seq, acc_preds, type = "l", col = rand_color, lwd =0.25,
           las = 1, xlab = "Tmax difference", ylab = "Marginal effect on RGR")
      already_plotted = TRUE
    } else if (acc %in% plot_n && already_plotted) {
      lines(tmax_seq, acc_preds, type = "l", col = rand_color, lwd =0.25)
    }
    
    # Plot points where progeny lie
    if(acc %in% plot_n){
      # sub_progeny <- sim_dat[sim_dat$accession == acc, ]
      # points(sub_progeny$tmax_dif, sub_progeny$tmax_dif_effect, pch = 19, col = rand_color)
    }
    
    x = x + 1
  }
    
    # Plot average species-wide response
    lines(tmax_seq, dnorm(tmax_seq, 
                          mean = mean(tmax_opt),
                          sd = tmax_opt_sd), lwd = 3)
    abline(v = 0, lwd = 2, lty = 2)
    
  
    summary(acc_prs$prs_true)
    
  
## Figure for supplement of how 'true' breeding value was calculated
  fam1 <- dnorm(tmax_seq, 
                mean = -5,
                sd = tmax_opt_sd)
  fam2 <- dnorm(tmax_seq, 
                mean = 1,
                sd = tmax_opt_sd)
  plot(tmax_seq, fam1, type = "l", col = "blue", lty = 2,
       las = 1, xlab = "Tmax difference", ylab = "Marginal effect on RGR")
  # Bold lines used in calculation
  lines(tmax_seq[tmax_seq > 0], fam1[tmax_seq > 0], col = "blue", lwd = 2)
  # Add second family
  lines(tmax_seq, fam2, type = "l", col = "red", lty = 2)
  lines(tmax_seq[tmax_seq > 0], fam2[tmax_seq > 0], col = "red", lwd = 2)
  abline(v = 0, lwd = 2, lty = 1)
  
  legend("topright", legend = c(paste0("BV = ", round(mean(fam1[tmax_seq > 0]), 3)),
                                paste0("BV = ", round(mean(fam2[tmax_seq > 0]), 3))),
                                lty = 1, 
         col = c("blue", "red"), lwd = 2)
  
  
  
# Simulate growth rates
  for(i in 1:nrow(sim_dat)){
    
    sim_dat$rgr[i] <- 1 +
                     block_effect[sim_dat$section_block[i]] + 
                     site_effect[sim_dat$site[i]] + 
                     accession_effect[sim_dat$accession[i]] + 
                     sim_dat$tmax_dif_effect_scaled[i] + 
                     rnorm(n = 1, mean = 0, sd = .25)
  }  
  
  
  hist(sim_dat$rgr)
  

## Gam version
  gam <- bam(rgr ~ section_block + s(accession, bs = "re") +  s(tmax_dif, bs = "cr"),
             data = sim_dat,
             discrete = FALSE,
             nthreads = 8,
             method = "fREML",
             # accession = "gaussian",
          #   accession = "tw",
             control = list(trace = FALSE))

  summary(gam)

 # visreg(gam, partial = F)
  


# Compare model estimate to truth for supplementary material
  par(mfrow = c(1, 2))
  plot(tmax_seq,dnorm(tmax_seq, 
                      mean = mean(tmax_opt),
                      sd = tmax_opt_sd), type = "l", col = "black", lwd = 3,
       las = 1, xlab = "Tmax difference", ylab = "Marginal effect on RGR",
       main = "Simulated species-level response")
  abline(v = 0, lwd = 2, lty = 2)
  visreg(gam, partial = F, xvar = "tmax_dif", ylab = "Relative Growth Rate", rug = F,
         xlab = "Tmax difference", main = "GAM estimated species-level response")
  abline(v = 0, lwd = 1.5)
  
  par(mfrow = c(1, 1))
  

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
  # summary(pca_gen)
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
  
  # Scale PCA axes
  dat_snp_all$PC1_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC1_gen)) / sd(dat_snp_all$PC1_gen)
  dat_snp_all$PC2_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC2_gen)) / sd(dat_snp_all$PC2_gen)
  dat_snp_all$PC3_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC3_gen)) / sd(dat_snp_all$PC3_gen)
  dat_snp_all$PC4_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC4_gen)) / sd(dat_snp_all$PC4_gen)
  dat_snp_all$PC5_gen <- (dat_snp_all$PC1_gen - mean(dat_snp_all$PC5_gen)) / sd(dat_snp_all$PC5_gen)
  
  # Save unscaled version  
  dat_snp_all_unscaled <- dat_snp_all
  
 # summary(dat_snp_all[, 500:510])
  
  
# Save means and sds for prediction
  scaled_snps_means_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                                 function(x) mean(x, na.rm = TRUE))
  scaled_snps_sds_all <- apply(dat_snp_all[, snp_col_names], MARGIN = 2,
                               function(x) sd(x, na.rm = TRUE))

  ## Scale genotypes
  dat_snp_all[, snp_col_names] <-   apply(dat_snp_all[,
                                                      snp_col_names], MARGIN = 2,
                                          function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))


  
# Set up prediction dataframe
  pred <-  expand.grid(height_2014 = 0,
                       accession = "1",
                       section_block = "1_1",
                       locality = "FHL",
                       PC1_clim =0, PC2_clim = 0,
                       PC1_gen = 0, PC2_gen = 0)
  
  pred <- left_join(pred, expand.grid(accession = "1",
                              tmax_dif = c(seq(-10, 10, length.out = 150), 0)))
  
  pred$tmax_dif_unscaled <- pred$tmax_dif
  
  head(pred)
  
  folds <- NA
  
  ### Run full model and save out residuals for cluster analysis
  fixed_effects <- paste0(paste0("rgr ~ section_block  + s(tmax_dif, bs=\"cr\") + s(accession, bs = \"re\")"))
  
  # Gam with interaction effect
  gam_snp_all = bam(formula = formula(fixed_effects),
                    data = dat_snp_all,
                    discrete = FALSE, 
                    nthreads = 8,
                    method = "fREML",
                    control = list(trace = FALSE))
  
  summary(gam_snp_all)
  
  dat_snp_all$rgr_resids <- resid(gam_snp_all)
  dat_snp_all_unscaled$rgr_resids <- resid(gam_snp_all)
  hist(dat_snp_all$rgr_resids, breaks = 50)
  
  
 # Save data to file that will be uploaded to cluster  
  # Save entire environment
  save.image(paste0("./gam_cluster_sim_", n_sites, "sites_workspace.Rdata"))
  

   

