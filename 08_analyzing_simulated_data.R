
# Analyzing cluster output ------------------------------------------------
library(rrBLUP)
library(visreg)
library(psych)
library(fdrtool)
library(qqman)
library(tidyverse)
library(readr)
library(mgcv, lib.loc = "/u/home/l/lukembro/R/x86_64-pc-linux-gnu-library/3.5") # Make sure to load latest version of mgcv


  # function to transform from raw to scaled variables
  forward_transform <- function(x, var, means, sds){
    ( x - means[var]) / sds[var]
  }
  
  # Function to back transform variables
  back_transform <- function(x, var, means, sds){
    (x * sds[var]) + means[var]
  }
  



# Read in command line arguments
args<-commandArgs(trailingOnly=TRUE)

print(args)

n_sites <- as.character(args[3])


# Path to results
path <- as.character(args[1])

path_to_summaries <- paste0(path, "/model_summaries/")
path_to_predictions <- paste0(path, "/model_predictions/")

# Load data

data_file <- as.character(args[2])

load(file = paste0("./", data_file))


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
  do(data.frame(height_change_warmer_absolute = mean(.$pred[.$tmax_dif_unscaled > 0]),
                height_change_warmer_absolute_training = mean(.$pred_training[.$tmax_dif_unscaled > 0])))

head(pred_df_long)
dim(pred_df_long)

summary(pred_df_long$height_change_warmer_absolute)
summary(pred_df_long$height_change_warmer_absolute_training)


# Join prediction dataframe to summary dataframe
sum_df_long_w_preds <- left_join(sum_df_long,
                                 pred_df_long)


# Correct for multiple testing


# F values of interaction term
hist(sum_df$f_val_gen_int, breaks = 50)

# fdr_fvals = fdrtool(c(sum_df$f_val_gen_int),
#                     statistic = "normal", plot = FALSE)

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
# Plot figure
# Need to hand-edit the axis labesl for q values to say Q instead of P
png( paste0("./figures/Qqplots of p values ", n_sites, "sites_", Sys.Date(), ".png"),
         width = 7, height = 7, units = "in", res = 300)

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


dev.off()

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




# Figure 2 - Manhattan plot of Q values -----------------------
## Using iterative filtering
snp_pos <-data.frame(chrom = "1", pos = 1:ncol(X.snp))


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
ggsave(paste0("./figures/manhattan plot interaction term wo labels_", n_sites, "sites_",
              Sys.Date(), ".png"),
       width = 19, height = 4, units = "cm", dpi = 600)




# Calculate PRS -----------------------------------------------------------

# Calculate polygenic risk score for each maternal genotype
gen_dat_moms <- dplyr::select(gen_dat, accession, snp_col_names)
dim(gen_dat_moms)

# Choose SNPs
n_snps_values <- round(c(seq(100, nrow(sum_df), length.out = 8)))

# Initialize dataframe for output
prs_df_summary <- tibble(n_snps = n_snps_values,
                         r2_training = NA)

# Order snps based on lowest Q values
snp_list_ordered <- sum_df$snp[order(sum_df$q_val)]

# Initialize
gen_dat_moms$prs <- NULL
gen_dat_moms$prs_count <- NULL
gen_dat_moms$prs_training <- NULL
gen_dat_moms$prs_training_count <- NULL

x = 1

for(snp_name in snp_list_ordered){
  
  if(x %% 500 == 0){
    cat("Working on snp number: ", x, "... \n")
  }
  
  # Subset to genotypes
  sub <- gen_dat_moms[ , snp_name]
  
  preds_sub <- dplyr::filter(sum_df_long_w_preds, snp == snp_name) %>%
    dplyr::select(genotype, height_change_warmer_absolute, 
                  height_change_warmer_absolute_training, n_gen, q_val)
  
  counted <- preds_sub$genotype
  
  # Add in median value for missing genotypes
  if(!"0" %in% preds_sub$genotype){
    preds_sub <- rbind(preds_sub, data.frame(genotype = "0",
                                             height_change_warmer_absolute = 0,
                                             height_change_warmer_absolute_training = 0,
                                             n_gen = NA,
                                             q_val = NA))
  }
  if(!"1" %in% preds_sub$genotype){
    preds_sub <- rbind(preds_sub, data.frame(genotype = "1",
                                             height_change_warmer_absolute = 0,
                                             height_change_warmer_absolute_training = 0,
                                             n_gen = NA,
                                             q_val = NA))
  }
  if(!"2" %in% preds_sub$genotype){
    preds_sub <- rbind(preds_sub, data.frame(genotype = "2",
                                             height_change_warmer_absolute = 0,
                                             height_change_warmer_absolute_training = 0,
                                             n_gen = NA,
                                             q_val = NA))
  }
  
  prs <-  case_when(
    sub == "0" ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "0"],  
    sub == "1"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "1"],
    sub == "2"  ~ preds_sub$height_change_warmer_absolute[preds_sub$genotype == "2"],
    TRUE ~ 0
  )
  
  prs_training <- case_when(
    sub == "0" ~  preds_sub$height_change_warmer_absolute_training[preds_sub$genotype == "0"],  
    sub == "1"  ~ preds_sub$height_change_warmer_absolute_training[preds_sub$genotype == "1"],
    sub == "2"  ~ preds_sub$height_change_warmer_absolute_training[preds_sub$genotype == "2"],
    TRUE ~ 0
  )
  
  prs_count <- ifelse(as.data.frame(sub)[, 1] %in% counted, yes = 1, no = 0)
  prs_training_count <- ifelse(as.data.frame(sub)[, 1] %in% counted, yes = 1, no = 0)
  
  if(x == 1){
    gen_dat_moms$prs <- prs
    gen_dat_moms$prs_count <- prs_count
    gen_dat_moms$prs_training <- prs_training
    gen_dat_moms$prs_training_count <- prs_training_count
  } else {
    gen_dat_moms$prs <- gen_dat_moms$prs + prs
    gen_dat_moms$prs_count <- gen_dat_moms$prs_count + prs_count
    gen_dat_moms$prs_training <- gen_dat_moms$prs_training + prs_training
    gen_dat_moms$prs_training_count <- gen_dat_moms$prs_training_count + prs_training_count
  }
  
  if(x %in% n_snps_values){
    cat("\n||--------------------------------||\n")
    cat("Running model for:", x, " SNPS... \n")
    
    # Adjust for missing data
    gen_dat_moms$prs_by_count <- gen_dat_moms$prs / gen_dat_moms$prs_count
    gen_dat_moms$prs_training_by_count <- gen_dat_moms$prs_training / gen_dat_moms$prs_training_count
    
    # Save prs score as its own variable for later use
    gen_dat_moms[, paste0("prs_", x)] <-  gen_dat_moms$prs_by_count
    gen_dat_moms[, paste0("prs_training_", x)] <-  gen_dat_moms$prs_training_by_count
    
    ### Join to testing dataframe
    #  dat_snp_testing_unscaled$accession <- as.numeric(as.character(dat_snp_testing_unscaled$accession))
    test2 <- dplyr::left_join(dat_snp_testing_unscaled, 
              dplyr::select(gen_dat_moms, accession, prs_training_by_count), by = "accession")
    
    # Scale prs
    test2$prs_training_scaled <- (test2$prs_training_by_count - mean(test2$prs_training_by_count)) / sd(test2$prs_by_count)
    

    # Set up model terms
    fixed_effects_resids <- paste0(paste0("rgr_resids ~ s(prs_training_scaled, bs=\"cr\") + 
                                          s(tmax_dif, by = prs_training_scaled, bs=\"cr\")"))
    
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
    visreg(gam_test_resids, xvar = "tmax_dif", by = "prs_training_scaled", partial = F,
           overlay = T, main = title)
    
    cat(paste(title, "\n"))
    
    ## Save results   
    prs_df_summary$r2_training[prs_df_summary$n_snps == x] <- summary(gam_test_resids)$r.sq
    
  } # End model loop
  
  x = x + 1 # Move onto next SNP
  
} # End snp loop


prs_df_summary

prs_df_summary %>%
  arrange(desc(r2_training))

## Plot relationship between # of SNPS used and variance explained in testing set
ggplot(prs_df_summary, aes(x = n_snps, y = r2_training)) + 
  geom_point(size = 2) + 
  geom_line() + 
  xlab("Number of SNPs in GEBV estimation") +
  #ylab("R2 adjusted in testing set") +
  ylab(bquote('R'^2[adj]~'training test')) +
  scale_x_continuous(breaks = prs_df_summary$n_snps) +
  # scale_y_continuous(breaks = seq(0, 0.02, by = 0.0025), limits = c(-0.0015, 0.015)) + 
  theme_bw(7) + 
  geom_hline(yintercept = 0, lty = 2) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = (margin(0, 0, 0, 0, "cm")),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  NULL


# # Save as file
ggsave(paste0("./figures/Testing variance explained by number of snps_", n_sites, "_sites_",
              Sys.Date(), ".png"),
       width = 10, height = 4, units = "cm", dpi = 300)



top_r2_index <- which(prs_df_summary$r2_training == max(prs_df_summary$r2_training,
                                               na.rm = T))

# Set top PRS
gen_dat_moms$prs_best <- gen_dat_moms[, paste0("prs_",
                                               prs_df_summary$n_snps[top_r2_index])]


prs_best_mean <- mean(gen_dat_moms$prs_best)
prs_best_sd <- sd(gen_dat_moms$prs_best)
gen_dat_moms$prs_best_unscaled <- gen_dat_moms$prs_best
gen_dat_moms$prs_best_scaled <- (gen_dat_moms$prs_best - prs_best_mean) / prs_best_sd


summary(gen_dat_moms$prs_best_scaled); hist(gen_dat_moms$prs_best_scaled, breaks = 50)








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

#  Save gam summary to file
sink(file = paste0("./figures/prs_all_gam_model_summary_",n_sites, "_sites_",
                   Sys.Date(), ".txt"))
summary(gam_prs_all)
anova(gam_prs_all)
sink()



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
#### Testing
dat_snp_testing2 <- dat_snp_testing %>%
  #  dplyr::mutate(accession = as.numeric(as.character(accession))) %>%
  dplyr::left_join(.,  dplyr::select(gen_dat_moms, accession, prs_best, prs_best_scaled))

dim(dat_snp_testing2)

# Gam with interaction effect
gam_prs_testing = bam(formula = formula(fixed_effects),
                      data = dat_snp_testing2,
                      discrete = TRUE,
                      nthreads = 8,
                      method = "fREML",
                      #   family = "tw",
                      control = list(trace = FALSE))

summary(gam_prs_testing)

## Save output

# Save gam summary to file
sink(file = paste0("./figures/prs_testing_gam_model_summary_",n_sites, "_sites_",
                   Sys.Date(), ".txt"))
summary(gam_prs_testing)
anova(gam_prs_testing)
sink()


## Plot predictions

# Set up dataframe for prediction
newdata <-  expand.grid(section_block = "1_1",
                        locality = "FHL",
                        height_2014 = 0,
                        accession = "1",
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
  #   scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5)) +
  #   scale_y_continuous(breaks = seq(0, 1.75, by = 0.25)) + 
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
ggsave(filename = paste0("./figures/prs growth ",n_sites, "_sites_", Sys.Date(), ".png"),
       units = "cm", width = 8, height = 6)





# Compare estimated to ‘true’ values ------------------------------


## Block effects - missing tmax effect
block_pred <- visreg(gam_snp_all, xvar = "section_block",
                     cond=list(tmax_dif=0))$fit

for(i in 1:nrow(block_pred)){
  site_effect_temp <- ifelse(grepl(pattern = "^1", 
                                   block_pred$section_block[i]),
                             yes = site_effect[1],
                             no = site_effect[2])
  
  block_pred$truth[i] <-1 +  block_effect[block_pred$section_block[i]] +
    site_effect_temp + 
    accession_effect[block_pred$accession[i]]
}

ggplot(block_pred, aes(section_block, visregFit)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = .1) + 
  geom_point(aes(section_block, truth), col = "red", size = 2) + 
  ylab("RGR") + 
  theme_bw(15)


# Accession predictions
accession_pred <- visreg(gam_snp_all, xvar = "accession",
                         cond=list(tmax_dif=0))$fit

for(i in 1:nrow(accession_pred)){
  site_effect_temp <- ifelse(grepl(pattern = "^1", 
                                   accession_pred$section_block[i]),
                             yes = site_effect[1],
                             no = site_effect[2])
  
  accession_pred$truth[i] <- 1 + 
    block_effect[accession_pred$section_block[i]] +
    site_effect_temp + 
    accession_effect[accession_pred$accession[i]]
}

ggplot(accession_pred, aes(accession, visregFit)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = .1) + 
  geom_point(aes(accession, truth), col = "red", size = 2) + 
  ylab("RGR") + 
  theme_bw(15)



# Plot true vs estimated breeding values  
summary(acc_prs$prs_true)
hist(acc_prs$prs, breaks = 50)

# Set colors for training and testing set
col <- ifelse(acc_prs$accession %in% training_moms, "steelblue", "coral")


png(paste0("./figures/GAM correlation with true GEBV ", 
                     n_sites, "_sites.png"), width = 8, height = 5, units = "in", res = 300)
plot(scale(acc_prs$prs_true), gen_dat_moms$prs_best_scaled, pch = 19,
     xlab = "GEBV True", ylab = "GEBV Estimated", col = col)
abline(a = 0, b = 1, lwd = 2)

cor_overall <- cor(acc_prs$prs_true, gen_dat_moms$prs_best_scaled)
cor_training <- cor(acc_prs$prs_true[acc_prs$accession %in% training_moms], 
                    gen_dat_moms$prs_best_scaled[gen_dat_moms$accession %in% training_moms])
cor_testing <- cor(acc_prs$prs_true[!acc_prs$accession %in% training_moms], 
                   gen_dat_moms$prs_best_scaled[!gen_dat_moms$accession %in% training_moms])


legend("bottomright", legend = c(paste0("R = ", round(cor_overall, 3)),
                                  paste0("R = ", round(cor_training, 3)),
                                  paste0("R = ", round(cor_testing, 3))),
        col = c("black", "steelblue", "coral"), pch = 19)

        
dev.off() 



## Only training moms
plot(scale(acc_prs$prs_true[acc_prs$accession %in% training_moms]),
     gen_dat_moms$prs_best_scaled[gen_dat_moms$accession %in% training_moms], pch = 19,
     xlab = "GEBV True", ylab = "GEBV Estimated")
abline(a = 0, b = 1, lwd = 2)
cor.test(acc_prs$prs_true[acc_prs$accession %in% training_moms], 
         gen_dat_moms$prs_best_scaled[gen_dat_moms$accession %in% training_moms])

## Only testing moms
plot(scale(acc_prs$prs_true[!acc_prs$accession %in% training_moms]),
     gen_dat_moms$prs_best_scaled[!gen_dat_moms$accession %in% training_moms], pch = 19,
     xlab = "GEBV True", ylab = "GEBV Estimated")
abline(a = 0, b = 1, lwd = 2)
cor.test(acc_prs$prs_true[!acc_prs$accession %in% training_moms], 
         gen_dat_moms$prs_best_scaled[!gen_dat_moms$accession %in% training_moms])


# Correlate p value with estimated effect size
snp_effect_df <- data.frame(snp = snp_col_names, 
                            snp_effect = u.snp + beta.snp)

sum_df_wtruth <- left_join(sum_df, snp_effect_df)

sum_df_wtruth <- sum_df_wtruth %>%
  mutate(snp_effect  = abs(snp_effect))

pairs.panels(select(sum_df_wtruth, snp_effect, dev_explained, rsq_adj,
                    p_val_gen_int, f_val_gen_int, p_val_gen_main,
                    f_val_gen_main, q_val))




# Testing out conventional approaches -------------------------------------

## Testing out different transfer functions per accession
gam_acc <- bam(rgr ~ section_block + s(tmax_dif, bs = "cr") + s(accession, bs = "re") +  s(tmax_dif, by = accession),
               data = sim_dat,
               discrete = TRUE,
               nthreads = 8,
               method = "fREML",
               # accession = "gaussian",
               #   accession = "tw",
               control = list(trace = FALSE))

# summary(gam_acc)

## Calculate RGR at tmax > 0 based on individual accession curves
# Make dataframe to save 'true' GEBV
acc_rgr_absolute <- data.frame(accession = NA, 
                               height_change_warmer = NA) 

plot = FALSE
x = 1
for(acc in names(tmax_opt)){
  
  cat("Working on accession:", acc, "... \n")
  
 # if(acc %in% testing_moms){
 #   acc_rgr_absolute[x, "accession"] <- acc
 #   acc_rgr_absolute[x, "height_change_warmer"] <- NA
 #   acc_rgr_absolute[x, "tmax_opt"] <- NA
 #   x = x + 1
 #   next
 # }
  
  # Build prediction dataframe
  pred_acc <-  expand.grid(height_2014 = 0,
                           accession = acc,
                           section_block = "1_1",
                           locality = "FHL",
                           PC1_clim =0, PC2_clim = 0,
                           PC1_gen = 0, PC2_gen = 0)
  
  pred_acc <- left_join(pred_acc, expand.grid(accession = acc,
                                              # tmax_dif= c(seq(min(sim_dat$tmax_dif),
                                              #                max(sim_dat$tmax_dif),
                                              #                length.out = 150))))
                                              
                                              tmax_dif= c(seq(-10,
                                                              10,
                                                              length.out = 150))))
  
  pred_acc$preds <- predict(gam_acc, newdata = pred_acc)
  
  ## Scaled plot
  if(acc %in% names(tmax_opt) && plot){
    plot(pred_acc$tmax_dif, 
         scales::rescale(pred_acc$preds, c(0, 1)), 
         type  = "l", lwd = 1.5, las = 1,
         ylab = "Scaled RGR effect", xlab = "Tmax dif",
         main = paste0("Accession: ", acc))
    abline(vline = 0, lty = 2)
    
    # Calculate 'true' transfer function
    truth <-  dnorm(pred_acc$tmax_dif,
                    mean = tmax_opt[acc],
                    sd = tmax_opt_sd)
    lines(pred_acc$tmax_dif, scales::rescale(truth, c(0, 1)), 
          lwd = 1.5, col = "blue")
    
    sub_progeny <- sim_dat[sim_dat$accession == acc, ]
    points(sub_progeny$tmax_dif, 
           scales::rescale(sub_progeny$tmax_dif_effect, c(0, 1)), 
           pch = 19, col = "blue")
  }
  
  acc_rgr_absolute[x, "accession"] <- acc
  acc_rgr_absolute[x, "height_change_warmer"] <- mean(pred_acc$preds[pred_acc$tmax_dif > 0])
  acc_rgr_absolute[x, "tmax_opt"] <- pred_acc$tmax_dif[which(pred_acc$preds == 
                                                               max(pred_acc$preds))]
  x = x + 1
}


acc_rgr_absolute

## Compare estimated and 'true' values of tmax_opt
# plot(tmax_opt, acc_rgr_absolute$tmax_opt, pch = 19)
# cor.test(tmax_opt, acc_rgr_absolute$tmax_opt)
# abline(a = 0, b = 1 , lwd = 2)

# Comapre with true values of gebv
plot(acc_prs$prs_true, acc_rgr_absolute$height_change_warmer, pch = 19)
cor.test(acc_prs$prs_true, acc_rgr_absolute$height_change_warmer)
abline(a = 0, b = 1 , lwd = 2)


plot(acc_rgr_absolute$tmax_opt, acc_rgr_absolute$height_change_warmer, pch = 19)
cor.test(acc_rgr_absolute$tmax_opt,
         acc_rgr_absolute$height_change_warmer)
abline(a = 0, b = 1 , lwd = 2)



### Calculate BLUPs based on BGLR
 library(BGLR)
 
 ind.names <- paste0("ind_", 1:length(y.total))
 
 # Computing the genomic relationship matrix for SNPs
 A.snp.temp <- scale(X.snp.unscaled,center=TRUE,scale=TRUE)
 A.snp <- tcrossprod(A.snp.temp)/ncol(A.snp.temp)
 
 
 # Fitting the model
 
 # Parameters    
 nIter <- 20000
 burnIn <- 10000
 verbose <- F
 
 y = acc_rgr_absolute$height_change_warmer
 summary(y)
 #  y = acc_rgr_absolute$tmax_opt
 
 ## Model with SNP and methylation data
 mod.snp <-BGLR( y = y,
                 ETA = list(list(K = A.snp, model = 'RKHS')),
                 nIter = nIter, burnIn = burnIn,
                 saveAt='./bglr/', 
                 verbose = verbose)
 
 str(mod.snp)
 
 # Variance explained
 1 - mod.snp$varE/var(y)
 
 
 ## Compare estimated breeding values
 bglr_gebv <- mod.snp$ETA[[1]]$u
 
 plot(scale(acc_prs$prs_true),  bglr_gebv, pch = 19,
      xlab = "GEBV True", ylab = "GEBV Estimated (BLUP)", las = 1)
 #  abline(a = 0, b = 1, lwd = 2)
 cor.test(acc_prs$prs_true, bglr_gebv)
 cor.test(acc_prs$prs_true, bglr_gebv, method = "spearman")
 
 
 
 
 # Set colors for training and testing set
 col <- ifelse(acc_prs$accession %in% training_moms, "steelblue", "coral")
 
 
 png(paste0("./figures/BLUP correlation with true GEBV ", 
                      n_sites, "_sites.png"), width = 8, height = 5, units = "in", res = 300)
 
 plot(scale(acc_prs$prs_true), scale(bglr_gebv), pch = 19,
      xlab = "GEBV True", ylab = "GEBV Estimated (BLUP)", col = col)
 abline(a = 0, b = 1, lwd = 2)
 
 cor_overall <- cor(acc_prs$prs_true, bglr_gebv)
 cor_training <- cor(acc_prs$prs_true[acc_prs$accession %in% training_moms], 
                     bglr_gebv[gen_dat_moms$accession %in% training_moms])
 cor_testing <- cor(acc_prs$prs_true[!acc_prs$accession %in% training_moms], 
                    bglr_gebv[!gen_dat_moms$accession %in% training_moms])
 
 
 legend("bottomright", legend = c(paste0("R = ", round(cor_overall, 3)),
                                  paste0("R = ", round(cor_training, 3)),
                                  paste0("R = ", round(cor_testing, 3))),
        col = c("black", "steelblue", "coral"), pch = 19)

 dev.off()     
 
 
 
 
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
 png(paste0("./figures/GAM and BLUP correlation with true GEBV ", 
                     n_sites, "_sites.png"), width = 10, height = 10, units = "in", res = 300)
 pairs.panels(data.frame(gebv_true = acc_prs$prs_true, 
                         gebv_blup_BGLR = bglr_gebv, 
                       #  gev_blup_GAPIT = gapit_out$BLUP,
                         gebv_gam = gen_dat_moms$prs_best_scaled),
              ellipses = F)
 dev.off()


  # Save entire environment
  save.image(paste0("./gam_cluster_sim_", n_sites, "sites_after_analysis.Rdata"))
  
 
