# TODO

# - reread Lasky 2015 to see if I'm missing anything
# - decide what type of Bayes model for markers
# - run BGLR models for real - longer run times and burn in, etc

# - make marker models where BV is predicted using only the top 100 (xxx) SNPs?

# - how to best present results of gam models
# - table with rows as models, columns with AIC and significance of smooth terms
# - choose to focus on models with highest AIC - best at explaining the response variable

# - run survival models?


# - Download more future climate variables



# Source files ------------------------------------------------------------

source("./01_clean-process-explore-data.R")
source("./03_adding-climate-data.R")


# Load libraries ----------------------------------------------------------

library(gamm4)
library(lmerTest)
library(caret)
library(usdm) # to check for vif
library(flashpcaR)
library(parallel)
library(cowsay)
library(visreg)

library(BGLR) # For genomic selection models - tutorials here: https://github.com/gdlc/BGLR-R
library(rrBLUP) # Calculate kinship matrix


# Genomic selection models to estimate breeding values  --------------------------------
  
  ## Notes about BGLR
  ## Can get 95% credible intervat from .dat files output to file but these include burn in samples
  ## Need to learn what RKHS output means - what us uStar parameter, etc?

  
# Formatting genetic data and calculating kinship matrix
  
  # Joining climate and genetic data
    gen_dat_clim <- left_join(climate_gbs_mom, dplyr::select(gen_dat, -accession), 
                              by = c("id" = "gbs_name"))
    dim(gen_dat_clim)

  # Getting genetic data
    snp_col_names <- colnames(gen_dat[, -c(1, 2)])
  
  # Scale data and impute missing before converting 012 to -101
    bglr_gen_scaled <- flashpcaR::scale2(gen_dat_clim[, snp_col_names])
    bglr_gen_scaled[1:10, 1:10]
  
  # Convert to -101 from 012
    bglr_gen101 <- gen_dat_clim[, snp_col_names]
    bglr_gen101[bglr_gen101 == 0] <- -1
    bglr_gen101[bglr_gen101 == 1] <- 0
    bglr_gen101[bglr_gen101 == 2] <- 1
    summary(bglr_gen101[, 1:10])
  
  # Calculate kinship matrix
    A <- rrBLUP::A.mat(X = bglr_gen101)
  
    
## Adding a random variable
  gen_dat_clim$random <- rnorm(n = nrow(gen_dat_clim))    
    
    
    
## Write a loop that goes through each climate variable and calculates marker and kinship models, calculates the CV accuracy, and saves output into usable formats    

mod_out <- list() # Save list that will hold model outputs and such
  
climate_vars <- c(climate_vars[1:2], "random")

for(var in climate_vars){
  
  # Initialize list
  mod_out[[var]] <- list()
  
  #cat("Working on:", var, "... \n")
  say(paste0("Working on: ", var, " ..."),  "trilobite")
  # Phenotypic response - cimate of origin
   y <- dplyr::pull(gen_dat_clim, var)
 
  # Setting the linear predictor
    # Options for model types - more info in manual
    # FIXED - flat prior
    # BL - Bayesian lasso
    # BRR - Bayesian ridge regression
    # BayesA - scaled-t prior
    # BL - Double-Exponential prior
    # BayesB - two component mixture prior with a point of mass at zero and a sclaed-t slab
    # BayesC - two component mixture prior with a point of mass at zero and a Gaussian slab
  
   ETA_kin <- list(list(K = A, model='RKHS')) # Kinship matrix
   
   ETA_marker <- list(list(X = bglr_gen_scaled, 
                            model='BayesC')) # Markers
   
   ETA_kinmarker <- list(list(K = A, model = "RKHS"), # Kinship matrix
                           list(X = bglr_gen_scaled, model = "BayesC")) # Markers
  
  ## Setting cross validation
    k_folds <- 10 # Set number of folds
    folds <- createFolds(y = 1:length(y), k = k_folds)
    
    ## Iterations include burn in
    n_iter <- 4000 
    n_burn <- 1000
    thin <- 5
    
  # Set output prefixes   
    prefix_kin <- paste0("./output/BGLR/", var, "_kin_")
    prefix_marker <- paste0("./output/BGLR/", var, "_marker_")  
    prefix_kinmarker <- paste0("./output/BGLR/", var, "_kinmarker_") 
    
  ## Cross validation in parallel  
    
    # First create a function that fits one fold
    fitFold <-  function(y, ETA, folds, fold){
     
      tst = folds[[fold]]
      yNA = y
      yNA[tst] = NA
      
      fm <- BGLR(y = yNA,
                 ETA = ETA,
                 nIter = n_iter,
                 burnIn = n_burn,
                 thin = thin,
                 verbose=F)
      
      cor_cv <- cor(y[tst], fm$yHat[tst])
      return(cor_cv)
    }
    
    # Kinship matrix
      cv_kin <- mclapply(FUN = fitFold, 
                    y = y,
                    ETA = ETA_kin,
                    folds = folds,
                    X = 1:length(folds),
                    mc.cores = 10)
      cv_kin <- unlist(cv_kin)
      
      # Save to list
      mod_out[[var]] <- append(mod_out[[var]], list("cv_kin" = cv_kin))
      
    
    # Marker scores
      cv_marker <- mclapply(FUN = fitFold, 
                            y = y,
                            ETA = ETA_marker,
                            folds = folds,
                            X = 1:length(folds),
                            mc.cores = 10)
      cv_marker <- unlist(cv_marker)
      
      # Save to list
      mod_out[[var]] <- append(mod_out[[var]], list("cv_marker" = cv_marker))
      
    # Combining kin & Marker scores
      cv_kinmarker <- mclapply(FUN = fitFold, 
                            y = y,
                            ETA = ETA_kinmarker,
                            folds = folds,
                            X = 1:length(folds),
                            mc.cores = 10)
      cv_kinmarker <- unlist(cv_kinmarker)
      
      # Save to list
      mod_out[[var]] <- append(mod_out[[var]], list("cv_kinmarker" = cv_kinmarker))
      
      
  ## Running full models
      
      ## Kinship model
      cat("Fitting full kinship model... \n")
      mod_kin <- BGLR(y =  y,  # The data-vector 
                      response_type = "gaussian",
                      ETA = ETA_kin,
                      nIter = n_iter, 
                      burnIn = n_burn,
                      thin = thin,
                      saveAt = prefix_kin, verbose=F)
      
      # Save to list
      mod_out[[var]] <- append(mod_out[[var]], list("mod_kin" = mod_kin))
      
      
      ## Marker model
      cat("Fitting full marker model... \n")
      mod_marker <- BGLR(y =  y,  # The data-vector 
                         response_type = "gaussian",
                         ETA = ETA_marker,
                         nIter = n_iter, 
                         burnIn = n_burn,
                         thin = thin,
                         saveAt = prefix_marker, verbose=F)
      
      mod_out[[var]] <- append(mod_out[[var]], list("mod_marker" = mod_marker))
      
      ## Kin-Marker model
      cat("Fitting full kinship & marker model... \n")
      mod_kinmarker <- BGLR(y =  y,  # The data-vector 
                         response_type = "gaussian",
                         ETA = ETA_kinmarker,
                         nIter = n_iter, 
                         burnIn = n_burn,
                         thin = thin,
                         saveAt = prefix_kinmarker, verbose=F)
      
      mod_out[[var]] <- append(mod_out[[var]], list("mod_kinmarker" = mod_kinmarker))
      
} # End climate variable loop   


# Save model output?
# save(mod_out, file = paste0("./output/mod_out_", Sys.Date(), ".Rdata"))

 load("./output/mod_out_2018-04-09.Rdata")

      
      
## Extract and compare CV scores
    
  # Kinship
    cv_kin <- as.data.frame(lapply(mod_out, "[[", "cv_kin"))
    # summary(cv_kin)
    colMeans(cv_kin)
    
  # Marker scores  
    cv_marker <- as.data.frame(lapply(mod_out, "[[", "cv_marker"))
    # summary(cv_marker)
    colMeans(cv_marker)
    
  # Kinmarker scores  
    cv_marker <- as.data.frame(lapply(mod_out, "[[", "cv_kinmarker"))
    # summary(cv_marker)
    colMeans(cv_marker)
    
  
# Check out trace plots for convergance
  samps <- list.files("./output/BGLR/", full.names = TRUE)
  for(samp in samps){
    tmp <- scan(samp)
    plot(tmp, type = "o", pch = 19, cex = 0.5,
         main = samp, xlab = "Iteration", xaxt = "n")
    axis(side = 1, at = seq(1, length(tmp), length.out = 10),
         labels =  round(seq(1, length(tmp), length.out = 10))*thin)
    abline(v = n_burn/thin, col =  "firebrick")
  }
    
# Extract model fits
  mod_fits_df <- data.frame(mod = "mod_kin", 
                            t(unlist(lapply(lapply(mod_out, "[[", "mod_kin"), 
                                            function(x) x$fit))))
  mod_fits_df <- rbind(mod_fits_df, 
                       data.frame(mod = "mod_marker", 
                      t(unlist(lapply(lapply(mod_out, "[[", "mod_marker"),
                                      function(x) x$fit)))))
  
  mod_fits_df <- rbind(mod_fits_df, 
                       data.frame(mod = "mod_kinmarker", 
                                  t(unlist(lapply(lapply(mod_out, "[[", "mod_kinmarker"),
                                                  function(x) x$fit)))))
  
  mod_fits_df
  
  
# Exploring marker effects for just marker model
  for(var in names(mod_out)){
    
    var_list <- mod_out[[var]]
    
    # Manhattan like plot
    plot(var_list$mod_marker$ETA[[1]]$b^2, 
         ylab = 'Estimated Squared-Marker Effect',
         pch = 19, cex = .5, col = "grey20", 
         main = paste0('Marker Effects:', var), las = 1,
         xlab = "Marker")
    
    # Plotting top SNPs
    top_snp <- names(sort(var_list$mod_marker$ETA[[1]]$b^2, decreasing = TRUE))[1]
    g = ggplot(gen_dat_clim, aes(y = pull(gen_dat_clim, top_snp), x = tmax_sum)) +
      ggtitle(paste0(var, ":", top_snp)) + ylab("Allele freq") +
      geom_jitter() + geom_smooth() + theme_bw()
    print(g)
  }  
  
  
# Calculating breeding values for just marker model
  for(var in names(mod_out)){

    # Matrix multiplaction - genotypes times marker effects 
     u <- as.matrix(bglr_gen_scaled) %*% mod_out[[var]]$mod_marker$ETA[[1]]$b
    
    # Save back into list
      mod_out[[var]]$mod_marker$ETA[[1]]$u <- u
  
    # Plot breeding values vs phenotype
      # Need to add in intercept (mu)
      plot(mod_out[[var]]$mod_marker$ETA[[1]]$u + mod_out[[var]]$mod_marker$mu, 
           mod_out[[var]]$mod_marker$y,
           col = "grey50", las = 1, main = paste0("Marker mod: ", var),
           ylab = var, xlab = "Breeding value", pch = 19) 
      abline(a=0,b=1,lwd=2)
  }
  

# Plot breeding values vs. phenotype for kinship models
  for(var in names(mod_out)){
    # Need to add in intercept (mu)
    plot(mod_out[[var]]$mod_kin$ETA[[1]]$u + mod_out[[var]]$mod_kin$mu, 
         mod_out[[var]]$mod_kin$y,
         col = "grey50", las = 1, main = paste0("Kinship mod: ", var),
         ylab = var, xlab = "Breeding value", pch = 19) 
    abline(a = 0, b = 1, lwd = 2)
  }
  
  
# Calculate breeding values for kin-marker models  
  for(var in names(mod_out)){
    
    # Matrix multiplaction - genotypes times marker effects 
    u <- as.matrix(bglr_gen_scaled) %*% mod_out[[var]]$mod_kinmarker$ETA[[2]]$b
    
    # Save back into list
    mod_out[[var]]$mod_kinmarker$ETA[[2]]$u <- u
    
    # Plot breeding values vs phenotype
    # Need to add in intercept (mu)
    plot(mod_out[[var]]$mod_kinmarker$ETA[[1]]$u + 
           mod_out[[var]]$mod_kinmarker$ETA[[2]]$u +
           mod_out[[var]]$mod_kinmarker$mu, 
           mod_out[[var]]$mod_kinmarker$y,
         col = "grey50", las = 1, main = paste0("Marker mod: ", var),
         ylab = var, xlab = "Breeding value", pch = 19) 
    abline(a=0,b=1,lwd=2)
  }
  
  
# Make a dataframe of breeding values that will match up growth data
  bvs <- data.frame(accession = gen_dat_clim$accession)
  
  for(var in names(mod_out)){
    bvs[ , paste0("bv.", var, ".kin")] <- c(mod_out[[var]]$mod_kin$ETA[[1]]$u)
    bvs[ , paste0("bv.", var, ".marker")] <- c(mod_out[[var]]$mod_marker$ETA[[1]]$u)
    # Calculating kinship & marker BV by summing them together
    bvs[ , paste0("bv.", var, ".kin.marker")] <- c(mod_out[[var]]$mod_kinmarker$ETA[[1]]$u + 
                                              mod_out[[var]]$mod_kinmarker$ETA[[2]]$u)
  }
  
  dim(bvs)
  summary(bvs)
  head(bvs)
  
  
  
# Correlating BVs with growth data ----------------------------------------
  
  dat_bv <- left_join(dat_all_scaled, bvs, by = "accession")
  dim(dat_bv)
  
  # Check for missing data in breeding values - there shouldn't be any
   dat_bv %>%
      filter(is.na(bv.tmax_sum.marker)) %>%
      dplyr::select(accession) %>%
      group_by(accession) %>%
      count()
   
   # Add in a random variable
   dat_bv$random <- rnorm(n = nrow(dat_bv))
   
   
   # Set up factors 
     dat_bv$section <- factor(dat_bv$section)
     dat_bv$section_block <- factor(dat_bv$section_block) 
     dat_bv$section_row <- factor(dat_bv$section_row)
     dat_bv$section_line <- factor(dat_bv$section_line)
     dat_bv$accession <- factor(dat_bv$accession)
     
  ## Check variance inflation factor
     vif(as.data.frame(dat_bv[, c("height_2014", "rgr", 
                                  "tmax_sum_dif", "bv.tmax_sum.kin")]))
     
# For each climate variable fit 5 game models
   # 1) Null model with just Height 2014
   # 2) Model with climate transfer
   # 3) Model with climate transfer x kinship breeding value
   # 4) Model with climate transfer x marker breeding value
   # 5) Model with climate transfer x kin_marker breeding value
   # 6) Model with climate transver x kinship interaction
   # 7) Model with climate transfer x marker interaction   
   # 8) Model with climate transfer x kin_marker interaction
     
 # Loop through climate variables
     
     gam_mods <- list()
     
     x = 1
     
     for(climate_var in climate_vars){
       
       say(paste0("Working on: ", climate_var, " ..."),  "trilobite")
       
       
       # Get variable names
       bv_var_kin <- paste0("bv.", climate_var, ".kin")
       bv_var_marker <- paste0("bv.", climate_var, ".marker")
       bv_var_kin_marker <- paste0("bv.", climate_var, ".kin.marker")
       
         if(climate_var != "random"){
           climate_var_dif <- paste0(climate_var, "_dif")
         } else {
           climate_var_dif = "random" # For random variable only
         }
       
       ## Set formulas for fixed effects / smoothed terms
       fixed_null <- "rgr ~ site + s(height_2014)"
       
       fixed_clim_dif <- paste0("rgr ~ site + s(height_2014) + s(", 
                                     climate_var_dif,")")
       
       fixed_bv_kin <- paste0("rgr ~ site + s(height_2014) + s(", 
                                      climate_var_dif,") +s(", bv_var_kin, ")")
       
       fixed_bv_marker <- paste0("rgr ~ site +  s(height_2014) + s(", 
                                         climate_var_dif,") +s(", bv_var_marker, ")")
       
       fixed_bv_kin_marker <-  paste0("rgr ~ site +  s(height_2014) + s(", 
                                              climate_var_dif,") +s(", bv_var_kin_marker, ")")
      ## Interaction models 
       fixed_bv_kin_int <- paste0("rgr ~ site + s(height_2014) + te(", 
                              climate_var_dif,") + te(", bv_var_kin, ") + ti(",
                              climate_var_dif, ", ", bv_var_kin, ")")
       fixed_bv_marker_int <- paste0("rgr ~ site + s(height_2014) + te(", 
                                  climate_var_dif,") + te(", bv_var_marker, ") + ti(",
                                  climate_var_dif, ", ", bv_var_marker, ")")
       fixed_bv_kin_marker_int <- paste0("rgr ~ site + s(height_2014) + te(", 
                                  climate_var_dif,") + te(", bv_var_kin_marker, ") + ti(",
                                  climate_var_dif, ", ", bv_var_kin_marker, ")")
       
       # Formula for random effects
       random <-  '+ s(section_block, bs="re")  + s(accession, bs = "re")'
       

     ## Convert fixed and random parts of formula to actual formula objetc
     # No idea why we need to call formula twice, but it works
     make_formula <- function(fixed, random){
       return(formula(formula(enquote(paste0(fixed, random)))))
     }

       # Fit models 
       gam_null <- bam(formula = make_formula(fixed_null, random),
                       data = dat_bv,
                       nthreads = 8,
                       method = "fREML")
       
       gam_clim_dif <- bam(formula = make_formula(fixed_clim_dif, random),
                        data = dat_bv,
                        nthreads = 8,
                        method = "fREML")
       
       gam_kin <- bam(formula = make_formula(fixed_bv_kin, random),
                      data = dat_bv,
                      nthreads = 8,
                      method = "fREML")
       
       gam_marker <- bam(formula = make_formula(fixed_bv_marker, random),
                         data = dat_bv,
                         nthreads = 8,
                         method = "fREML")
       
       gam_kin_marker <- bam(formula = make_formula(fixed_bv_kin_marker, random),
                             data = dat_bv,
                             nthreads = 8,
                             method = "fREML")
       ## Interaction models
       gam_kin_int <- bam(formula = make_formula(fixed_bv_kin_int, random),
                      data = dat_bv,
                      nthreads = 8,
                      method = "fREML")
       
       gam_marker_int <- bam(formula = make_formula(fixed_bv_marker_int, random),
                         data = dat_bv,
                         nthreads = 8,
                         method = "fREML")
       
       gam_kin_marker_int <- bam(formula = make_formula(fixed_bv_kin_marker_int, random),
                             data = dat_bv,
                             nthreads = 8,
                             method = "fREML")
       
       list_out <- list(list(gam_null = gam_null, 
                             gam_clim_dif = gam_clim_dif, 
                             gam_kin = gam_kin,
                             gam_marker = gam_marker,
                             gam_kin_marker = gam_kin_marker,
                             gam_kin_int = gam_kin_int,
                             gam_marker_int = gam_marker_int,
                             gam_kin_marker_int = gam_kin_marker_int))
       
       names(list_out) <- climate_var
       
       gam_mods[[x]] <- list_out
       
       x = x +1
     } # End climate vars loop 
     
     gam_mods <- unlist(gam_mods, recursive = FALSE)
     
 # Save gam mods
   save(gam_mods, file = paste0("./output/gam_mods_out_", Sys.Date(), ".Rdata"))
 
   
 
 # Extract AIC values from each model - the Mer part
   AIC_scores <- numeric()
   for(mod in names(gam_mods[[1]])){
     tmp <- unlist(lapply(lapply(gam_mods, "[[", mod), function(x) AIC(x)))
     names(tmp) <- paste0(mod, "_", names(tmp))
     AIC_scores <- c(AIC_scores, tmp)
   }
   
   AIC_df <- data.frame(climate_var = rep(names(gam_mods), times = length(gam_mods[[1]])),
                        model = rep(names(gam_mods[[1]]), each = length(climate_vars)),
                        AIC = AIC_scores)
   
   AIC_df
   
   AIC_df %>%
     arrange(climate_var, AIC)   
 
   
   
 # Extract P values for terms 
   mod_stats <- data.frame(NULL)
   for(var in names(gam_mods)){
     for(mod in names(gam_mods[[1]])){
       
       cat("Working on..", var, "...", mod, "...\n")
       
       # Make dataframe out of summary table for each gam mod
         sum <- data.frame(summary(gam_mods[[var]][[mod]])$s.table)
         sum$parameter <- rownames(sum) # Make column for parameter
         rownames(sum) <- NULL # remove row names
         
      # Bind to mod stats   
       mod_stats <- rbind(mod_stats, 
                           data.frame(climate_var = var, model = mod, sum))
     }
   }   
   
   mod_stats2 <- mod_stats %>%
            dplyr::select(climate_var, model, parameter, everything()) %>%
            dplyr::select(-Ref.df) %>%
            mutate(significance = ifelse(p.value <= 0.05, "***", "")) %>%
            mutate_at(vars(c("edf", "F", "p.value")), funs(round(., 3))) %>%
            left_join(., AIC_df, by = c("climate_var", "model"))

   View(mod_stats2)
   
# Output gam model results to table
   
   
   mod_stats2 <- mod_stats2 %>%
     mutate(parameter_general = case_when(
       grepl(pattern = "ti\\(", parameter) ~ "interaction",
       grepl(pattern = "_dif", parameter) ~ "climate_transfer",
       grepl(pattern = "random)", parameter) ~ "climate_transfer",
       grepl(pattern = "\\.kin)", parameter) ~ "kinship",
       grepl(pattern = "\\.kin\\.marker)", parameter) ~ "kin_marker",
       grepl(pattern = "\\.marker)", parameter) ~ "marker",
       grepl(pattern = "section_block", parameter) ~ "section_block",
       grepl(pattern = "accession", parameter) ~ "accession",
       TRUE                      ~  "random_effect")
     )
   
   View(mod_stats2)
   
  # Munging to get in right format 
  mod_stats3 <-  mod_stats2 %>% 
     dplyr::select(-edf, - F, -parameter, -significance) %>%
     spread(parameter_general, p.value) %>%
    mutate(AIC = round(AIC, 2)) %>% # Round AIC
    arrange(climate_var, AIC) # Order by lowest AIC score
  
  mod_stats3
  dim(mod_stats3)
  
  View(mod_stats3)
   
   ## SHould be 24 rows (one for each model)
   library(rtf)
  
  rtf <- RTF(file="./output/mod_stats.doc")
  addTable(rtf, mod_stats3)
  done(rtf)
   

    
  ## Print out model outputs
    pdf("./output/gam_mod_out.pdf")   
    for(var in names(gam_mods)){
      for(mod in names(gam_mods[[1]])){
        
        visreg::visreg(gam_mods[[var]][[mod]],
                       main = paste0(var, " - ", mod))
      }
    }    
    dev.off()
      

    
    
    
    
  

# Simulating data ---------------------------------------------------------
  
  # To see if able to pick up interaction terms, etc with only 2 garden sites  
    
  nSites <- 10  
    
  sim <- expand.grid(accession = 1:350, site = 1:nSites)  
  
  sim <- sim %>%
    left_join(., data.frame(accession = 1:350, 
                            tmax_sum = rnorm(n = 350, mean = 0, sd = 1)))
  
  tmax_sites <- seq(-3, 3, length.out = nSites)
  
  for(site in unique(sim$site)){
    sim$tmax_sum_site[sim$site == site] <- tmax_sites[site]
  }
  
  sim$tmax_sum_dif <- sim$tmax_sum_site - sim$tmax_sum
  
  
  plot(sim$tmax_sum_dif, sim$tmax_sum)
  
  ## simulate growth rates
  
  intercept = 0.5
  
  sim$rgr <- intercept + 
              #  .25  * sim$tmax_sum_dif + -.25 * sim$tmax_sum_dif^2 + # Polynomial transfer
                    -.5 * sim$tmax_sum_dif +
                    -.5 * sim$tmax_sum +  # Linear effect
                    -.15 *sim$tmax_sum * sim$tmax_sum_dif + # Interaction effect
                    rnorm(n = nrow(sim), mean = 0, sd = 1) # Residuals
  
  # plot(sim$tmax_sum_dif, y = .25  * sim$tmax_sum_dif + -.25 * sim$tmax_sum_dif^2 +
  #        -.25 * sim$tmax_sum +  # Linear effect
  #        -.15 *sim$tmax_sum * sim$tmax_sum_dif)
  # 
  
  pairs.panels(sim[, -c(1,2)])
  
  
  gam1 <- bam(rgr ~ s(tmax_sum_dif) + 
                    s(tmax_sum),
              
              data = sim[sim$site %in% c(3,5),])
  
  gam1 <- bam(rgr ~ ti(tmax_sum_dif) + 
                ti(tmax_sum) +
                te(tmax_sum_dif, tmax_sum),
              
              data = sim[sim$site %in% c(3,8),])
  
  summary(gam1)
  gam.check(gam1)
  visreg(gam1)
  visreg(gam1, xvar = "tmax_sum", by = "tmax_sum_dif", overlay = FALSE)
  
  
  plot(sim[sim$site %in% c(3, 5),]$tmax_sum_dif, sim[sim$site %in% c(3, 5),]$tmax_sum)
  
    
  
 
  

# Testing out rrBLUP ------------------------------------------------------

# library(rrBLUP) 
#   
#   pairs.panels(dat_all_scaled[ , c("rgr", "height_2014", "height_2015", "height_2017", "alive_2017")])
#   
#  
#   ## Calculate slopes of intercepts for each accession
#   fit_height_blup <- lmer(height_2017~ 
#                           #  tmax_sum_dif + 
#                             height_2014 + 
#                            #   (1 | section_block) +
#                             #  (1 | accession) + 
#                               (0 + tmax_sum_dif | accession) +
#                             (1 | section_row) + 
#                             (1 | section_line),
#                           REML = FALSE,
#                           data = dat_all_scaled)
#   
#   summary(fit_height_blup)
#   coef(fit_height_blup)
#   
#   int_coef <- coef(fit_height_blup)$accession
#   int_coef$accession <- as.numeric(rownames(int_coef))
#   int_coef$tmax_sum_dif_slope <- int_coef$tmax_sum_dif
#   int_coef$tmax_sum_dif <- NULL
#   int_coef$height_2014 <- NULL
#   
#   int_coef$tmax_sum_dif
#   summary(int_coef$tmax_sum_dif)
#   
#   plot(int_coef$tmax_sum_dif_slope, int_coef$`(Intercept)`)
#   
#   
#   dat_all_scaled %>%
#     filter(accession < 10) %>%
#     ggplot(aes(x = tmax_sum_dif, y = rgr, group = accession)) +
#     facet_wrap(~accession) +
#     geom_point() + geom_smooth() + theme_bw()
#   
#   # Predictions
#   newdat <- expand.grid(tmax_sum_dif = c(0, 2),
#                         accession = unique(dat_all_scaled$accession),
#                         section_row = "Block1-A",
#                         section_line = "Block1-1",
#                         height_2014 = 0)
#   newdat$pred <- predict(fit_height_blup,
#                       newdata = newdat, re.form = NULL)
#   newdat %>%
#     filter(accession < 10) %>%
#     ggplot(aes(x = tmax_sum_dif, y = pred, col = factor(accession))) +
#     facet_wrap(~accession) +
#     geom_point() + geom_line() + theme_bw() + theme(legend.position = "none")
#   
#   
#   ### Getting genetic data
#   snp_col_names <- colnames(gen_dat[, -c(1, 2)])
#   
#   dat_all_gen <- inner_join(dat_all_scaled, gen_dat)
#   dim(dat_all_gen)
#   dat_all_gen <- as.data.frame(dat_all_gen)
#   
#   dat_all_gen <- left_join(dat_all_gen, int_coef, by = "accession")
#   
#   
#   
#   ## Calculate additive relationship matrix
#   
#   ## Filter down to unique accessions
#   blup_gen <- dat_all_gen %>%
#     dplyr::distinct(accession, .keep_all = TRUE) %>%
#     dplyr::filter(accession %in% dat_all_gen$accession) %>%
#     dplyr::select(accession, snp_col_names)
#   dim(blup_gen)
#   
#   rownames(blup_gen) <- blup_gen$accession
#   blup_gen <- blup_gen[, -1]
#   
#   blup_gen[blup_gen == 0] <- -1
#   blup_gen[blup_gen == 1] <- 0
#   blup_gen[blup_gen == 2] <- 1
#   
#   summary(blup_gen[, 1:10])
#   
#   A <- A.mat(X = blup_gen)
#   
#   
#   ## data frame with columns for phenotype, genotype identify, and any environmental variables
#   dat <- data.frame(rgr = dat_all_gen$rgr,
#                   #  rgr_resids = dat_all_gen$rgr_resids,
#                     height_2017 = dat_all_gen$height_2017,
#                     accession = dat_all_gen$accession,
#                     site = dat_all_gen$site,
#                     section = dat_all_gen$section,
#                     section_block = dat_all_gen$section_block,
#                     height_2015 = dat_all_gen$height_2015,
#                     height_2014 = dat_all_gen$height_2014,
#                     tmax_sum_dif = dat_all_gen$tmax_sum_dif,
#                     tmax_sum = dat_all_gen$tmax_sum,
#                     tmax_sum_dif_slope = dat_all_gen$tmax_sum_dif_slope,
#                     tmin_winter = dat_all_gen$tmin_winter,
#                     cwd = dat_all_gen$cwd,
#                     latitude = dat_all_gen$latitude,
#                     bioclim_18 = dat_all_gen$bioclim_18,
#                     random = rnorm(nrow(dat_all_gen)),
#                   elevation = dat_all_gen$elevation)
#   
#   ## Subset to just seedlings planted in warmer climates
#     tmax_thresh <- ((0 - scaled_var_means[["tmax_sum_dif"]]) / scaled_var_sds[["tmax_sum_dif"]])
#     tmax_thresh
#   #  dat <- dat[dat$tmax_sum_dif > tmax_thresh , ]
#     dim(dat)
#   
#   ## Subset to just down in Chico
#    # dat <- dat[dat$section == section, ]
#     dim(dat)
#   
# 
#   ## Average over accession
#   dat <- dat %>%
#     group_by(accession) %>%
#     summarise_all("mean")
#   dim(dat)
# 
#   
#   ## Set phenotype 
#   pheno = "rgr"
#   
#   # Do 10 fold cross validation
#   folds <- createFolds(y = 1:nrow(dat), k = 5)
#   cors <- NA
#   x=1
#     for(fold in folds){
#       dat_train <- dat
#       dat_train[fold, pheno] <- NA
#   
#   
#      test = kin.blup(data = as.data.frame(dat_train),
#                       geno = "accession",
#                       pheno = pheno,
#                   #   fixed = c("section"),
#                       covariate = c("height_2014"),
#                       K = A)  
#   
#     test
#     
#     test$Vg
#     
#     # Broad sense heritability
#     cat("BSH == ", test$Vg / (test$Vg + test$Ve), "\n") ## Broad sense heritability
#     
#     ## Evaluate accuracy on test set
#     cv_cor <- cor(as.data.frame(dat)[fold, pheno],
#         test$g[ match(as.character(dat$accession[fold]), names(test$g)  )  ])
#     plot(y = as.data.frame(dat)[fold, pheno],
#         x =  test$g[  match(as.character(dat$accession[fold]), names(test$g)  )  ],
#         ylab = pheno, xlab = "breeding value" )
#     
#     cat("CV correlation == ", cv_cor, "... \n")
#     
#     cors[x] <- cv_cor
#     x = x + 1
# 
#     }
#   
#   summary(cors)
#   
#   
#   test2 <- data.frame(g = test$g,   
#                       accession = as.numeric(names(test$g)))
#   dim(test2)
#   test2 <- left_join(test2, dat)
#   dim(test2)
#   plot(test2$tmax_sum, test2$g)
#   cor.test(test2$tmax_sum, test2$g)
#   
#   test2 %>%
#   ggplot(aes(x = pull(test2, pheno), y  = g))  +
#     geom_point() + geom_smooth() + theme_bw()
#   
#   
#   test2 %>%
#     filter(accession %in% dat$accession) %>%
#     # filter(tmax_sum > tmax_thresh) %>% 
#     ggplot(aes(x = g, y  = height_2017))  +
#     geom_point() + geom_smooth() + theme_bw()
#   
#   
#   ## Matching back to garden data
#   test2 <- data.frame(g = test$g,   
#                       accession = as.numeric(names(test$g)))
#   
#   dat_all_gen2 <- left_join(dat_all_gen, test2)
#   
#   
#   dat_all_gen2 %>%
#     select(accession, section, g, rgr, height_2017, tmax_sum, tmax_sum_dif, bioclim_18) %>%
#     group_by(accession, section) %>%
#     summarise_all("mean") %>%
#   ggplot(aes(x = g, y = rgr, color = tmax_sum_dif)) + 
#     geom_point() + 
#     facet_wrap(~section) +
#     theme_bw() + geom_smooth()
#   
#  t = lmer(rgr~ g + height_2014 + (1|section_block), dat_all_gen2)
#  summary(t)
#  
#  ## Check predictive power of breeding value on growth in common gardens
#  fit_g <- gamm4(rgr ~ s(height_2014)
#                        + s(g),
#                        random = ~ (1 | section) +
#                          (1 | accession) +
#                          (1 | section_row) + 
#                          (1 | section_line),
#                        data = dat_all_gen2)
#  check_model(fit_g)
#  plot_gam(fit_g$gam, var = "g", var_index = 2, 
#           xlab = "g", plot_raw_pts = TRUE)
#  
#  
#  newdat <- expand.grid(height_2014 = 0,
#                        section_block = "IFG-1",
#                        g = seq(-3, 3, by = .1),
#                        tmax_sum= c(-2, 0, 2))
#  
#  newdat$predictions <- predict(t, newdat = newdat, re.form = NA)
#  
#  ggplot(newdat, aes(x = g, y = predictions, colour = factor(tmax_sum))) +
#    geom_line() + theme_bw()
#  
#  
#  pairs.panels(dat_all_gen2[ , c("rgr", "g", climate_vars)])
#    
#  

## BGLR section
    
    # #1# Estimated Marker Effects & posterior SDs - for marker effects
    # bHat<- bglr_test$ETA[[1]]$b
    # SD.bHat<- bglr_test$ETA[[1]]$SD.b
    # plot(bHat^2, ylab='Estimated Squared-Marker Effect',
    #      type='o',cex=.5,col=4,main='Marker Effects')
    # #2# Predictions
    # # Total prediction
    # yHat<-bglr_test$yHat
    # tmp<-range(c(y,yHat))
    # plot(yHat~y,xlab='Observed',ylab='Predicted',col=2,
    #      xlim=tmp,ylim=tmp); abline(a=0,b=1,col=4,lwd=2)
    # # Just the genomic part
    # gHat<-as.matrix(X)%*%bglr_test$ETA[[1]]$b
    # plot(gHat~y,xlab='Phenotype',
    #      ylab='Predicted Genomic Value',col=2); abline(a=0,b=1,col=4,lwd=2)
    # #3# Godness of fit and related statistics
    # bglr_test$fit
    # bglr_test$varE # compare to var(y)
    # #4# Trace plots
    # list.files()
    # # Residual variance
    # varE<-scan('./BGLR/test_varE.dat')
    # plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
    # abline(h=bglr_test$varE,col=4,lwd=2);
    # abline(v=bglr_test$burnIn/bglr_test$thin,col=4)
    # # lambda (regularization parameter of the Bayesian Lasso)
    # lambda<-scan('./BGLR/test_ETA_1_lambda.dat')
    # plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
    # abline(h=bglr_test$ETA[[1]]$lambda,col=4,lwd=2);
    # abline(v=bglr_test$burnIn/bglr_test$thin,col=4)
    # 
 
    
    
### Testing out environment interactions
    
   #  
   #  # Do CV with BGLR gxe models
   #  # Figure out how to interpret gxe interaction terms 
   #  # Figure out how to link to climate - simple correlations with breeding value and climate variables?
   #  
   #  
   #  ## Set up data so that rows are accessions and columns are environments
   #  dat_all_gen_int <- dat_all_gen
   #  
   #  Y <- dat_all_gen_int %>%
   #                    dplyr::select(accession, site, rgr, climate_vars)  %>% 
   #                    group_by(accession, site) %>%
   #                    summarise_all("mean") %>%
   #                    tidyr::spread(site, rgr)
   #  
   #  Y <- as.data.frame(Y)
   #  
   #  # Make DF of climate differences in chico
   #  chico_dif <- dat_all_gen_int %>%
   #              dplyr::select(accession, site, rgr, climate_vars_dif, height_2014) %>%
   #              dplyr::filter(site == "Chico") %>%
   #              group_by(accession, site) %>%
   #              summarise_all("mean") %>%
   #              left_join(Y, ., by = "accession")
   #  
   #  ifg_dif <- dat_all_gen_int %>%
   #              dplyr::select(accession, site, rgr, climate_vars_dif, height_2014) %>%
   #              dplyr::filter(site == "IFG") %>%
   #              group_by(accession, site) %>%
   #              summarise_all("mean") %>%
   #              left_join(Y, ., by = "accession")
   #  
   #  
   #  # Genotype matrix - can't have missing data
   #  ## Filter down to unique accessions
   #    bglr_gen <- dat_all_gen %>%
   #       dplyr::distinct(accession, .keep_all = TRUE) %>%
   #       dplyr::filter(accession %in% dat_all_scaled$accession) %>%
   #      dplyr::select(accession, snp_col_names)
   #    dim(bglr_gen)
   #    
   #    ## Make sure in the same order as Y
   #    bglr_gen <- bglr_gen[match(Y$accession, bglr_gen$accession), ]
   # 
   #    sum(bglr_gen$accession != Y$accession)
   #    
   #    dim(bglr_gen)
   #    
   #    rownames(bglr_gen) <- bglr_gen$accession
   #    bglr_gen <- bglr_gen[, -1]
   #    
   #    bglr_gen[bglr_gen == 0] <- -1
   #    bglr_gen[bglr_gen == 1] <- 0
   #    bglr_gen[bglr_gen == 2] <- 1
   #    
   #    summary(bglr_gen[, 1:10])
   #    
   #    A <- A.mat(X = bglr_gen)
   # 
   #  ## Choose environmental variables  
   #  env <- c("Chico", "IFG") # choose any set of environments from 1:ncol(Y)
   #  
   #  nEnv <- length(env)
   #  
   #  y <- as.vector(as.matrix(Y[,env]))
   #  
   #  ## Adjusting rgr for differences in initial sizes
   #  y_adj <- residuals(lm(y ~  c(chico_dif$height_2014, ifg_dif$height_2014),
   #                        na.action = "na.exclude"))
   #  plot(y, y_adj, pch = 19)
   #  
   #  
   #  # Fixed effect (env-intercepts)
   #    envID <- rep(env,each=nrow(Y))
   # 
   #  ## Start model descriptions  
   #  ETA <- list(list(~factor(envID)-1, model="FIXED"))
   # #   ETA <- list(list(X = tmax_dif, model="FIXED"))
   #  
   #  # Main effects of markers
   #    G0 <- kronecker(matrix(nrow=nEnv,
   #                           ncol=nEnv,
   #                           1), A)
   #    
   #    # is equivalent to fitting a genomic regression model using the average performance of each line across environments as a phenotype.
   #    ETA[[2]] <- list(K=G0, model='RKHS')
   #  
   #  # Adding interaction terms
   #    for(i in 1:nEnv){
   #      tmp <- rep(0,nEnv) ; tmp[i] <- 1
   #      G1 <- kronecker(diag(tmp),A)
   #      ETA[[(i+2)]] <- list(K=G1, model='RKHS')
   #    }
   #  
   #  ## Adding in fixed effect
   #    
   #    tmax_dif <- c(chico_dif$tmax_sum_dif, ifg_dif$tmax_sum_dif)
   #     ## Mean impute NAs
   #     tmax_dif[is.na(tmax_dif)] <- 0
   #    
   #    ETA[[length(ETA) + 1]] <- list(X = as.data.frame(tmax_dif), model = 'FIXED')
   #    
   #    
   #  # Start cross validation
   #  prefix = "./BGLR/inter_"  
   #  
   #  ## There are two basic cross-validation schemes used in genome-enabled prediction: (1) predicting the performance of certain proportion of lines that have not been evaluated in any of the observed environments (CV1), and (2) predicting the performance of a proportion of lines that have been evaluated in some environments, but not in others, also called sparse testing (CV2). Another prediction problem that does not involve random cross-validation is predicting one environment using another environment (pairwise environment). The fourth prediction problem consists of predicting one environment (i.e., site-year combination) that was not included in the usual set of testing environments in the evaluation system (leave-one-environment-out); the only available information on this untested environment could be certain char- acteristics that would have been previously collected such as soil type, altitude, longitude, maximum and mini- mum temperature, precipitation during other cropping cycles, etc.
   #    
   #  folds <- createFolds(y = 1:length(y_adj), k = 10)
   #  x = 1
   #  cv_cors <- NA
   #  
   #  for(fold in folds){ 
   #      yNA <- y_adj
   #      yNA[fold] <- NA
   #    
   #    # Model Fitting
   #      # Paper did 55,000 iters and 5000 burnin
   #    fm <- BGLR(y=yNA, ETA=ETA,
   #               nIter=5000,burnIn=2000,saveAt=prefix) 
   #    
   #    
   #    #5# Assesment of correlation in TRN and TST data sets
   #    cor(x = fm$yHat[-fold], y = y_adj[-fold], use = "complete.obs") #TST
   #    cv_cors[x] <- cor(fm$yHat[fold], y_adj[fold], use = "complete.obs") #TRN
   #    
   #    cat("Correlation is:", cv_cors[x], "... \n")
   #    
   #    x = x+1
   #  
   #  }
   #  
   #  summary(cv_cors)
   # 
   #  
   #  ## Model fit
   #  fm$fit
   #  
   #  # Extracting estimates of variance parameters
   #    fm$varE # residual variance
   #    fm$ETA[[2]]$varU # genomic variance (main effect)
   #    vGInt <- rep(NA,nEnv)
   #    for(i in 1:nEnv){ # interaction variances
   #      vGInt[i] <- fm$ETA[[(i+2)]]$varU
   #    }
   #    vGInt
   #    
   #   # Model R2 was computed as the ratio of the sum of the main and interaction variance, relative to the total variance (residual + main effect + interaction). Env, environment.
   #    (fm$ETA[[2]]$varU) /  (fm$ETA[[2]]$varU + fm$varE)
   #    (fm$ETA[[2]]$varU + sum(vGInt)) / (fm$ETA[[2]]$varU + sum(vGInt) + fm$varE)
   #    
   #  # Predictions (this is all within training)
   #    tmpEnv <- 2
   #    plot(y_adj[envID==env[tmpEnv]]~fm$yHat[envID==env[tmpEnv]])
   #    
   #    plot(y_adj ~ fm$yHat)
   #    
   #    # Samples
   #    varE <- scan(paste(prefix,'varE.dat',sep=''))
   #    plot(varE,type='o',cex=.5,col=4)
   #    varU0 <- scan(paste(prefix,'ETA_2_varU.dat',sep=''))
   #    plot(varU0,type='o',cex=.5,col=4)
   #    
   #    varU1 <- matrix(nrow=length(varU0),ncol=nEnv,NA)
   #    for(i in 1:nEnv){
   #      varU1[,i] <- scan(paste(prefix,'ETA_',i+2,'_varU.dat',sep=''))
   #    }
   #    
   #    tmpEnv <- 2
   #    plot(varU1[,tmpEnv],type='o',col=4,cex=.5)
   #    
   #    plot(fm$ETA[[2]]$u, y, col = factor(envID), pch = 19)
   #    
   #    ### ## Plot both gardens
   #    pairs.panels(cbind(fm$ETA[[2]]$u, rbind(chico_dif, ifg_dif)))
   #    
   #    plot(fm$ETA[[2]]$u, c(chico_dif$tmax_sum_dif, ifg_dif$tmax_sum_dif))
   #    cor.test(fm$ETA[[2]]$u, c(chico_dif$tmax_sum_dif, ifg_dif$tmax_sum_dif))
   #    
   #    plot(y = c(fm$ETA[[3]]$u[envID == "Chico"],
   #           fm$ETA[[4]]$u[envID == "IFG"]), 
   #         x = c(chico_dif$tmax_sum, ifg_dif$tmax_sum), pch = 19)
   #    cor.test(c(fm$ETA[[3]]$u[envID == "Chico"],
   #               fm$ETA[[4]]$u[envID == "IFG"]), c(chico_dif$tmax_sum, ifg_dif$tmax_sum))
   #    
   #    
   #    test = cbind(fm$ETA[[2]]$u, rbind(chico_dif, ifg_dif))
   #    
   #    test <- test %>%
   #      dplyr::select(-site) %>%
   #      group_by(accession) %>%
   #      summarise_all("mean")
   #    pairs.panels(test)
   #    
   #    
   #    
   #    ## Plot variance in chico and climate vars dif
   #    pairs.panels(cbind(y_adj[envID == "Chico"], 
   #                       fm$ETA[[3]]$u[envID == "Chico"], 
   #                       chico_dif[, climate_vars_dif]))
   #    
   #    
   #    ## Plot variance in chico and tmax_sum_dif
   #    pairs.panels(cbind(y_adj[envID == "IFG"],
   #                       fm$ETA[[4]]$u[envID == "IFG"], 
   #                       ifg_dif[, climate_vars_dif]))
   #    
   #    test <- read_tsv("./BGLR/inter_ETA_5_b.dat")
   #    test
   #    head(test)
   #    dim(test)
   #    hist(test)
   #    summary(test)
   #    
   #    quantile(pull(test, tmax_dif), c(0.025, 0.975))
   #    
      
 
 

# Testing out RDA  --------------------------------------------------------

  # library(vegan)  
  # 
  # rgr_rda <- rda(dat_all_scaled$rgr ~ height_2015 +  tmax_sum_dif + tmin_winter_dif + cwd_dif + bioclim_04_dif + 
  #                  bioclim_15_dif + bioclim_18_dif + bioclim_19_dif, data = as.data.frame(dat_all_scaled))
  # 
  # rgr_rda
  # 
  # summary(rgr_rda)
  # 
  # ## Look at variance inflation factor
  #  vif.cca(rgr_rda) 
  #  
  # ## Significance of overall model 
  # anova.cca(rgr_rda, parallel=getOption("mc.cores")) 
  # 
  # ## Significance by axis
  # anova.cca(rgr_rda, by="axis", parallel=getOption("mc.cores"))
  # 
  #       
  # ## Plotting RDA
  # plot(rgr_rda, scaling = 3)
  # 
  # 
  # 
  # ## With SNP data
  # dat_all_gen <- left_join(dat_all_scaled, gen_dat)
  # 
  # rgr_rda <- rda(dat_all_gen$rgr ~ ., data = as.data.frame(dplyr::select(dat_all_gen, `1_76174`:`8_59007667`)))
  # 
  # rgr_rda
  # 
  # summary(rgr_rda)
  # 
  # ## Look at variance inflation factor
  # vif.cca(rgr_rda) 
  # 
  # ## Significance of overall model 
  # anova.cca(rgr_rda, parallel=getOption("mc.cores")) 
  # 
  # ## Significance by axis
  # anova.cca(rgr_rda, by="axis", parallel=getOption("mc.cores"))
  # 
  # 
  # ## Plotting RDA
  # plot(rgr_rda, scaling = 3)
  # 
  
  
  


# Models with PC of climate variables -------------------------------------

    # # Null model with only random effects
    # fit_null<- gamm4(rgr  ~ 1,
    #                  random = ~(1|block) + (1|locality) + (1|accession) + (1|row) + (1|line),
    #                  data = dat_all_scaled)
    # 
    # # Height in 2016 only
    # fit_height_only <- gamm4(rgr  ~ s(height_2015),
    #                          random = ~ (1|block) + (1|locality) + (1|accession) + (1|row) + (1|line),
    #                          data = dat_all_scaled)
    # 
    # ## PC Climate variables
    # fit_clim_pca <- gamm4(rgr ~ s(height_2015)
    #                  + PC1_clim_dif
    #                  + s(PC2_clim_dif)
    #                  + s(PC3_clim_dif),
    #                  random = ~(1|block) + (1|locality) + (1|accession) + (1|row) + (1|line),
    #                  data = dat_all_scaled)
    # 
    # fit_clim_pca_clust <- gamm4(rgr ~ s(height_2015)
    #                       + s(PC2_clim_dif, by = PC1_gen_avg),
    #                      # + s(PC2_clim_dif, by = PC1_gen_avg)
    #                       #+ s(PC3_clim_dif, by = PC1_gen_avg)
    #                       #+ s(PC4_clim_dif, by = PC1_gen_avg),
    #                       random = ~(1|block) + (1|locality) + (1|accession) + (1|row) + (1|line),
    #                       data = dat_all_scaled)
    # 
    # 
    # 
    # ## Check PC climate variables
    # check_model(fit_clim_pca_clust)
    # 
    # ## These don't work right if not a smoothing factor!!
    # 
    # plot_gam(fit_clim_pca$gam, var = "PC1_clim_dif", var_index = 2, 
    #          xlab = "PC1_clim_dif", plot_raw_pts = F, PCA = TRUE)
    # 
    # plot_gam(fit_clim_pca$gam, var = "PC2_clim_dif", var_index = 3, 
    #          xlab = "PC2_clim_dif", plot_raw_pts = F, PCA = TRUE)
    # 
    # plot_gam(fit_clim_pca$gam, var = "PC3_clim_dif", var_index = 4, 
    #          xlab = "PC3_clim_dif", plot_raw_pts = F, PCA = TRUE)
    # 
    # plot_gam(fit_clim_pca$gam, var = "PC4_clim_dif", var_index = 5, 
    #          xlab = "PC4_clim_dif", plot_raw_pts = F, PCA = TRUE)
    # 

    
      

# Graveyard ---------------------------------------------------------------



# Spatial autocorrelation -------------------------------------------------



## Spatial autocorrelation in residuals
# 
#     resids <- data.frame(y = dat_all_pca$row[!is.na(dat_all_pca$rgr)], 
#                          x = dat_all_pca$line[!is.na(dat_all_pca$rgr)], 
#                          resids = residuals(fit))
#     coordinates(resids) <- c("x", "y")
#     
#     bubble(resids, zcol = "resids")
#     
#     
#     ## Variogram
#     library(gstat)
#     
#     vario <- variogram(resids ~ 1, data = resids)
#     plot(vario)
#     
#     library(ncf)
#     
#     plot(spline.correlog(x = resids$x, y = resids$y, z = resids$resids, resamp = 0))
#     
#     
#     ## Moran's I
#     dists <- as.matrix(dist(cbind(resids$x, resids$y)))
#     
#     dists.inv <- 1/dists
#     diag(dists.inv) <- 0
#     
#     library(ape)
#     
#     Moran.I(resids$resids, dists.inv)



# GLMM --------------------------------------------------------------------

  # ### PC model
  # fit_nocluster = lmer(rgr  ~ PC1 + I(PC1^2) +
  #                        + PC2 + I(PC2^2)  
  #                      + PC3 + I(PC3^2) 
  #                      + height_2015_scaled
  #                      + (1 | block) # + (1 | site) 
  #                      + (1 | locality) + (1 | accession), 
  #                      data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # ### PC model
  # fit = lmer(rgr  ~ PC1*cluster_assigned + I(PC1^2)*cluster_assigned 
  #            + PC2*cluster_assigned + I(PC2^2)*cluster_assigned  
  #            + PC3*cluster_assigned + I(PC3^2) *cluster_assigned
  #            + cluster_assigned
  #            + height_2015_scaled
  #            + (1 | block)
  #            + (1 | locality) + (1 | accession), 
  #            data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # 
  # fit_pca = lmer(rgr  ~ PC1*PC2_avg + I(PC1^2)*PC2_avg 
  #                + PC2*PC2_avg + I(PC2^2)*PC2_avg  
  #                + PC3*PC2_avg + I(PC3^2) *PC2_avg
  #                + height_2015_scaled
  #                + (1 | block)
  #                + (1 | locality) + (1 | accession), 
  #                data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # 
  # anova(fit, fit_nocluster, fit_pca)
  # 
  # ranef(fit)
  # 
  # ## Model summary
  # 
  # summary(fit)
  # print(r.squaredGLMM(fit))
  # sjp.lmer(fit, type = "fe", p.kr = FALSE)
  # sjp.lmer(fit, type = "eff")
  # 
  # sjp.int(fit)



# Random forest models ----------------------------------------------------

  library(randomForest)
  library(forestFloor)

  ran_for_data <- dat_all_scaled %>%
   # dplyr::filter(site == "chico") %>%
    dplyr::filter(!is.na(rgr)) %>%
    dplyr::select(rgr, tmax_sum_dif, height_2015, height_2017, site, section_block, accession, cwd_dif, latitude,
                  section_row, section_line)

  height_rf = randomForest(height_2017 ~ tmax_sum_dif + 
                             latitude + 
                             height_2015 +
                             as.numeric(factor(ran_for_data$site)) +
                             as.numeric(factor(ran_for_data$section_block)) +
                             as.numeric(factor(ran_for_data$section_row)) +
                             as.numeric(factor(ran_for_data$section_line)) +
                             as.numeric(factor(ran_for_data$accession)),
                           data = ran_for_data,
                           keep.inbag= T,
                           importance = TRUE)

  height_rf

  varImpPlot(height_rf)


  ff = forestFloor(
    rf.fit = height_rf, # mandatory
    X = ran_for_data, # mandatory
    calc_np = FALSE, # TRUE or FALSE both works, makes no difference
    binary_reg = FALSE # takes no effect here when rfo$type="regression"
  )

  #plot partial functions of most important variables first
  plot(ff, # forestFloor object
       # plot_seq = 1:7, # optional sequence of features to plot
       orderByImportance=TRUE # if TRUE index sequence by importance, else by X line
  )
  
  plot(ff, plot_seq = 3, ylim = c(-50, 50))



  library(ranger)

  height_ranger <- ranger(rgr ~ ., data = ran_for_data, importance = "permutation")

  sort(height_ranger$variable.importance, decreasing = T)

  importance_pvalues(height_ranger)

  height_ranger$r.squared




# Gradient boosted machines -----------------------------------------------

  # library(gbm)
  # library(dismo)
  
  
  # gbm_data <- dat_all_pca %>%
  #   dplyr::filter(!is.na(rgr)) %>%
  #   #  dplyr::filter(!is.na(cluster_assigned)) %>%
  #   dplyr::select(rgr, 
  #                 #cluster1_avg:cluster3_avg,
  #                 height_2015,
  #                 PC1_avg:PC20_avg, 
  #                 #PC1:PC5, 
  #                 block, site,
  #                 #   locality, accession,
  #                 climate_vars_dif,
  #   )
  # names(gbm_data)[-1]
  # 
  # 
  # # height_gbm = gbm(rgr ~ ., data = gbm_data, n.trees = 1000, shrinkage = 0.01,
  # #                  interaction.depth = 5, bag.fraction = 0.5, train.fraction = 1.0,
  # #                  cv.folds = 5, keep.data = TRUE, distribution = "gaussian")
  # 
  # ## Step function
  # 
  # response <-  which(colnames(gbm_data) == "rgr")
  # predictors <- which(colnames(gbm_data) %in% setdiff(names(gbm_data), c("rgr", "name")))
  # 
  # ## Step function
  # height_gbm <- gbm.step(data = as.data.frame(gbm_data),
  #                        gbm.y = response, gbm.x = predictors,
  #                        tree.complexity = 5,
  #                        learning.rate = 0.001, bag.fraction = 0.5,
  #                        family = "gaussian")
  # 
  # summary(height_gbm)
  # 
  # ## Simplify by dropping variables
  # 
  # height_gbm_simp <- gbm.simplify(height_gbm, n.drops = 5)
  # 
  # 
  # gbm.plot(height_gbm, n.plots = 12)
  # 
  # gbm.plot.fits(height_gbm, v = 1:6)
  # 
  # find.int <- gbm.interactions(height_gbm)
  # 
  # find.int$interactions
  # 
  # find.int$rank.list
  # 
  # plot(gbm_data$rgr ~ gbm_data$accession)
  # 
  # 
  # 
  # print(height_gbm)
  # 
  # summary(height_gbm)
  # 
  # gbm.perf(height_gbm, method = "OOB")
  # gbm.perf(height_gbm, method = "test")
  # 
  # best.iter <- gbm.perf(height_gbm, method="cv")
  # print(best.iter)
  # 
  # summary(height_gbm, best.iter)
  # 
  # plot(height_gbm, 15, best.iter, las = 1, continuous.resolution = 1000)
  # 
  # interact.gbm(height_gbm, gbm_data, i.var = c(5, 17), n.trees = best.iter)
  # 
  # test = partial(height_gbm, "height_2015_scaled", n.trees = best.iter, grid.resolution = 100,
  #                smooth = TRUE)
  # 
  # plotPartial(test, smooth = TRUE)
  # 
  # 
  # test = partial(height_gbm, c("PC5_avg", "height_2015_scaled"), n.trees = best.iter, grid.resolution = 100)
  # plotPartial(test)
  


  # 
  # # library(xgboost)
  # # library(pdp)
  # 
  # previous_na_action <- options('na.action')
  # options(na.action='na.pass')
  # # Do your stuff...
  # 
  # boost_dat = sparse.model.matrix(rgr ~.-1, data = gbm_data, na.action = "na.pass")
  # 
  # dim(boost_dat)
  # 
  # options(na.action=previous_na_action$na.action)
  # 
  # xg <- xgb.cv(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear", 
  #              nfold = 3, nrounds = 2)
  # 
  # xg <- xgboost(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear",
  #               nrounds = 50, max_depth = 4, fill = TRUE)
  # 
  # xg
  # 
  # xg_import <- xgb.importance(feature_names = colnames(boost_dat), model = xg, )
  # 
  # xg_import
  # 
  # xgb.plot.importance(importance_matrix = xg_import)
  # 
  # xgb.plot.tree(model = xg)
  # 
  # explain <- buildExplainer(xg, trainingData = boost_dat, type = "regression")
  # 
  # partial(xg)
  # 
  # ### 
  # # Tune an XGBoost model using 10-fold cross-validation
  # #library(caret) # functions related to classification and regression training
  # 
  # boston.xgb <- train(x = data.matrix(boost_dat),
  #                     y = gbm_data$rgr, method = "xgbTree", metric = "Rsquared",
  #                     trControl = trainControl(method = "cv", number = 10),
  #                     tuneLength = 10)



# Using h2o for RF and GBM ------------------------------------------------

  # 
  # #library(h2o)
  # 
  # 
  # ## Tutorial: https://github.com/h2oai/h2o-3/blob/master/h2o-docs/src/product/tutorials/gbm/gbmTuning.Rmd
  # 
  # ## Initialize h2o cluster
  # h2o.shutdown()
  # h2o.init(nthreads = -1, enable_assertions = FALSE) 
  # 
  # ## Write data to csv file first
  # gbm_data <- dat_all_pca %>%
  #   dplyr::filter(!is.na(rgr)) %>%
  #   dplyr::select(rgr, cluster1_avg:cluster_assigned, height_2015_scaled,
  #                 PC1_avg:PC5_avg, 
  #                 #  lat:RH, Tmax_seas:Tmin,
  #                 PC1:PC5, 
  #                 locality, accession, block, 
  #                 climate_vars_dif)
  # 
  # ## Change to h2o format
  # gbm_data_h2o <- as.h2o(x = gbm_data, key = gbm_data_h2o)
  # 
  # summary(gbm_data_h2o) ## Make sure factors are read correctly
  # 
  # 
  # 
  # ### SPLIT FIT PREDICT 
  #
  # 
  # response <-  "rgr"
  # predictors <- setdiff(names(gbm_data_h2o), c(response, "name"))
  # 
  # ## Split data for training, testing, and validation
  # ## May not be correct since data is not independent
  # 
  # splits <- h2o.splitFrame(
  #   data = gbm_data_h2o, 
  #   ratios = c(0.7,0.15),   ## only need to specify 2 fractions, the 3rd is implied
  #   destination_frames = c("train.hex", "valid.hex", "test.hex"), seed = 1234
  # )
  # train <- splits[[1]]
  # valid <- splits[[2]]
  # test  <- splits[[3]]
  # 
  # 
  # print(paste("Training data has", ncol(train), "lines and", nrow(train), "rows, valid has",
  #             nrow(valid), "rows, test has", nrow(test)))
  # 
  # ### Run gbm model
  # gbm <- h2o.gbm(x = predictors,
  #                y = response,
  #                training_frame = h2o.rbind(train, valid),
  #                # training_frame =  train,
  #                # validation_frame = valid,
  #                nfolds = 5,
  #                
  #                model_id = "gbm",
  #                
  #                ## more trees is better if the learning rate is small enough 
  #                ## here, use "more than enough" trees - we have early stopping
  #                ntrees = 10000,   
  #                
  #                max_depth = 10,
  #                
  #                ## smaller learning rate is better (this is a good value for most datasets, but see below for annealing)
  #                learn_rate=0.01,                                                         
  #                
  #                ## early stopping once the validation MSE doesn't improve by 
  #                # at least 0.01% for 5 consecutive scoring events
  #                stopping_rounds = 5, stopping_tolerance = 1e-4, stopping_metric = "MSE", 
  #                
  #                ## sample 80% of rows per tree
  #                sample_rate = 0.8,                                                       
  #                
  #                ## sample 80% of columns per split
  #                col_sample_rate = 0.8,                                                   
  #                
  #                ## fix a random number generator seed for reproducibility
  #                seed = 1234,                                                             
  #                
  #                ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
  #                score_tree_interval = 10 
  #                
  # )
  # 
  # ###### TESTING
  # 
  # ## Depth 10 is usually plenty of depth for most datasets, but you never know
  # hyper_params = list( 
  #   ## restrict the search to the range of max_depth established above
  #   max_depth = seq(1,30,2),                                      
  #   
  #   ## search a large space of row sampling rates per tree
  #   sample_rate = seq(0.2,1,0.01),                                             
  #   
  #   ## search a large space of column sampling rates per split
  #   col_sample_rate = seq(0.2,1,0.01),                                         
  #   
  #   ## search a large space of column sampling rates per tree
  #   col_sample_rate_per_tree = seq(0.2,1,0.01),                                
  #   
  #   ## search a large space of how column sampling per split should change as a function of the depth of the split
  #   col_sample_rate_change_per_level = seq(0.9,1.1,0.01),                      
  #   
  #   ## search a large space of the number of min rows in a terminal node
  #   min_rows = 2^seq(0,log2(nrow(train))-1,1),                                 
  #   
  #   ## search a large space of the number of bins for split-finding for continuous and integer columns
  #   nbins = 2^seq(4,10,1),                                                     
  #   
  #   ## search a large space of the number of bins for split-finding for categorical columns
  #   nbins_cats = 2^seq(4,12,1),                                                
  #   
  #   ## search a few minimum required relative error improvement thresholds for a split to happen
  #   min_split_improvement = c(0,1e-8,1e-6,1e-4),                               
  #   
  #   ## try all histogram types (QuantilesGlobal and RoundRobin are good for numeric columns with outliers)
  #   histogram_type = c("UniformAdaptive","QuantilesGlobal","RoundRobin")       
  # )
  # 
  # search_criteria = list(
  #   ## Random grid search
  #   strategy = "RandomDiscrete",      
  #   
  #   ## limit the runtime to 60 minutes
  #   max_runtime_secs = 3600,         
  #   
  #   ## build no more than 100 models
  #   max_models = 100,                  
  #   
  #   ## random number generator seed to make sampling of parameter combinations reproducible
  #   seed = 1234,                        
  #   
  #   ## early stopping once the leaderboard of the top 5 models is converged to 0.1% relative difference
  #   stopping_rounds = 5,                
  #   stopping_metric = "MSE",
  #   stopping_tolerance = 1e-3
  # )
  # 
  # grid <- h2o.grid(
  #   ## hyper parameters
  #   hyper_params = hyper_params,
  #   
  #   ## full Cartesian hyper-parameter search
  #   search_criteria = search_criteria,
  #   
  #   ## which algorithm to run
  #   algorithm="gbm",
  #   
  #   ## identifier for the grid, to later retrieve it
  #   grid_id="final_grid",
  #   
  #   ## standard model parameters
  #   x = predictors, 
  #   y = response, 
  #   training_frame = train, 
  #   validation_frame = valid,
  #   
  #   ## more trees is better if the learning rate is small enough 
  #   ## here, use "more than enough" trees - we have early stopping
  #   ntrees = 10000,                                                            
  #   
  #   ## smaller learning rate is better
  #   ## since we have learning_rate_annealing, we can afford to start with a bigger learning rate
  #   learn_rate = 0.05,                                                         
  #   
  #   ## learning rate annealing: learning_rate shrinks by 1% after every tree 
  #   ## (use 1.00 to disable, but then lower the learning_rate)
  #   learn_rate_annealing = 0.99,                                               
  #   
  #   ## fix a random number generator seed for reproducibility
  #   seed = 1234,                                                             
  #   
  #   ## early stopping once the validation AUC doesn't improve by at least 0.01% for 5 consecutive scoring events
  #   stopping_rounds = 5,
  #   stopping_tolerance = 1e-4,
  #   stopping_metric = "MSE", 
  #   
  #   ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
  #   score_tree_interval = 10                                                
  # )
  # 
  # ## by default, display the grid search results sorted by increasing logloss (since this is a classification task)
  # grid                                                                       
  # 
  # ## sort the grid models by decreasing AUC
  # sortedGrid <- h2o.getGrid("final_grid", sort_by="MSE", decreasing = FALSE)    
  # sortedGrid
  # 
  # gbm <- h2o.getModel(sortedGrid@model_ids[[1]])
  # 
  # 
  # # Run DRF
  # drf <- h2o.randomForest(x = predictors,
  #                         y = response,
  #                         training_frame = h2o.rbind(train, valid),
  #                         # training_frame =  train,
  #                         # validation_frame = valid,
  #                         nfolds = 5,
  #                         model_id = "rf",
  #                         ntrees            = 250,
  #                         max_depth         = 30)
  # 
  # 
  # # 4- Score on holdout set & report
  # train_rmse_gbm  <- h2o.rmse(gbm, train = TRUE)
  # xval_rmse_gbm   <- h2o.rmse(gbm, xval = TRUE)
  # test_perf_gbm <- h2o.performance(model = gbm, newdata = test)
  # test_rmse_gbm   <- h2o.rmse(object = test_perf_gbm)
  # print(paste0("GBM rmse TRAIN = ", train_rmse_gbm, ", rmse XVAL = ", xval_rmse_gbm, ", rmse TEST = ",
  #              test_rmse_gbm))
  # 
  # train_rmse_drf  <- h2o.rmse(drf, train = TRUE)
  # xval_rmse_drf   <- h2o.rmse(drf, xval = TRUE)
  # test_perf_drf <- h2o.performance(model = drf, newdata = test)
  # test_rmse_drf   <- h2o.rmse(object = test_perf_drf)
  # print(paste0("DRF rmse TRAIN = ", train_rmse_drf, ", rmse XVAL = ", xval_rmse_drf, ", rmse TEST = ",
  #              test_rmse_drf))
  # 
  # 
  # ## Show detailed model summary
  # gbm
  # 
  # drf
  # 
  # 
  # ## Get the Mean Squared Error on the validation set
  # h2o.mse(h2o.performance(gbm, newdata = valid)) ## Best is 0.009
  # 
  # h2o.varimp(gbm)
  # h2o.varimp_plot(gbm, ncol(train))
  # 
  # h2o.varimp(drf)
  # h2o.varimp_plot(drf, ncol(train))
  # 
  # h2o.partialPlot(gbm, h2o.rbind(train, valid))
  # 
  # h2o.saveModel(gbm, "h2o_model")
  # 
  # 
  # 
  # 
  # ### Make predictions
  # 
  # pred_gbm <- h2o.predict(object = gbm,
  #                         newdata = test)
  # 
  # plot(as.data.frame(test)$rgr, as.data.frame(pred_gbm)$predict, pch = 19, las = 1,
  #      xlab = "Observed", ylab = "Predicted")
  # abline(a = 0, b = 1, lwd = 3, col = "steelblue")
  # 
  # 
  # ## Partial plots
  # 
  # h2o.partialPlot(gbm, h2o.rbind(train, valid), c("height_2015_scaled", "DD_0_dif"))
  # h2o.partialPlot(gbm, h2o.rbind(train, valid), "prov", nbins = 100)








