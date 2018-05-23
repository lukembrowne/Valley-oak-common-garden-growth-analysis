
# Make good name for gam_list - unique to task
# Rename hello R based on cliamte variable or something?

# Read in command line arguments
	args<-commandArgs(trailingOnly=TRUE)

  print(args)

	task_id <- as.numeric(args[1])
  interval <- as.numeric(args[2])
  climate_var_dif <- as.character(args[3])

# Load libraries
	library(mgcv)
	library(tidyverse)

# Load data
 load("../gam_cluster.Rdata")


# Initialize variables
  gam_list <- list()
  x = 1

## Convert fixed and random parts of formula to actual formula objetc
     # No idea why we need to call formula twice, but it works
     make_formula <- function(fixed, random){
       return(formula(formula(enquote(paste0(fixed, random)))))
     }

# Formula for fixed effects
     fixed_effects <- paste0(paste0("height_2017 ~ site + s(height_2014) + te(", 
                              climate_var_dif,") + te(snp_dbl, k = 3) + ti(",
                              climate_var_dif, ", snp_dbl, k = 3)"))

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
 snp <- snp_col_names[snp_index]
 
 cat("Working on: ", snp, "...\n" )
  
   #dat_snp$snp_factor <- factor(pull(dat_snp, snp))
   dat_snp$snp_dbl <- pull(dat_snp, snp)
  
   gam_snp = bam(formula = make_formula(fixed_effects, random_effects),
                 data = dat_snp, 
                 discrete = TRUE,
                 nthreads = 8,
                 method = "fREML", family = "tw")

  # Attach data so that we can use visreg later if we need
  gam_snp$data <- dat_snp

  # Save into list
   gam_list[[x]] <- gam_snp
   names(gam_list)[x] <- snp

   x = x + 1

   
 # Save individual .Rdata file for each model
    cat("Saving gam model to file... \n")

     gam_name <- paste0("gam_", climate_var_dif, "_", snp)
     assign(gam_name, gam_snp) # Assign model the name

     save(list = gam_name, file = paste0("./gam_mods_data/", gam_name, ".Rdata"))

 # Timing of loop
  end_time <- Sys.time()

   cat("This loop took: ")
   print(end_time - start_time)   
     
}


# Output model summaries


# Run summary for each model
 gam_list_summary <- lapply(gam_list, summary)
 
# Save into dataframe 
 gam_mod_summary_df <- data.frame(
            # climate var
            climate_var = climate_var_dif,
   
            # Snp name
            snp = names(gam_list),
            # Deviance explained
            dev_explained = unlist(lapply(gam_list_summary, function(x) x$dev.expl)),
            
            # Estimated degrees of freedom
            edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
            
            # AIC
            aic = unlist(lapply(gam_list, function(x) x$aic)),
            
            # P value of interaction term
            p_val_int = unlist(lapply(gam_list_summary, function(x) x$s.table[paste0("ti(", climate_var_dif, ",snp_dbl)"), "p-value"]))
            
            )
 
 head(gam_mod_summary_df)
 
 write.csv(gam_mod_summary_df, file = paste0("./model_summaries/gam_summaries_", task_id, ".csv"),
                               row.names = FALSE)


###### LAGNIAPPE

## Creating dat_snp dataframe

 # dat_bv2 <- dat_bv
 # dat_bv2$accession <- as.numeric(as.character(dat_bv$accession))
 # 
 # dat_snp = left_join(dat_bv2, bind_cols(gen_dat_clim[, "accession"],
 #                                        bglr_gen_scaled), 
 #                     by = "accession")
 # 
 # dat_snp$accession <- factor(dat_snp$accession)








