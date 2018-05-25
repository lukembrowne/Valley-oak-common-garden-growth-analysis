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
 load("../gam_cluster_2018-05-25.Rdata")


# Initialize variables
  gam_list <- list()
  x = 1

## Convert fixed and random parts of formula to actual formula object
     # No idea why we need to call formula twice, but it works
     make_formula <- function(fixed, random){
       return(formula(formula(enquote(paste0(fixed, random)))))
     }

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
 snp <- paste0("snp_", snp_col_names[snp_index])
 
 cat("Working on: ", snp, "... number: ", x, " ...\n" )

  # Make sure there's at least 3 levels
    if(length(table(dat_snp[, snp])) < 3){
      cat("Skipping because not at least 3 levels... \n")
      next
    }


# Formula for fixed effects
    fixed_effects <- paste0(paste0("height_2017 ~ site + 
                                    s(height_2014) + te(", 
                                   climate_var_dif,") + te(",snp,",k = 3) + ti(",
                                   climate_var_dif, ",", snp,", k = 3) +
                                   s(PC1_gen) + s(PC2_gen) + s(PC3_gen)"))

  
  
   gam_snp = bam(formula = make_formula(fixed_effects, random_effects),
                 data = dat_snp, 
                 discrete = TRUE,
                 nthreads = 8,
                 method = "fREML", family = "tw")

  # Attach data so that we can use visreg later if we need
   # gam_snp$data <- dat_snp

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

            # Rsquared adjusted
            r_sq_adj = unlist(lapply(gam_list_summary, function(x) x$r.sq)),
            
            # Estimated degrees of freedom
            edf = unlist(lapply(gam_list, function(x) sum(x$edf))),
            
            # AIC
            aic = unlist(lapply(gam_list, function(x) AIC(x))),
            
            # P value of interaction term
            ## IF ORDER OF VARIABLES IN MODEL EVER CHANGES - THIS NEEDS TO CHANGE TOO!!
             p_val_int = unlist(lapply(gam_list_summary, function(x) x$s.table[4, "p-value"]))

            # Percent change in height with 3 degree increase
            # height_change_3_deg = unlist(lapply(lapply(gam_list, function(x) predict(x, 
            #                                      newdata = predictions_3deg  %>% 
            #                                        dplyr::distinct(get(climate_var_dif), 
            #                                                        .keep_all = TRUE), 
            #                                      type = "response")),
            #                   function(x) ((x[2] - x[1]) / x[1])*100))
            
            )
 
 head(gam_mod_summary_df)
 
 write.csv(gam_mod_summary_df, file = paste0("./model_summaries/gam_summaries_", task_id, ".csv"),
                               row.names = FALSE)


###### LAGNIAPPE








