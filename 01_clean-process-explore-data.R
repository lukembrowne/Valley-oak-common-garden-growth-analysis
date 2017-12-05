
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(maptools)
library(scatterpie)


rm(list = ls()) ## Clean workspace

# Load data ---------------------------------------------------------------

  ## Load in garden data 
  dat_all <- read_csv("./data/cleaned_data/Growth data  2016_2017 census 2017_11_15.csv")
  dat_all$lat <- NULL; dat_all$lon <- NULL; dat_all$elev <- NULL ## Clean up before merge
  
  ## Read in climate data of individuals and merge with 
  climate_ind <- read_csv("./data/cleaned_data/Climate data for all individuals 2017_09_27.csv")
  
  ## Merge data
  dat_all <- left_join(dat_all, climate_ind, by = c("prov", "mom", "progeny"))
  
  ## Load in climate of common garden locations - 2010 to 2015
  garden_climate <- read_csv("./data/cleaned_data/CommonGardenClimate 2017_09_28.csv")
  
  ## Average across 2010 - 2015
  garden_climate <- garden_climate %>%
    group_by(site) %>%
    summarise_all(mean)

  ## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))



# Spatial locations in garden ---------------------------------------------

  
  
  ### Converting row and column numbers to XY locations
  
  ## Check for duplicate locations within each section - should all be 0
  sum(table(dat_all$row[dat_all$section == "block1"], 
            dat_all$column[dat_all$section == "block1"]) > 1)
  sum(table(dat_all$row[dat_all$section == "north"], 
            dat_all$column[dat_all$section == "north"]) > 1 )
  sum(table(dat_all$row[dat_all$section == "north_annex"], 
            dat_all$column[dat_all$section == "north_annex"]) > 1)
  sum(table(dat_all$row[dat_all$section == "placerville"], 
            dat_all$column[dat_all$section == "placerville"]) > 1)
  
  ## With proper ordering, reversed so plotting looks OK
  dat_all$row_orig <- dat_all$row
  dat_all$row <- as.numeric(factor(dat_all$row, 
                                   levels = rev(c( "A","B","C","D","E","F","G","H","I","J","K","L","M",
                                                   "N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
                                                   "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ",
                                                   "AK", "AL","AM","AN","AO","AP")))) 

  ## Adding a row of 3 tree sizes in placerville for the road
  dat_all$row[dat_all$row_orig == "AD"]
  
  dat_all$row[dat_all$row <= 13 & dat_all$section == "placerville"] <- dat_all$row[dat_all$row <= 13 
                                                             & dat_all$section == "placerville"] - 3
  
  
  table(dat_all$section)
  

  # Add set values so that xy coordiates don't overlap across sections
  
  dat_all$row[dat_all$section == "north"] <- dat_all$row[dat_all$section == "north"] + 50
  dat_all$row[dat_all$section == "north_annex"] <- dat_all$row[dat_all$section == "north_annex"] + 100
  dat_all$row[dat_all$section == "placerville"] <- dat_all$row[dat_all$section == "placerville"] + 250
  
  plot(dat_all$column, dat_all$row, pch = 22)


  
  

# Formatting height -------------------------------------------------------
    
  ### Setting max height to main height if secondary stem is larger
  for(i in 1:nrow(dat_all)){
    if(!is.na(dat_all$height_max_2017[i])){
      if(dat_all$height_max_2017[i] > dat_all$height_2017[i]){
        dat_all$height_2017[i] <- dat_all$height_max_2017[i]
      }
    }
    
    if(!is.na(dat_all$height_max_2016[i])){
      if(dat_all$height_max_2016[i] > dat_all$height_2016[i]){
        dat_all$height_2016[i] <- dat_all$height_max_2016[i]
      }
    }
  }
  
  
  #### Round height values to the nearest 5
  dat_all$height_2016 <- round(dat_all$height_2016 / 5) * 5
  dat_all$height_2017 <- round(dat_all$height_2017 / 5) * 5
  
  
  ## Calculate relative growth rates
  dat_all$rgr <- (log(dat_all$height_2017) - log(dat_all$height_2016))/1
  
  summary(dat_all$rgr)
  
  ## Plots of 2016 vs 2017
  plot(jitter(dat_all$height_2016), jitter(dat_all$height_2017),
       pch = 19, cex = 0.5, las = 1)
  abline(a = 0, b = 1, lwd = 2, col = "steelblue")
  
  ## Plot of RGR
  plot(dat_all$height_2016, dat_all$rgr, pch = 19, cex = 0.5)
  cor.test(dat_all$height_2016, dat_all$rgr)
  
  ## View 'outliers
    # View(dat_all[dat_all$rgr < -1 & !is.na(dat_all$rgr), ])
     
    # View(dat_all[dat_all$rgr > 1 & !is.na(dat_all$rgr), ])


# Reading in GBS data -----------------------------------------------------

  
  ## Admixture results - Read in Admixture ancestry estimates
  admix_results <- read.table("./data/GBS_data/gbs451.GDP4.AN50.biallelic.QD10.filtered.final.recode.pruned_ld.3.Q")
  
  ## Read in family names output from PLINK/Admixture
  admix_moms <- read_delim("./data/GBS_data/gbs451.GDP4.AN50.biallelic.QD10.filtered.final.recode.pruned_ld.fam", delim = " ", col_names = FALSE)
  
  admix_results <- cbind(admix_moms$X1, admix_results)
  colnames(admix_results) <- c("mom", paste("cluster",
                                            1:(ncol(admix_results)-1), sep = ""))
  admix_results$mom <- as.character(admix_results$mom)
  head(admix_results)
  
  ## Read in datasheet that has provenance names of each individual with GBS data
  ## This datasheet was made by hand
  admix_moms_provs <- read_delim("./data/GBS_data/gbs451 samples provenance list.txt", delim = "\t")
  
  
  ## Join with admix_results
  admix_results <- left_join(admix_results, admix_moms_provs)
  sum(is.na(admix_results$prov)) ## Should equal to 0 - might have changed since NOV prov format gets messed up by excel and changed to date
  # View(admix_results[is.na(admix_results$prov), ])
  
  
  ### Average cluster assignment at the provenance level
  ## So that we are able to use information from samples across the provenance 
  ## instead of just samples that have maternal trees in the common garden
  
  admix_results_avg = admix_results %>% 
    group_by(prov) %>%
    dplyr::select(prov, contains("cluster")) %>%
    summarise_all(mean, na.rm = TRUE)
  
  ## Rename columns
  colnames(admix_results_avg)[grep("cluster", colnames(admix_results_avg))] <- paste(colnames(admix_results_avg)[grep("cluster", colnames(admix_results_avg))], "_avg", sep = "" ) ## Rename columns
  
  
  ## Barplot of ancestry compositions
  barplot(t(as.matrix(admix_results[, c(2,3,4)])), col=rainbow(3),
          xlab="Individual #", ylab="Ancestry", names.arg = admix_results$prov)


  ## Join back with main dataset
  dat_all <- left_join(dat_all, admix_results_avg, by = "prov")
  
  
  ## Assign cluster based on majority assignment
  dat_all$cluster_assigned <- NA
  dat_all$cluster_assigned[dat_all$cluster1_avg > 0.50 & !is.na(dat_all$cluster1_avg)] <- "1"
  dat_all$cluster_assigned[dat_all$cluster2_avg > 0.50 & !is.na(dat_all$cluster2_avg)] <- "2"
  dat_all$cluster_assigned[dat_all$cluster3_avg > 0.50 & !is.na(dat_all$cluster3_avg)] <- "3"
  dat_all$cluster_assigned[dat_all$cluster4_avg > 0.50 & !is.na(dat_all$cluster4_avg)] <- "4"
  dat_all$cluster_assigned[dat_all$cluster5_avg > 0.50 & !is.na(dat_all$cluster5_avg)] <- "5"
  
  
  table(dat_all$cluster_assigned)
  
  dat_all$cluster_assigned <- factor(dat_all$cluster_assigned) # Refactor
  
  
  ## Plot of RGR by cluster
  boxplot(rgr ~ cluster_assigned, dat_all)
  
  
  ## How many seedlings not assigned?
  sum(is.na(dat_all$cluster_assigned)) 
  
  
  
  ### Plot circle graph onto map of population assignment
  cali_fort <- fortify(cali_outline) ## Fortify to print california outline
  
  dat_all_prov_avg <- dat_all %>% ## Average across provenances
    mutate(cluster_assigned = as.numeric(cluster_assigned)) %>%
    group_by(prov) %>%
    summarise_all(mean, na.rm = T) %>%
    mutate(cluster_assigned = factor(cluster_assigned))
  
  ## Make plot
  ggplot() + geom_scatterpie(aes(x = lon, y = lat, group = prov, r = .2),
                             data = dat_all_prov_avg,
                             cols = colnames(dplyr::select(dplyr::select(dat_all_prov_avg, 
                                                                         contains("cluster")),
                                                           contains("_avg"))))+
    coord_equal() + 
    geom_path(data = cali_fort, aes(x = long, y = lat, group = group )) +
    theme_bw()
  
  ggplot(dat_all_prov_avg[dat_all_prov_avg$cluster_assigned != "NaN",]) + 
    geom_point(aes(x = lon, y = lat, col = cluster_assigned), size = 3) +
    coord_equal() + 
    geom_path(data = cali_fort, aes(x = long, y = lat, group = group )) +
    theme_bw()
  
  
  

# Read in PCA of SNPs -----------------------------------------------------

  #### Read in data
  gen_pca <- read_table("./data/GBS_data/plink.eigenval", col_names = FALSE)
  
  gen_pca / sum(gen_pca) ## Percent contribution to variation
  
  
  gen_pca_vals <- read_tsv("./data/GBS_data/plink.eigenvec")
  
  plot(gen_pca_vals$PC1, gen_pca_vals$PC2, pch = 19)
  text(gen_pca_vals$PC1, gen_pca_vals$PC2, gen_pca_vals$FID)
  
  
  ## Join with provenance names and average across provenances - same as what's done with admixture
  pca_results <- left_join(gen_pca_vals, admix_moms_provs, by = c("FID" =  "mom"))
  sum(is.na(pca_results$prov)) ## Should equal to 0 - might have changed since NOV prov format gets messed up by excel and changed to date
  
  ### Average PCA Values for each provenance
  ## So that we are able to use information from samples across the provenance instead of just samples that have maternal trees in the common garden
  
  pca_results_avg = pca_results %>% 
    group_by(prov) %>%
    dplyr::select(prov, PC1:PC20) %>%
    summarise_all(mean, na.rm = TRUE) 
  colnames(pca_results_avg)[grep("PC", colnames(pca_results_avg))] <- paste(colnames(pca_results_avg)[grep("PC", colnames(pca_results_avg))], "_avg", sep = "" ) ## Rename columns
  
  ## Join back with main dataset
  dat_all <- left_join(dat_all, pca_results_avg)




# Calculate seasonal climate variables ------------------------------------
  
  ## Seedlings
  dat_all <- dat_all %>%
    mutate(Tmax_seas = apply( dplyr::select(., Tmax01:Tmax12), 1, sd))
  
  dat_all <- dat_all %>%
    mutate(Tave_seas = apply( dplyr::select(., Tave01:Tave12), 1, sd))
  
  dat_all <- dat_all %>%
    mutate(Tmin_seas = apply( dplyr::select(., Tmin01:Tmin12), 1, sd))
  
  ## Garden sites
  garden_climate <- garden_climate %>%
    mutate(Tmax_seas = apply( dplyr::select(., Tmax01:Tmax12), 1, sd))
  
  garden_climate <- garden_climate %>%
    mutate(Tave_seas = apply( dplyr::select(., Tave01:Tave12), 1, sd))
  
  garden_climate <- garden_climate %>%
    mutate(Tmin_seas = apply( dplyr::select(., Tmin01:Tmin12), 1, sd))
  
  
  ### Calculate max and min temperatature of the warmest month
  dat_all <- dat_all %>%
    mutate(Tmax = apply( dplyr::select(., Tmax01:Tmax12), 1, max))
  
  dat_all <- dat_all %>%
    mutate(Tmin = apply( dplyr::select(., Tmin01:Tmin12), 1, min))
  
  garden_climate <- garden_climate %>%
    mutate(Tmax = apply( dplyr::select(., Tmax01:Tmax12), 1, max))
  
  garden_climate <- garden_climate %>%
    mutate(Tmin = apply( dplyr::select(., Tmin01:Tmin12), 1, min))
  
  
  

# Choose climate variables ------------------------------------------------

  ## All climate variables - both averaged across year and by season, excluding radiaton  
  # climate_vars = colnames(
  #   #  dplyr::select(dat_all, Tmax_wt:RH,  ## Seasonal and annual variables
  #   dplyr::select(dat_all, MAT:RH, ## Just annual variables
  #                 #  dplyr::select(dat_all, Tmax_wt:RH_at, ## Seasonal variables
  #                 -contains("Rad"), - contains("MAR"), ## Getting rid of solar radiation
  #                 -contains("DD_0_sm"), -contains("PAS_sm"))) ## Do not vary at all..
  
  
  ### Climate vars from previous sork papers
  climate_vars = c("CMD", "MWMT", "MCMT", "DD5") 
  
  climate_vars_dif <- paste(climate_vars, "_dif", sep = "")
  
  x = 1
  
  ## Calculate difference in climate as Garden - Provenance
  for(var in climate_vars){
    
    dat_all[dat_all$site == "chico", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "chico", var]) - dat_all[dat_all$site == "chico", var]
    
    dat_all[dat_all$site == "placerville", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "placerville", var]) - dat_all[dat_all$site == "placerville", var]
    
    x = x + 1
  }
  
  
  ## Convert categorical variables to factors
  dat_all$site <- factor(dat_all$site)
  table(dat_all$site)
  
  dat_all$prov <- factor(dat_all$prov)
  table(dat_all$prov)
  
  dat_all$mom <- factor(dat_all$mom)
  table(dat_all$mom)
  
  dat_all$section <- factor(dat_all$section)
  table(dat_all$section)
  
  dat_all$block <- factor(dat_all$block)
  table(dat_all$block)
  
  ## Make sure values look OK
  # 
  # options(max.print = 10000)
  # 
  # numSummary(dat_all)
  # 
  # charSummary(dat_all)
  
  ## Correlations between variables
  
  # library(qtlcharts)
  # grpd <- dat_all %>%
  #   group_by(mom) %>% summarise_all(mean)
  # 
  # iplotCorr(dat_all[, climate_vars], reorder=FALSE)
  # 
  # 
  
  
  # Filtering data ----------------------------------------------------------
  
  
  ## Exclude trees that died of things like gophers or are shaded by trees
  dat_all <- dat_all[-which(dat_all$exclude == "y"), ]
  
  ## Change to NA outlier values in relative growth rate
  dat_all$rgr[dat_all$rgr <= quantile(dat_all$rgr, 0.01, na.rm = TRUE) |
                dat_all$rgr >= quantile(dat_all$rgr, 0.99, na.rm = TRUE)] <- NA
  
  ## Filter out individuals without an estimated RGR
  dat_all <- dplyr::filter(dat_all, !is.na(rgr))
  
  ## Filter out individuals without genetic data at the provenance level
  dat_all <- dplyr::filter(dat_all, !is.na(cluster_assigned))
  

# Scaling predictor variables -------------------------------------------------------

  dat_all_scaled <- dat_all
  
  x = 1

  climate_var_means <- NA
  climate_var_sds <- NA
  
  to_scale = c(climate_vars_dif, climate_vars, "height_2016")

  for(var in to_scale){

    climate_var_means[x] <- mean(dat_all_scaled[[var]], na.rm = TRUE)
    climate_var_sds[x] <- sd(dat_all_scaled[[var]], na.rm = TRUE)

    dat_all_scaled[var] <- (dat_all_scaled[var] - climate_var_means[x]) / climate_var_sds[x]

    x = x + 1
  }

  names(climate_var_means) <- to_scale
  names(climate_var_sds) <- to_scale


  summary(dat_all_scaled[, to_scale]) ## Means should all be 0


