# Load libraries ----------------------------------------------------------

library(tidyverse)
library(raster)
library(sp)
library(maptools)
library(httr)
library(stringr)

#install.packages("climates",,'http://rforge.net/',type='source', dependencies = TRUE)
library(climates)



# Download future climate scenarios ---------------------------------------

  ## Future climate files downloaded from:   https://geo.pointblue.org/commonsmap/index.php?ds=1129
    
    # Choose climate scenario
     scenario = "CCSM4_rcp85"
    
    # Climate variables in bcm format
      bcm_clim_vars <- c("aet", "cwd", "ppt", "tmn", "tmx")
    
    # Time periods - including water year (wy), months, and seasons (jja, djf)
      time_periods <- c("wy", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep",
                        "oct", "nov", "dec", "jja", "djf")
  
    # Directory to save output files  
    output_dir <- paste0("./data/gis/climate_data/BCM/future/", scenario, "/")
    
    errors <- NULL
  
  # Loop through climate vars    
  for(climate_var in bcm_clim_vars){
    
    # Loop through time periods / months / seasons
    for(time_period in time_periods){
      
      cat("\n Working on:", climate_var, " in ", time_period, " ...\n")
      
    ## Combinations to skip  
      if((climate_var == "aet" | climate_var == "cwd") & time_period != "wy"){
        cat("Skipping this combination...\n")
        next
      }
      
      if(climate_var == "tmn" & time_period == "jja"){
        cat("Skipping this combination...\n")
        next
      }
      
      if(climate_var == "tmx" & time_period == "djf"){
        cat("Skipping this combination...\n")
        next
      }
      
      if(climate_var == "ppt" & time_period %in% c("jja", "djf")){
        cat("Skipping this combination...\n")
        next
      } 
      
      
   # Since water years are labeled a bit differently, need to set up two separate labels for water year   
    if(time_period != "wy"){  
      time_period1 <- time_period
      time_period2 <- time_period
    } else {
      time_period1 <- ""
      time_period2 <- "wy"
    }
      
      # Make URL to file - example from the interactive BCM map (https://geo.pointblue.org/commonsmap/index.php?ds=1129) and then copying link address on the "Download data" button
      url <- paste0("https://geo.pointblue.org/commonsmap/add2.php?ds=1129&commonsuser=lukembrowne::lukembrowne@gmail.com&guid=bb94c436-12bd-4b82-a2ee-479b2e01205e&raster_path=/mnt/data/tif/BCM2014/&srid=3310&active_raster=", climate_var, "2070_2099", time_period1, "_ave_CCSM4_rcp85&width=null&height=null&llx=-433493.711264374&lly=-623577.7967801094&urx=654964.4482971333&ury=613350.9509733636&variable=flint_30yr&stat=ave&bounds=-13869908.243,3803287.923,-12585021.03,5366378.461&model=", scenario, "&year=2070_2099&climatevar=", climate_var, "&bcm_style=flint_", climate_var, "_ave&legend=flint_", climate_var, "_ave&month=", time_period2, "&raster=", climate_var, "2070_2099", time_period1, "_ave_", scenario)
      
    # Download file
      file <- GET(url, write_disk(paste0(output_dir, "tmp.file"), 
                                  overwrite = TRUE), progress())
      
    ## Check to make sure it's a big enough file
      size = as.numeric(headers(file)$`content-length`) / 1000000
      cat("Size of file is: ", round(size, 2), " mb... \n")
      if(size < 5){
        cat("**** Error: size of file is < 5mb... \n\n")
      errors <- c(errors, paste0(climate_var, "--", time_period))
      }
      
    ## Rename file  
      filename <- str_match(headers(file)$`content-disposition`, "\"(.*)\"")[2]
      # rename
      file.rename(paste0(output_dir, "tmp.file"), paste0(output_dir, filename)) 
      
    ## Unzip the file to folder
      unzip(paste0(output_dir, "/", filename),
            exdir = paste0(output_dir, gsub(pattern = ".zip", "", filename)))
    
    ## Remove .zip file
      file.remove(paste0(output_dir, filename))
      
    } # End time period loop
    
  } # End climate variable loop 
    
  errors
   
  
  
   
  
# Read in and process and clean raster files ------------------------------

  
  ## Directory path to where 1951-1980 historical BCM climate data is located
  dir_name_hist <- "./data/gis/climate_data/BCM/historical/1951-1980/"
  
  
  dir_name_future <- "./data/gis/climate_data/BCM/future/"
  

## Creates vector with list of file paths to all .tif raster files
  raster_files_hist <- list.files(dir_name_hist, full.names = TRUE, recursive = TRUE)
  raster_files_hist <- raster_files_hist[grep("*[[:digit:]].tif$", raster_files_hist)] # Only those with .tif extension
  raster_files_hist
  
  raster_files_future <- list.files(dir_name_future, full.names = TRUE, recursive = TRUE)
  raster_files_future <- raster_files_future[grep("*[[:digit:]].tif$", raster_files_future)] # Only those with .tif extension
  raster_files_future
  
## Load into raster stacks
  raster_hist <- raster::stack(raster_files_hist)
  raster_hist
  
  raster_future <- raster::stack(raster_files_future)
  raster_future
  
## Project to lat long - maybe not necessary until the end
  # raster_hist <- raster::projectRaster(raster_hist, crs="+proj=longlat +datum=WGS84") 
  # raster_future <- raster::projectRaster(raster_future, crs="+proj=longlat +datum=WGS84")
  

## Rename variables
  
  ## HISTORICAL
  
  ## Remove second half of colname - e.g. "_ave_HST_1513043823"
  names(raster_hist) <- unlist(lapply(strsplit(x = names(raster_hist), split = "_ave_"), '[[', 1))
  
  ## Take out years
  names(raster_hist) <- gsub("1951_1980", "", names(raster_hist))
  
  ## Add month number to monthly variables
  names(raster_hist) <- gsub("jan", "_jan", names(raster_hist))
  names(raster_hist) <- gsub("feb", "_feb", names(raster_hist))
  names(raster_hist) <- gsub("mar", "_mar", names(raster_hist))
  names(raster_hist) <- gsub("apr", "_apr", names(raster_hist))
  names(raster_hist) <- gsub("may", "_may", names(raster_hist))
  names(raster_hist) <- gsub("jun", "_jun", names(raster_hist))
  names(raster_hist) <- gsub("jul", "_jul", names(raster_hist))
  names(raster_hist) <- gsub("aug", "_aug", names(raster_hist))
  names(raster_hist) <- gsub("sep", "_sep", names(raster_hist))
  names(raster_hist) <- gsub("oct", "_oct", names(raster_hist))
  names(raster_hist) <- gsub("nov", "_nov", names(raster_hist))
  names(raster_hist) <- gsub("dec", "_dec", names(raster_hist))
  
  ## For seasonal variables
  names(raster_hist) <- gsub("djf", "_winter", names(raster_hist))
  names(raster_hist) <- gsub("jja", "_sum", names(raster_hist))
  
  ## Replace tmn with tmin and tmx with tmax
  names(raster_hist) <- gsub("tmn", "tmin", names(raster_hist))
  names(raster_hist) <- gsub("tmx", "tmax", names(raster_hist))
  
  names(raster_hist)
  
  
  ## FUTURE
  ## Remove second half of colname - e.g. "_ave_HST_1513043823"
  names(raster_future) <- unlist(lapply(strsplit(x = names(raster_future), split = "_ave_"), '[[', 1))
  
  ## Take out years
  names(raster_future) <- gsub("2070_2099", "", names(raster_future)) ## Will have to change d
  
  ## Add month number to monthly variables
  names(raster_future) <- gsub("jan", "_jan", names(raster_future))
  names(raster_future) <- gsub("feb", "_feb", names(raster_future))
  names(raster_future) <- gsub("mar", "_mar", names(raster_future))
  names(raster_future) <- gsub("apr", "_apr", names(raster_future))
  names(raster_future) <- gsub("may", "_may", names(raster_future))
  names(raster_future) <- gsub("jun", "_jun", names(raster_future))
  names(raster_future) <- gsub("jul", "_jul", names(raster_future))
  names(raster_future) <- gsub("aug", "_aug", names(raster_future))
  names(raster_future) <- gsub("sep", "_sep", names(raster_future))
  names(raster_future) <- gsub("oct", "_oct", names(raster_future))
  names(raster_future) <- gsub("nov", "_nov", names(raster_future))
  names(raster_future) <- gsub("dec", "_dec", names(raster_future))
  
  ## For seasonal variables
  names(raster_future) <- gsub("djf", "_winter", names(raster_future))
  names(raster_future) <- gsub("jja", "_sum", names(raster_future))
  
  ## Replace tmn with tmin and tmx with tmax
  names(raster_future) <- gsub("tmn", "tmin", names(raster_future))
  names(raster_future) <- gsub("tmx", "tmax", names(raster_future))
  
  names(raster_future)
  

## Will need a section here to calculate bioclimatic variables in future climates
  
  
## Calculate difference between future and historical  
  raster_dif <- raster_future # Initialize
  
  for(var in names(raster_dif)){
    cat("Working on", var, "...\n")
    raster_dif[[var]] <- raster_future[[var]] - raster_hist[[var]]
  }
  
  names(raster_dif) <- paste0(names(raster_dif), "_dif") ## Add _dif to raster names
  

# Extract and scale values ------------------------------------------------
  
  ## Extract values from raster into dataframe
  extracted_vals = raster::extract(raster_dif, 1:ncell(raster_dif))
  extracted_vals <- as.data.frame(extracted_vals)
  
  ### Scale extracted values to match with scaled variables
  extracted_vals_scaled <- extracted_vals
  
  for(var in climate_vars_dif){
    if(var %in% colnames(extracted_vals_scaled)){
    extracted_vals_scaled[, var] <- (extracted_vals[, var] - scaled_var_means[var]) / 
      scaled_var_sds[var]
    }
  }
  
  summary(extracted_vals_scaled)
  
  
  
## Testing out predicting breeding values from raster maps
  
  ## Extract values from raster into dataframe
  extracted_vals = raster::extract(raster_hist[["tmax_sum"]], 
                                   1:ncell(raster_hist[["tmax_sum"]]))
  extracted_vals <- as.data.frame(extracted_vals)
  
  ### Scale extracted values to match with scaled variables
  extracted_vals_scaled <- extracted_vals
  
  var = "tmax_sum"
      extracted_vals_scaled <- (extracted_vals - scaled_var_means[var]) / 
        scaled_var_sds[var]
  
  summary(extracted_vals_scaled)
  
  bvs <- bglr_test$ETA[[1]]$u
  gebv_lm <- lm(bvs ~ y)
  
  
  
  gebvs <- predict(gebv_lm, newdata = data.frame(y = extracted_vals_scaled[, 1]))
  
  raster_gebvs <- raster_hist[["tmax_sum"]]
  values(raster_gebvs) <- gebvs
  
  levelplot(raster_gebvs)
  
  

  cor.test(values(raster_hist[["tmax_sum"]]), values(raster_dif[["tmax_sum_dif"]]), na.rm = T)

  
  ran <- sample(1:1500000, size = 300000)
  
  plot(values(raster_hist[["tmax_sum"]])[ran],
       values(raster_dif[["tmax_sum_dif"]])[ran], pch = ".")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
