
## This script loads a series of raster files in a directory, crops them to california outline (more or less), projects them into lat / long, and saves them as a geotiff format

  ## Load required packages
  library(tidyverse)
  library(raster)
  library(rasterVis)
  library(sp)
  library(maptools)
  
  #install.packages("climates",,'http://rforge.net/',type='source', dependencies = TRUE)
  library(climates)
  


#### Run code below if needed to extract and write climate data to file


  
# Read in Historical BCM data --------------------------------------------------------
  
  ## California Basin Characterization Model  
  ## https://ecologicalprocesses.springeropen.com/articles/10.1186/2192-1709-2-25
  

## Get locations of maternal trees in common garden
  dat_mom <- gs_read(gs_key("1DUEV-pqV28D6qJl6vJM6S1VLWbxc9E71afufulDRbIc"), ws = 2)
  
  ## Subset down to those with accession numbers
  dat_mom <- dat_mom[!is.na(dat_mom$Accession),]
  dim(dat_mom) ## Should be 659
  
 
 ## For calculating bioclim variables across full rasters
  # Read in raster, convert to lat long, and the make DF with lat long coordinates - HUGE DF
  tmax_rast <- raster("./data/gis/climate_data/BCM/historical/1951-1980/tmx1951_1980jja_ave_HST_1513103038/tmx1951_1980jja_ave_HST_1513103038.tif")
  tmax_rast <- raster::projectRaster(tmax_rast, crs="+proj=longlat +datum=WGS84") # Project
  dat_mom <- xyFromCell(tmax_rast, 1:ncell(tmax_rast)) # Extract lat long coordinates
  colnames(dat_mom) <- c("Longitude", "Latitude") # Rename cell names
  dat_mom <- as.data.frame(dat_mom) # Convert to dataframe
  dat_mom$cell_id <- 1:nrow(dat_mom) # Make column for cell ids
  dat_mom <- dat_mom[!is.na(values(tmax_rast)), ] # Filter down to not NAs
  
  head(dat_mom)
  dim(dat_mom)
  
  ## Directory path to where 1951-1980 historical BCM climate data is located
  dir_name <- "./data/gis/climate_data/BCM/historical/1951-1980/"
  
  ## Directory path to where 1921-1950 historical BCM climate data is located
  # dir_name <- "./data/gis/climate_data/BCM/historical/1921-1950/"
  
    ## Creates vector with list of file paths to all .tif raster files
    raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    raster_files <- raster_files[grep("*[[:digit:]].tif$", raster_files)] # Only those with .tif extension
    

  ## Loop through raster files and add climate values to dat_mom dataframe
  for(file in raster_files){
    
    cat("Working on file:", file, "...\n")
    
    ## Load in raster
    raster_temp <- raster::stack(file)
    
    ## For 1921-1950 data because some files are missing their CRS
    # raster_temp <- raster::raster(x = file, 
    #                              crs = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 
    
    ## Project to lat long
   raster_latlong <- raster::projectRaster(raster_temp, crs="+proj=longlat +datum=WGS84") 
    
    ## Extract climate values to dataframe
    col_name <-  names(raster_temp) ## Use for adding back to dataframe

    dat_mom[, col_name] <-  raster::extract(raster_latlong, dat_mom[, c("Longitude", "Latitude")])
    
  } ## End file loop
 
  head(dat_mom)
  
  
  
  
  
  
  ## Save as backup just in case
    dat_mom2 <- dat_mom
  
  ## Rename columns
  
  ## Remove second half of colname - e.g. "_ave_HST_1513043823"
  colnames(dat_mom) <- unlist(lapply(strsplit(x = colnames(dat_mom), split = "_ave_"), '[[', 1))
  
  ## Take out years
  colnames(dat_mom) <- gsub("1951_1980", "", colnames(dat_mom))
  
  ## Add month number to monthly variables
  colnames(dat_mom) <- gsub("jan", "_jan", colnames(dat_mom))
  colnames(dat_mom) <- gsub("feb", "_feb", colnames(dat_mom))
  colnames(dat_mom) <- gsub("mar", "_mar", colnames(dat_mom))
  colnames(dat_mom) <- gsub("apr", "_apr", colnames(dat_mom))
  colnames(dat_mom) <- gsub("may", "_may", colnames(dat_mom))
  colnames(dat_mom) <- gsub("jun", "_jun", colnames(dat_mom))
  colnames(dat_mom) <- gsub("jul", "_jul", colnames(dat_mom))
  colnames(dat_mom) <- gsub("aug", "_aug", colnames(dat_mom))
  colnames(dat_mom) <- gsub("sep", "_sep", colnames(dat_mom))
  colnames(dat_mom) <- gsub("oct", "_oct", colnames(dat_mom))
  colnames(dat_mom) <- gsub("nov", "_nov", colnames(dat_mom))
  colnames(dat_mom) <- gsub("dec", "_dec", colnames(dat_mom))
  
  ## For seasonal variables
  colnames(dat_mom) <- gsub("djf", "_winter", colnames(dat_mom))
  colnames(dat_mom) <- gsub("jja", "_sum", colnames(dat_mom))
  
  ## Replace tmn with tmin and tmx with tmax
  colnames(dat_mom) <- gsub("tmn", "tmin", colnames(dat_mom))
  colnames(dat_mom) <- gsub("tmx", "tmax", colnames(dat_mom))
  
  names(dat_mom)
  
  ## Reorder columns so monthly variables are sequential
  dat_mom <- dat_mom[, c(
                      #  "Locality full name", "Locality", "Sample #","Accession", 
                         "cell_id", "Latitude", "Longitude", # If working with rasters
                      #  "ID", "Latitude", "Longitude", "Elevation", 
                         "aet",
                         "cwd", "ppt", "tmin_winter", "tmax_sum", "tmax", "tmin",
                         "ppt_jan","ppt_feb","ppt_mar","ppt_apr","ppt_may","ppt_jun",
                         "ppt_jul","ppt_aug","ppt_sep","ppt_oct","ppt_nov","ppt_dec",
                         "tmax_jan","tmax_feb","tmax_mar","tmax_apr","tmax_may","tmax_jun",
                         "tmax_jul","tmax_aug","tmax_sep","tmax_oct","tmax_nov","tmax_dec",
                         "tmin_jan","tmin_feb","tmin_mar","tmin_apr","tmin_may","tmin_jun",
                         "tmin_jul","tmin_aug","tmin_sep","tmin_oct","tmin_nov","tmin_dec")]
  

# Calculate Bioclimatic variables for each maternal tree ------------------

  bioclim_vars <-  climates::bioclim2(tmin = as.data.frame(dplyr::select(dat_mom, tmin_jan:tmin_dec)),
                     tmax = as.data.frame(dplyr::select(dat_mom, tmax_jan:tmax_dec)),
                     prec = as.data.frame(dplyr::select(dat_mom, ppt_jan:ppt_dec)),
                     files.as.inputs = FALSE)
  
  dat_mom <- bind_cols(dat_mom, as.data.frame(bioclim_vars))
  
  str(dat_mom)
  

  ## For common garden moms
  ## Reorder columns
  dat_mom_final <- dat_mom %>%
    dplyr::select(# `Locality`,Accession, Latitude, Longitude, `Elevation (m)`, 
                  "cell_id", "Latitude", "Longitude", # If working with rasters
                  #`ID`, Latitude, Longitude, Elevation, # For pauls samples
                  tmax, tmax_sum, tmin, tmin_winter,
                  ppt, cwd, aet, bioclim_01:bioclim_19, tmax_jan:tmax_dec, 
                  tmin_jan:tmin_dec, ppt_jan:ppt_dec)
  
  ## Write to file
  
  # write_csv(dat_mom_final, path = "./data/cleaned_data/GBS tree climate data BCM 1951-1980 2018_03_08.csv")
  
  
## Save back as raster

  # Bioclim_04 - temperature seasonality
    rast_bioclim_04 <- tmax_rast
    values(rast_bioclim_04) <- NA
  
    values(rast_bioclim_04)[dat_mom_final$cell_id] <- dat_mom_final$bioclim_04
    levelplot(rast_bioclim_04, margin = FALSE)
    
    writeRaster(rast_bioclim_04, "./data/gis/climate_data/BCM/historical/1951-1980/bioclim_04", 
                format = "GTiff")
    
  # Bioclim_15 - Precipitation seasonality
    rast_bioclim_15 <- tmax_rast
    values(rast_bioclim_15) <- NA
    
    values(rast_bioclim_15)[dat_mom_final$cell_id] <- dat_mom_final$bioclim_15
    levelplot(rast_bioclim_15, margin = FALSE)
    
    writeRaster(rast_bioclim_15, "./data/gis/climate_data/BCM/historical/1951-1980/bioclim_15", 
                format = "GTiff")
    
  # Bioclim_18 - Precipitation of warmest quarter
    rast_bioclim_18 <- tmax_rast
    values(rast_bioclim_18) <- NA
    
    values(rast_bioclim_18)[dat_mom_final$cell_id] <- dat_mom_final$bioclim_18
    levelplot(rast_bioclim_18, margin = FALSE)
    
    writeRaster(rast_bioclim_18, "./data/gis/climate_data/BCM/historical/1951-1980/bioclim_18", 
                format = "GTiff")
    
  # Bioclim_19 - Precipitation of coolest quarter
    rast_bioclim_19 <- tmax_rast
    values(rast_bioclim_19) <- NA
    
    values(rast_bioclim_19)[dat_mom_final$cell_id] <- dat_mom_final$bioclim_19
    levelplot(rast_bioclim_19, margin = FALSE)
    
    writeRaster(rast_bioclim_19, "./data/gis/climate_data/BCM/historical/1951-1980/bioclim_19", 
                format = "GTiff")  
  

    

# Read in current BCM climate data for common garden sites ----------------------------------------

    dat_garden <- data.frame(site = c("Chico", "IFG"),
                          lat = c(39.710764, 38.740382),
                          lon = c(-121.786064, -120.735832))
    
    ## Directory path to where monthly current BCM climate data is located
    dir_name <- "./data/gis/climate_data/BCM/current/"
    
    
    ## If files haven't already been projected to lat long
    
    # ## Creates vector with list of file paths to all .tif raster files
    # raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    # raster_files <- raster_files[grep("*.tif$", raster_files)] # Only those with .tif extension
    # 
    # ## Loop through raster files, convert to lat long if needed so we can extract values
    # for(file in raster_files){
    #   
    #   cat("Working on file:", file, "...\n")
    #   
    #   ## Load in raster
    #   raster_temp <- raster::stack(file) 
    #   
    #   ## Project to lat long and save to file
    #   raster::projectRaster(raster_temp, crs="+proj=longlat +datum=WGS84",
    #                        file = gsub(".tif", "_latlong.tif", file)) 
    # 
    # } ## End file loop
    # 
    
    
    ## Creates vector with list of file paths to all .tif raster files projected to lat long
    raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    raster_files <- raster_files[grep("*_latlong.tif$", raster_files)] # Only those with .tif extension

    ## Loop through raster files
    for(file in raster_files){
      
      cat("Working on file:", file, "...\n")
      
      raster_temp <- raster::stack(file) 
      
      ## Extract climate values to dataframe
      dat_garden <- cbind(dat_garden, raster::extract(raster_temp, dat_garden[, c("lon", "lat")]))
    
    }
    
    str(dat_garden)
    
    ## Remove latlong from column names
    colnames(dat_garden) <- gsub(pattern = "latlong.", "", colnames(dat_garden))
    colnames(dat_garden)
    
    ## Months / bands are based on water years, which begin Oct 1st and end Sep 30
    ## So band/month 1 = Oct, 2 = Nov, 3 = Dec, etc
    ## Need to convert column names
    colnames(dat_garden) <- gsub(pattern = "_1$", "_oct", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_2", "_nov", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_3", "_dec", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_4", "_jan", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_5", "_feb", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_6", "_mar", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_7", "_apr", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_8", "_may", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_9", "_jun", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_10", "_jul", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_11", "_aug", colnames(dat_garden))
    colnames(dat_garden) <- gsub(pattern = "_12", "_sep", colnames(dat_garden))
    
    ## Replace tmn with tmin and tmx with tmax
    colnames(dat_garden) <- gsub("tmn", "tmin", colnames(dat_garden))
    colnames(dat_garden) <- gsub("tmx", "tmax", colnames(dat_garden))
    
    colnames(dat_garden)

    
    
    ## Average values across years

    ## Function that takes variable abbreviation, average across years within the 
    ## same month and outputs table
    avg_data <- function(dat_garden, var, method){
      
        month_order <- c("jan", "feb", "mar", "apr", "may", "jun", 
                         "jul", "aug", "sep", "oct", "nov", "dec")
    
       temp =  dat_garden %>%
          dplyr::select(site, contains(var)) %>% ## Select only one variable
          tidyr::gather(key = year_mo, value = climate_var, contains(var)) %>% ## Turn into long format
          dplyr::mutate(year_mo =  gsub(pattern = paste0(var,"_wy"),
                                        "", .$year_mo)) %>% ## Remove extra text
          tidyr::separate(year_mo, into = c("year", "month")) %>% ## Separate year mo into cols
         dplyr::mutate(month = forcats::fct_relevel(month, month_order))
       
       
       ## Separate out values by month, then average across years
       if(method == "years_by_month"){
         
         temp <- temp %>%
              dplyr::group_by(site, month) %>%
              dplyr::select(-year) %>% ## Remove year column
              dplyr::summarise_all("mean") %>% ## Average across years by month
              tidyr::spread(key = month, value = climate_var)
       
          colnames(temp)[-1] <- paste0(var, "_", colnames(temp)[-1])
          return(temp)
       }
       
       ## Sum up value within each year and then average across years
       if(method == "sum_across_year"){
         temp <- temp %>%
           dplyr::select(-month) %>%
           dplyr::group_by(site, year) %>%
           dplyr::summarise_all("sum") %>%
           dplyr::select(-year) %>%
           dplyr::group_by(site) %>%
           dplyr::summarise_all("mean")
         
         colnames(temp)[2] <- var
         return(temp)
         
       }
       
       ## Average values across months within year, then average across years
       if(method == "avg_across_year"){
         
         temp <- temp %>%
           dplyr::select(-month) %>%
           dplyr::group_by(site, year) %>%
           dplyr::summarise_all("mean") %>%
           dplyr::select(-year) %>%
           dplyr::group_by(site) %>%
           dplyr::summarise_all("mean")
         
         colnames(temp)[2] <- var
         
         return(temp)
         
       }
       
       ## Average values across June-August, then average across years
       if(method == "avg_across_summer"){
         
         temp <- temp %>%
           dplyr::filter(month %in% c("jun", "jul", "aug")) %>%
           dplyr::select(-month) %>%
           dplyr::group_by(site, year) %>%
           dplyr::summarise_all("mean") %>%
           dplyr::select(-year) %>%
           dplyr::group_by(site) %>%
           dplyr::summarise_all("mean")
         
         colnames(temp)[2] <- paste0(var, "_sum")
         
         return(temp)
         
       }
       
       ## Average values across December - February then average across years
       if(method == "avg_across_winter"){
         
         temp <- temp %>%
           dplyr::filter(month %in% c("dec", "jan", "feb")) %>%
           dplyr::select(-month) %>%
           dplyr::group_by(site, year) %>%
           dplyr::summarise_all("mean") %>%
           dplyr::select(-year) %>%
           dplyr::group_by(site) %>%
           dplyr::summarise_all("mean")
         
         colnames(temp)[2] <- paste0(var, "_winter")
         
         return(temp)
         
       }
       
       
     
    }
    
   ## Max temp 
   tmax_by_month =  avg_data(dat_garden, "tmax", method = "years_by_month") 
   tmax_by_month
   
   ## Min temp
   tmin_by_month =  avg_data(dat_garden, "tmin", method = "years_by_month")
   tmin_by_month
   
   ## Precipitation
   ppt_by_month = avg_data(dat_garden, "ppt", method = "years_by_month")
   ppt_by_month
   
   
   ## Averaged across year
   
   ## Tmax
   tmax_year =  avg_data(dat_garden, "tmax", method = "avg_across_year") 
   tmax_year
   
   ## Tmax
   tmin_year =  avg_data(dat_garden, "tmin", method = "avg_across_year") 
   tmin_year
   
   
  ## Average across summer
   ## Tmax
   tmax_sum = avg_data(dat_garden, "tmax", method = "avg_across_summer") 
   tmax_sum
   
   ## Average across winter
   ## Tmax
   tmin_winter = avg_data(dat_garden, "tmin", method = "avg_across_winter") 
   tmin_winter
   
   
   ## Summed across year
   
     ## PPT
     ppt_year <- avg_data(dat_garden, "ppt", method = "sum_across_year")
     ppt_year
     
     ## CWD
     cwd_year <- avg_data(dat_garden, "cwd", method = "sum_across_year")
     cwd_year
     
     ## AET
     aet_year <- avg_data(dat_garden, "aet", method = "sum_across_year")
     aet_year
     
  
  
# Calculate Bioclimatic variables for each garden site ------------------
    
    bioclim_vars_garden <-  climates::bioclim2(tmin = as.data.frame(tmin_by_month[, -1]),
                                               ## Remove site column
                                        tmax = as.data.frame(tmax_by_month[, -1]),
                                        prec = as.data.frame(ppt_by_month[, -1]),
                                        files.as.inputs = FALSE)
    
   
     
     ### Bring in all variables into one dataframe
    dat_garden_final <- dplyr::bind_cols(dat_garden[, c("site", "lat", "lon")], 
                                   as.data.frame(bioclim_vars_garden))   
    
    dat_garden_final <- left_join(dat_garden_final, tmax_by_month ) %>%
                       left_join(., tmin_by_month) %>%
                       left_join(., ppt_by_month) %>%
                       left_join(., tmax_year) %>%
                        left_join(., ppt_year) %>%
                        left_join(., cwd_year) %>%
                        left_join(., aet_year) %>%
                        left_join(., tmin_year) %>%
                        left_join(., tmax_sum) %>%
                        left_join(., tmin_winter)
   
    glimpse(dat_garden_final)

    ## Reorder columns
    dat_garden_final <- dat_garden_final %>%
      dplyr::select(site, lat, lon, tmax, tmax_sum, tmin, tmin_winter,
                    ppt, cwd, aet, bioclim_01:bioclim_19, tmax_jan:tmax_dec, 
                    tmin_jan:tmin_dec, ppt_jan:ppt_dec)
    
    ## Write to file
    
    # write_csv(dat_garden_final, 
    # path = "./data/cleaned_data/common garden climate data BCM 2014-2016 2018_03_08.csv")

  

# Extract values for future climate scenarios -----------------------------

    ## Get locations of maternal trees in common garden
    dat_mom <- gs_read(gs_key("1DUEV-pqV28D6qJl6vJM6S1VLWbxc9E71afufulDRbIc"), ws = 2)
    
    ## Subset down to those with accession numbers
    dat_mom <- dat_mom[!is.na(dat_mom$Accession),]
    dim(dat_mom) ## Should be 659
    
    
    ## Directory path where future climate scenarios are located
    dir_name <- "./data/gis/climate_data/BCM/future/"
    
    ## Creates vector with list of file paths to all .tif raster files
    raster_files <- list.files(dir_name, full.names = TRUE, recursive = TRUE)
    raster_files <- raster_files[grep("*[[:digit:]].tif$", raster_files)] # Only those with .tif extension
    raster_files <- raster_files[grep("jja_ave", raster_files)] # Only tmax sum files
    
    raster_files
    
    
    ## Loop through raster files and add climate values to dat_mom dataframe
    for(file in raster_files){
      
      cat("Working on file:", file, "...\n")
      
      ## Load in raster
      raster_temp <- raster::stack(file)
      
      ## Project to lat long
      raster_latlong <- raster::projectRaster(raster_temp, crs="+proj=longlat +datum=WGS84") 
      
      ## Extract climate values to dataframe
      col_name <-  names(raster_temp) ## Use for adding back to dataframe
      
      dat_mom[, col_name] <-  raster::extract(raster_latlong, dat_mom[, c("Longitude", "Latitude")])
      
    } ## End file loop
    
    head(dat_mom)
    
    ## Save as backup just in case
    dat_mom2 <- dat_mom
    
    ## Rename columns
    
    ## Replace tmn with tmin and tmx with tmax
    colnames(dat_mom) <- gsub("tmx", "tmax_sum", colnames(dat_mom))
    
    names(dat_mom)
    
    ## Write to file
    
    # write_csv(dat_mom, path = "./data/cleaned_data/maternal tree FUTURE climate data BCM 2018_06_18.csv")
    


