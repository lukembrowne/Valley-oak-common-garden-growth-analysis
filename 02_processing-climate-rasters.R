

## This script loads a series of raster files in a directory, crops them to california outline (more or less), projects them into lat / long, and saves them as a geotiff format

## Load required packages
library(raster)
library(sp)
library(maptools)


## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  
  
  
  
# Read in BCM data --------------------------------------------------------
  
  ## California Basin Characterization Model  
  ## https://ecologicalprocesses.springeropen.com/articles/10.1186/2192-1709-2-25
  
  
  ## Tmax
  bcm_tmax = raster("./data/gis/climate_data/BCM/tmx1951_1980_ave_HST_1513035565/tmx1951_1980_ave_HST_1513035565.tif")
  bcm_tmax = projectRaster(bcm_tmax, crs="+proj=longlat +datum=WGS84")
  
  ## Tmin
  bcm_tmin = raster("./data/gis/climate_data/BCM/tmn1951_1980_ave_HST_1513039895/tmn1951_1980_ave_HST_1513039895.tif")
  bcm_tmin = projectRaster(bcm_tmin, crs="+proj=longlat +datum=WGS84")
  
  ## CWD
  bcm_cwd = raster("./data/gis/climate_data/BCM/cwd1951_1980_ave_HST_1513035274/cwd1951_1980_ave_HST_1513035274.tif")
  bcm_cwd = projectRaster(bcm_cwd, crs="+proj=longlat +datum=WGS84")
  
  ## MAP
  bcm_map = raster("./data/gis/climate_data/BCM/ppt1951_1980_ave_HST_1513043823/ppt1951_1980_ave_HST_1513043823.tif")
  bcm_map = projectRaster(bcm_map, crs="+proj=longlat +datum=WGS84")
  
  ## Add variables to data frame
  dat_all_mom_avg <- dat_all %>%
    dplyr::select(mom, lat, lon) %>%
    group_by(mom) %>%
    summarise_all(mean) %>%
    mutate(bcm_tmax = raster::extract(bcm_tmax, .[, c("lon", "lat")]),
           bcm_tmin = raster::extract(bcm_tmin, .[, c("lon", "lat")]),
           bcm_cwd = raster::extract(bcm_cwd, .[, c("lon", "lat")]),
           bcm_map = raster::extract(bcm_map, .[, c("lon", "lat")]))
  
  dim(dat_all_mom_avg)
  
  

# Read in climateWNA data -----------------------------------------
  
  clim_wna <- read_csv("./data/cleaned_data/climateWNA data for all provenances 2017_12_11_Normal_1951_1980Y.csv")
  
  clim_wna <- clim_wna %>%
    rename(prov = id1, mom = ID2 ) %>%
    dplyr::select(-Latitude, -Longitude, -Elevation) %>%
    mutate(Tmax_monthly = rowMeans(dplyr::select(., Tmax01:Tmax12)),
           Tmin_monthly = rowMeans(dplyr::select(., Tmin01:Tmin12)))

  
#  clim_wna$mom <- as.factor(clim_wna$mom)
  
  ## Merge with bcm dataset
  
  dat_all_mom_avg <- left_join(dat_all_mom_avg, clim_wna)

  
  
  ## Tmax
  ggplot(dat_all_mom_avg, aes(x = Tmax_monthly, y = bcm_tmax)) + 
    geom_point(cex = 4, pch = 21, fill = "steelblue") + 
    xlab("climateWNA Tmax") + ylab("Flint BCM Tmax") + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, lty = 2) +
    theme_bw() + ggtitle("Tmax") + 
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  cor.test(dat_all_mom_avg$Tmax_monthly, dat_all_mom_avg$bcm_tmax)
  
  ## Tmin
  ggplot(dat_all_mom_avg, aes(x = Tmin_monthly, y = bcm_tmin)) + 
    geom_point(cex = 4, pch = 21, fill = "steelblue") + 
    xlab("climateWNA Tmin") + ylab("Flint BCM Tmin") + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, lty = 2) +
    theme_bw() + ggtitle("Tmin") + 
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  cor.test(dat_all_mom_avg$Tmin_monthly, dat_all_mom_avg$bcm_tmin)
  
  
  
  
  ## CMD
  ggplot(dat_all_mom_avg, aes(x = CMD, y = bcm_cwd)) + 
    geom_point(cex = 4, pch = 21, fill = "steelblue") + 
    xlab("climateWNA CMD") + ylab("Flint BCM CMD") + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, lty = 2) +
    theme_bw() +  ggtitle("CMD") + 
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  cor.test(dat_all_mom_avg$CMD, dat_all_mom_avg$bcm_cwd)
  
  ## Precipitation
  ggplot(dat_all_mom_avg, aes(x = MAP, y = bcm_map)) + 
    geom_point(cex = 4, pch = 21, fill = "steelblue") + 
    xlab("climateWNA MAP") + ylab("Flint BCM MAP") + 
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, lty = 2) +
    theme_bw() +  ggtitle("MAP") + 
    theme(plot.margin = unit(c(1,1,1,1), "cm"), 
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  cor.test(dat_all_mom_avg$MAP, dat_all_mom_avg$bcm_map)
  

# Crop climateWNA rasters -------------------------------------------------


## project to lamber conformal conic to use for cropping below - transforming this and then cropping speeds up operations a lot rather than transforming raster and cropping that

## CRS comes from read in .nc file with raster() and copying the auto-read CRS, adjusting for a typo
  cali_outline_lcc <- spTransform(cali_outline, CRS("+proj=lcc +lat_1=49.0 +lat_2=77.0 +lat_0=0.0 +lon_0=-95.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


## Read in climate variables - https://adaptwest.databasin.org/pages/adaptwest-climatena

  dir_name <- "./data/gis/climate_data/NA_NORM_8110_Monthly_netCDF/"
  
  raster_files <- list.files(dir_name, full.names = TRUE)
  raster_files <- raster_files[grep("*.nc", raster_files)] # Only those with nc extension


## Loop over all rasters in the directory
  for(file in raster_files){
    
    cat("Working on file:", file, "... \n")
    
    ## Load in raster file
    raster_temp <- raster(file)
    
    ## Assign correct CRS since there is a typo in the NetCDF file
    crs(raster_temp) <- "+proj=lcc +lat_1=49.0 +lat_2=77.0 +lat_0=0.0 +lon_0=-95.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    ## Crop
    raster_temp_cropped <- raster::crop(raster_temp, cali_outline_lcc)
    
    ## Masking to get exact shape
    raster_temp_cropped <- mask(raster_temp_cropped, cali_outline_lcc)
    
    ## Project to latlong
    raster_cropped_latlong <- projectRaster(raster_temp_cropped, crs="+proj=longlat +datum=WGS84")
    
    ## Write to file - geotiff format set by setting file extension
    
    file_new = gsub(".nc", "_cropped.tif", file)
    
    writeRaster(raster_cropped_latlong, file_new, overwrite = TRUE)
    
  }
  
  


  
  

