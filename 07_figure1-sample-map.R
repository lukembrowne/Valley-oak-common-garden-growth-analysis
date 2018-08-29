# Source files ------------------------------------------------------------

  source("./01_clean-process-explore-data.R")
  source("./03_adding-climate-data.R")


# Load libraries ----------------------------------------------------------

  library(raster)
  library(rasterVis)



# Loading in rasters ------------------------------------------------------

## Read in elevation dem and create a hillshade for mapping
  dem <- raster("./data/gis/dem/ClimateNA_DEM_cropped.tif")
  
  slope = raster::terrain(dem, opt='slope')
  aspect = raster::terrain(dem, opt='aspect')
  hill = hillShade(slope, aspect, 40, 270)

## Load in california outline
  cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range <- readShapePoly("./data/gis/valley_oak_range/qlobata_refined.shp",
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  lobata_range_rough <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                      proj4string = CRS("+proj=longlat +datum=WGS84")) 
  
  
## Tmax summer
  tmax_rast_cropped <- raster("./data/gis/climate_data/BCM/cropped_tmax/tmax_sum_1951-1980_cropped.tif", proj4string = CRS("+proj=longlat +datum=WGS84"))
  plot(tmax_rast_cropped)
  
  
# Make SP points for samples  
  all_moms <- SpatialPoints(climate_garden_mom[, c("longitude", "latitude")], 
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  plot(all_moms)
  
  # Moms with gbs samples
  gbs_moms_sub <- climate_gbs_mom %>%
                dplyr::filter(climate_gbs_mom$accession %in% dat_all_scaled$accession)
  
  
  gbs_moms <- SpatialPoints(gbs_moms_sub[, c("longitude", "latitude")],
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  plot(gbs_moms)
  
  
  # All gbs samples
  gbs_moms_all <- SpatialPoints(climate_gbs_mom[, c("longitude", "latitude")],
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  plot(gbs_moms_all)
  
  # Sites of common gardens
  garden_sites <- SpatialPoints(garden_climate[, c("longitude", "latitude")],
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  plot(garden_sites)
  



# Make and save Plot ---------------------------------------------------------------

     pixel_num = 1e6
     hsTheme <- modifyList(GrTheme(), list(regions=list(alpha= .15)))
     tmaxTheme <- modifyList(BuRdTheme(), list(regions=list(alpha=1)))
    
     
  ## First plot showing garden moms   
    p1 = levelplot(tmax_rast_cropped, 
                  contour = FALSE, 
                  margin = FALSE, 
                  par.settings = tmaxTheme, 
                  maxpixels = pixel_num, 
                  colorkey = TRUE, 
                  main = "Maximum temperature summer")  +
      levelplot(hill,  margin = FALSE, par.settings=hsTheme) +
      latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
      latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 2)) +
      latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.75)) + 
      latticeExtra::layer(sp.points(all_moms, pch = 21, cex = 1,
                                    fill = "forestgreen", col = "grey10", alpha = .75)) +
      latticeExtra::layer(sp.points(gbs_moms, pch = 21, cex = 1,
                                    fill = "steelblue2", col = "grey10", alpha = .75))  +
      latticeExtra::layer(sp.points(garden_sites, pch = 22, cex = 2, 
                                    col = "grey10", fill = "orange")) 
    
    print(p1)
    
    png(paste0("./figs_tables/Figure S1_sample_map_garden_moms ", Sys.Date(), ".png"), res = 300,
        height = 2400, width = 1600)
    print(p1)
    dev.off()
    
    ## Smaller version for the legend
    png(paste0("./figs_tables/Figure S1_sample_map_garden_moms for legend", Sys.Date(), ".png"), res = 300,
        height = 800, width = 1600)
    print(p1)
    dev.off()
    
    
  ## Second plot showing GBS samples
    
    
    p2 = levelplot(tmax_rast_cropped, 
                   contour = FALSE, 
                   margin = FALSE, 
                   par.settings = tmaxTheme, 
                   maxpixels = pixel_num, 
                   colorkey = TRUE, 
                   main = "")  +
      levelplot(hill,  margin = FALSE, par.settings=hsTheme) +
      latticeExtra::layer(sp.polygons(cali_outline, lwd=1.5, col = "grey10")) + 
      latticeExtra::layer(sp.polygons(lobata_range_rough, col = "black", lwd = 2)) +
      latticeExtra::layer(sp.polygons(lobata_range, col = "grey50", lwd = 0.75)) + 
      latticeExtra::layer(sp.points(gbs_moms_all, pch = 21, cex = 1,
                                    fill = "mediumpurple1", col = "grey10", alpha = .75))
    
    print(p2)
    
    png(paste0("./figs_tables/Figure S1_sample_map_gbs_moms ", Sys.Date(), ".png"), res = 300,
        height = 2400, width = 1600)
    print(p2)
    dev.off()
    
    ## Smaller version for the legend
    png(paste0("./figs_tables/Figure S1_sample_map_gbs_moms for legend", Sys.Date(), ".png"), res = 300,
        height = 800, width = 1600)
    print(p2)
    dev.off()
  
    
    
    ## Histogram of temp transfer distances
    ggplot(dat_all_clim, aes(x = tmax_sum_dif)) +
      geom_histogram(fill = "grey95", col = "black") + 
      xlab("Tmax transfer distance") + ylab("Count") +
      geom_vline(aes(xintercept = 0), lty = 2) +
      geom_vline(xintercept = 1.1) + 
      geom_vline(xintercept = 4.88, lwd = 1.5) +
      theme_bw(15) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
