## TODO


# Load libraries ----------------------------------------------------------
  library(tidyverse)
  library(googlesheets)
  library(maptools)
  library(qtlcharts)


# Load and quality check data ---------------------------------------------------------------

    ## Load in california outline
    cali_outline <- readShapePoly("./data/gis/california_outline/california_outline.shp",
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    lobata_range <- readShapePoly("./data/gis/valley_oak_range/querloba.shp",
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
    


  ## Load in garden data from 2017 - "Qlobata census measurements 2017"

    ## If token is stale: gs_auth(new_user = TRUE)
    dat_17_raw <- gs_read(gs_key("1CvaUjE7qBvKbbc8QaG_B4vxFZlCNe79z_q9JwYhMhhw"), ws = 2)
    colnames(dat_17_raw) <- paste0(colnames(dat_17_raw), "_2017")
    
    ## Remove if we don't know the accession / progeny
    dat_17_raw <- dat_17_raw[!is.na(dat_17_raw$Accession_progeny_2017),]
    
    glimpse(dat_17_raw) 
    dim(dat_17_raw)
    
  ## Load in garden data from 2015 - "Qlobata census measurements 2015"
    dat_15_raw <- gs_read(gs_key("1ljUeS4j3tte1qaXQFNJGmeOL7wTIlBRC8N5vJ1HlN8k"), ws = 2)
    colnames(dat_15_raw) <- paste0(colnames(dat_15_raw), "_2015")
    
    ## Remove if we don't know the accession / progeny
    dat_15_raw <- dat_15_raw[!is.na(dat_15_raw$Accession_progeny_2015),]
    
    glimpse(dat_15_raw)
    dim(dat_15_raw)
    
    
  ## Join 2017 and 2015 data
    dat_all <- left_join(dat_17_raw, dat_15_raw, 
                         by = c("Accession_2017" = "Accession_2015", 
                                "Progeny_2017" = "Progeny_2015"))
    dim(dat_all)
    
    
    
  ## Make sure 2015 and 2017 data match up - should all sum to 0
    
    sum(dat_all$Site_2015 != dat_all$Site_2017, na.rm = TRUE)
    sum(dat_all$Section_2015 != dat_all$Section_2017, na.rm = TRUE)
    sum(dat_all$Block_2015 != dat_all$Block_2017, na.rm = TRUE)
    sum(dat_all$Row_2015 != dat_all$Row_2017, na.rm = TRUE)
    sum(dat_all$Line_2015 != dat_all$Line_2017, na.rm = TRUE)
    sum(dat_all$Accession_progeny_2015 != dat_all$Accession_progeny_2017, na.rm = TRUE)
    
    
  ## Filter down and rename columns
    dat_all <- dat_all %>%
      select(contains("Problem"), Site_2017, Section_2017, Block_2017, Row_2017, Line_2017, 
             Locality_2017, Accession_2017, Progeny_2017, Accession_progeny_2017, Alive_2017,
             contains("Height"), contains("Comments")) %>%
      rename(Site = Site_2017, Section = Section_2017, Block = Block_2017, Row = Row_2017,
             Line = Line_2017, Locality = Locality_2017, Accession = Accession_2017,
             Progeny = Progeny_2017, Accession_progeny = Accession_progeny_2017,
             Height_2017 = `Height (cm)_2017`, Height_max_2017 = `Height max (cm)_2017`,
             Height_2015 = `Height (cm)_2015`) %>%
      select_all(., tolower)
    
    

# Formatting rows and lines -----------------------------------------------

    ## Make rows and columns unique to use as random effect latter
    dat_all$section_row <- paste0(dat_all$section, "-", dat_all$row)
    dat_all$section_line <- paste0(dat_all$section, "-", dat_all$line)
    
    ## Make sure they are unique
      count =  dat_all %>%
          group_by(section_row, section_line) %>%
          mutate(count = n()) 
      table(count$count) ## Should all be 1
      which(count$count > 1) ## Should return nothing
    
    

# Formatting height and RGR -----------------------------------------------

      ### Setting max height to main height if secondary stem is larger
      for(i in 1:nrow(dat_all)){
        if(!is.na(dat_all$height_max_2017[i])){
          if(dat_all$height_max_2017[i] > dat_all$height_2017[i]){
            dat_all$height_2017[i] <- dat_all$height_max_2017[i]
          }
        }
      }
      
      
      #### Round height values to the nearest 5
      dat_all$height_2015 <- round(dat_all$height_2015 / 5) * 5
      dat_all$height_2017 <- round(dat_all$height_2017 / 5) * 5
      
      ## Calculate height difference
      dat_all$height_dif <- dat_all$height_2017 - dat_all$height_2015
      
      summary(dat_all$height_dif); hist(dat_all$height_dif)
      
      plot(dat_all$height_2015, dat_all$height_dif, pch = 19, cex = 0.5)
      
      
    ## Calculate relative growth rates: over 2 years
      dat_all$rgr <- (log(dat_all$height_2017) - log(dat_all$height_2015))/ 2
      
      summary(dat_all$rgr)
      
      ## Plots of 2015 vs 2017
      plot(jitter(dat_all$height_2015), jitter(dat_all$height_2017),
           pch = 19, cex = 0.5, las = 1,
           xlab = "Height 2015 (cm)", ylab = "Height 2017 (cm)")
      abline(a = 0, b = 1, lwd = 2, col = "steelblue")
      
      ## Plot of RGR
      plot(dat_all$height_2015, dat_all$rgr, pch = 19, cex = 0.5, las = 1)
      cor.test(dat_all$height_2015, dat_all$rgr)
      
      ## Histogram of RGR
      ggplot(dat_all, aes(rgr)) + geom_histogram(fill = "steelblue", 
                                                 col = "black") + theme_bw() +
        xlab("Relative growth rate (rgr)")
      
      ## View 'outliers
      # View(dat_all[dat_all$rgr < -1 & !is.na(dat_all$rgr), ])
      
      # View(dat_all[dat_all$rgr > 1 & !is.na(dat_all$rgr), ])
      
      
  
# Filtering data ----------------------------------------------------------
      
    ## Remove inviduals without an accession or progeny number
      dat_all <- dplyr::filter(dat_all, !is.na(accession))
      
  
    # ## Filter out individuals without an estimated RGR
    # dat_all <- dplyr::filter(dat_all, !is.na(rgr))
    # 
    # ## Change to NA outlier values in relative growth rate
    # dat_all$rgr[dat_all$rgr <= quantile(dat_all$rgr, 0.01, na.rm = TRUE) |
    #               dat_all$rgr >= quantile(dat_all$rgr, 0.99, na.rm = TRUE)] <- NA
    # dat_all <- dplyr::filter(dat_all, !is.na(rgr))
    # 
  

