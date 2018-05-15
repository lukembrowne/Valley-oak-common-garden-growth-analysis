
# Load libraries ----------------------------------------------------------
  library(tidyverse)
  library(googlesheets)
  library(maptools)
  library(qtlcharts)
  


# Load and quality check data ---------------------------------------------------------------

    
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
    
    
  ## Load in 2014 data which has heights before planted into garden
    dat_14_raw <- gs_read(gs_key("10Yss2ciOo0cazkT9Udu55QUr37T7zXeQbnGuLpTnAMM"), ws = 2)
    colnames(dat_14_raw) <- paste0(colnames(dat_14_raw), "_2014")
    dat_14_raw$Progeny_2014 <- as.character(dat_14_raw$Progeny_2014)
    
    # Remove duplicates
    dat_14_raw <- dat_14_raw[ -which(dat_14_raw$Problem_2014 == "tree number duplicated"), ]
    
    dim(dat_14_raw)
    
    # Remove NAs
    dat_14_raw <- dat_14_raw[!is.na(dat_14_raw$Accession_2014), ]
    dat_14_raw <- dat_14_raw[!is.na(dat_14_raw$Progeny_2014), ]
    
    dim(dat_14_raw)
    
  ## Join 2017 and 2015 data
    dat_all <- left_join(dat_17_raw, dat_15_raw, 
                         by = c("Accession_2017" = "Accession_2015", 
                                "Progeny_2017" = "Progeny_2015"))
    dim(dat_all)
    
  ## Add in 2014 height
    dat_all <- left_join(dat_all, 
                         dplyr::select(dat_14_raw, Accession_2014, Progeny_2014, 
                                       `Height (cm)_2014`), 
                         by = c("Accession_2017" = "Accession_2014", 
                                "Progeny_2017" = "Progeny_2014"))
    
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
      dplyr::select(contains("Problem"), Site_2017, Section_2017, Block_2017, Row_2017, Line_2017, 
             Locality_2017, Accession_2017, Progeny_2017, Accession_progeny_2017, Alive_2017,
             contains("Height"), contains("Comments")) %>%
      rename(Site = Site_2017, Section = Section_2017, Block = Block_2017, Row = Row_2017,
             Line = Line_2017, Locality = Locality_2017, Accession = Accession_2017,
             Progeny = Progeny_2017, Accession_progeny = Accession_progeny_2017,
             Height_2017 = `Height (cm)_2017`, Height_max_2017 = `Height max (cm)_2017`,
             Height_2015 = `Height (cm)_2015`, Height_2014 = `Height (cm)_2014`) %>%
      select_all(., tolower)
    

# Formatting block --------------------------------------------------------

  # Which individuals have missing block information?
    dat_all %>%
      filter(!is.na(accession)) %>%
      filter(is.na(block))  
    
  dat_all$block[dat_all$section == "North" & dat_all$row == "A"] <- 3
  dat_all$block[dat_all$section == "North" & dat_all$row == "B"] <- 3
  dat_all$block[dat_all$section == "North" & dat_all$row == "U"] <- 4
  dat_all$block[dat_all$section == "North" & dat_all$row == "W"] <- 4
  dat_all$block[dat_all$section == "North" & dat_all$row == "Y"] <- 4
  
  # Set so that all rows have the same block in North annex
  dat_all$block[dat_all$section == "North annex" ] <- 2
 
  #
  dat_all$section_block <- paste0(dat_all$section, "_", dat_all$block)
  table(dat_all$section_block)
    

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
        ggplot(dat_all, aes(rgr)) + geom_histogram(fill = "steelblue2",
                                                   col = "black") + theme_bw() +
          xlab("Relative growth rate (rgr)")
      
      ## View 'outliers
      # View(dat_all[dat_all$rgr < -1 & !is.na(dat_all$rgr), ])
      
      # View(dat_all[dat_all$rgr > 1 & !is.na(dat_all$rgr), ])
      
      

# Reading in GBS data -----------------------------------------------------

     
  ## Reading in 012 genotype matrix created with vcftools, processing and cleaning in R
  ## Should be 15,627 Loci based on file gbs451.GDP4.AN50.biallelic.QD10.MAF10.MIS10.vcf from Sorel received March 2018    
      
  ## Read in position of SNPs - first column is chromosome number, second column is position
  snp_pos <- readr::read_tsv("./data/GBS_data/cleaned_data/gbs.all.012.pos", 
                             col_names = c("chrom", "pos"))
  snp_pos
  
  ## Read in genotype cols
  
  ## provide column names since they are not included in vcftools output
  ## First column is sequence of numbers 0-450, which we will remove later
  col_names <- c("NUM", paste0(snp_pos$chrom, "_", snp_pos$pos))
  
  ## Set missing values (-1) to NA
  genotypes <- readr::read_tsv("./data/GBS_data/cleaned_data/gbs.all.012",
                               col_names = col_names, na = "-1")
  
  ## Remove first column
  genotypes <- dplyr::select(genotypes, -NUM)
  dim(genotypes) # Changes depends on filtering options
  genotypes
  
  ## Read in individual
  indv <- readr::read_tsv("./data/GBS_data/cleaned_data/gbs.all.012.indv",
                          col_names = "ID")
  dim(indv)
  indv 
  
  ## Read in file that has gbs_name, locality, and sork lab sample number, in the other of GBS samples
    sample_key <- readr::read_tsv("./data/GBS_data/gbs451.locality_sample_key.txt")
    sample_key
    
    ## Link up with accession number
    accession_key <- gs_read(gs_key("1DUEV-pqV28D6qJl6vJM6S1VLWbxc9E71afufulDRbIc"), ws = 2)
    accession_key
    
    sample_key <- dplyr::left_join(sample_key, accession_key, 
                                   by = c("locality" = "Locality", "sample" = "Sample #")) %>%
                  rename(accession = Accession)
    
    sample_key
    
    ## How many GBS samples have accession numbers? 
    length(table(sample_key$accession))
    
      ## Gbs samples NOT in common garden to check for typos etc that may be wrongly excluding samples
      # View(sample_key[is.na(sample_key$accession), ])
    
    ## Output list of samples in common garden with GBS data
      #write_tsv(data.frame(gbs_name = sample_key$gbs_name[!is.na(sample_key$accession)]),
      #  path = "./data/GBS_data/samples_w_accession_number.txt",
      #  col_names = FALSE)
    
      
    ## Join individual names and genotype data and subset to those that are in the common garden

    gen_dat <- dplyr::left_join(dplyr::select(sample_key, gbs_name, accession),
                                dplyr::bind_cols(indv, genotypes),
                                                 by = c("gbs_name" = "ID")) 
    
    dim(gen_dat)
    gen_dat
    
  ## Subset gen_dat to just individuals we have GBS data for
    gen_dat <- gen_dat %>%
      filter(gbs_name %in% indv$ID)
  
      
      
  
# Filtering data ----------------------------------------------------------
      
    ## Remove inviduals without an accession or progeny number
      dat_all <- dplyr::filter(dat_all, !is.na(accession))
    
    ## Remove individuals with mechanical damage
      mech <- dplyr::filter(dat_all, grepl('mechanical', comments_2015) |
                              grepl('mechanical', comments_2017))
      mech$comments_2017
      dat_all <- dplyr::filter(dat_all, !(accession_progeny %in% mech$accession_progeny))

    ## Filtering based on height
      
      if(response_variable == "height"){
      
        ## Filter out NA heights
        dat_all <- dplyr::filter(dat_all, !is.na(height_2017))
        dat_all <- dplyr::filter(dat_all, !is.na(height_2014))  
  
        # dat_all$height_2017[dat_all$height_2017 <= quantile(dat_all$height_2017, 0.01,
        #                                                     na.rm = TRUE) |
        #               dat_all$height_2017 >= quantile(dat_all$height_2017, 0.99,
        #                                               na.rm = TRUE)] <- NA
      
      }

    
    ## Filtering based on RGR

      if(response_variable == "rgr"){

    ## Individuals with negative growth rates
      # neg_rgr <- dplyr::filter(dat_all, rgr < quantile(dat_all$rgr, 0.025, na.rm = TRUE)) %>% arrange(rgr)
      #View(dplyr::select(neg_rgr, rgr, height_2015, height_2017) )
      
    ## Individuals with high positive growth rates
      # pos_rgr <- dplyr::filter(dat_all, rgr > quantile(dat_all$rgr, 0.99, na.rm = TRUE)) %>% arrange(rgr)
      # View(dplyr::select(pos_rgr, accession_progeny, rgr, height_2015, height_2017) )
      
      # Change to NA outlier values in relative growth rate
      dat_all$rgr[dat_all$rgr <= quantile(dat_all$rgr, 0.025, na.rm = TRUE) |
                    dat_all$rgr >= quantile(dat_all$rgr, 0.975, na.rm = TRUE)] <- NA
  
  
      ## Filter out individuals without an estimated RGR
       dat_all <- dplyr::filter(dat_all, !is.na(rgr))
      }
      
      
  ## Filtering out individuals without genetic data
      ## This happens again in 03 adding climate script to filter out moms who don't have climate data!
      dat_all <- dplyr::filter(dat_all, accession %in% gen_dat$accession)    
      
      
    ## Filtering based on survival  
      

    dim(dat_all)
  

