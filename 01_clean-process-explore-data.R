
# Load libraries ----------------------------------------------------------
  library(tidyverse)
  library(googlesheets)
  library(maptools)
  
  library(gamm4)
  library(visreg)
  library(rrBLUP)
  library(qvalue) # devtools::install_github("jdstorey/qvalue")
  library(MuMIn)
  library(vegan)
  library(beepr)
  library(patchwork) # devtools::install_github("thomasp85/patchwork")
  library(rasterVis)
  library(pdp)
  library(ggrepel)
  library(caret)
  library(flashpcaR) # devtools::install_github("gabraham/flashpca/flashpcaR")

  library(psych)
  library(factoextra)

  library(foreach)
  library(doParallel)


  # Function to transform from raw to scaled variables
  forward_transform <- function(x, var, means, sds){
    ( x - means[var]) / sds[var]
  }
  
  # Function to back transform variables
  back_transform <- function(x, var, means, sds){
    (x * sds[var]) + means[var]
  }


# Choose model type -------------------------------------------------------

    response_variable = "rgr"

# Load and quality check data ---------------------------------------------------------------
    
  ## Load in garden data from 2017 - "Qlobata census measurements 2017"

    dat_17_raw <- read_csv("./data/cleaned_data/Qlobata census measurements 2017 - Measurement data 2017.csv", guess_max = 50000)
    colnames(dat_17_raw) <- paste0(colnames(dat_17_raw), "_2017")
    
    ## Remove if we don't know the accession / progeny
    dat_17_raw <- dat_17_raw[!is.na(dat_17_raw$Accession_progeny_2017),]
    
    glimpse(dat_17_raw) 
    dim(dat_17_raw)
    
  ## Load in garden data from 2015 - "Qlobata census measurements 2015"
    dat_15_raw <- read_csv("./data/cleaned_data/Qlobata census measurements 2015 - Measurement Data 2015.csv", guess_max = 50000)
    colnames(dat_15_raw) <- paste0(colnames(dat_15_raw), "_2015")
    
    ## Remove if we don't know the accession / progeny
    dat_15_raw <- dat_15_raw[!is.na(dat_15_raw$Accession_progeny_2015),]
    
    glimpse(dat_15_raw)
    dim(dat_15_raw)
    
    
  ## Load in 2014 data which has heights before planted into garden
    dat_14_raw <- read_csv("./data/cleaned_data/Qlobata census measurements 2014 - Measurement data 2014.csv", guess_max = 50000)
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
             Locality_2017, Accession_2017, Progeny_2017, Accession_progeny_2017, Alive_2017, `Border tree_2017`,
             contains("Height"), contains("Comments")) %>%
      rename(Site = Site_2017, Section = Section_2017, Block = Block_2017, Row = Row_2017,
             Line = Line_2017, Locality = Locality_2017, Accession = Accession_2017,
             Progeny = Progeny_2017, Accession_progeny = Accession_progeny_2017, 
             Border_tree = `Border tree_2017`,
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
  
    ## Add coordinates for rows and lines
      #  dat_all$row <- factor(dat_all$row, levels = c(LETTERS,
      #                                               "AA", "AB", "AC", "AD", "AE",
      #                                               "AF", "AG", "AH", "AI", "AJ",
      #                                               "AK", "AL", "AM", "AN", "AO",
      #                                               "AP"))
      # table(dat_all$row)
      # dat_all$row <- as.numeric(dat_all$row)
      # 
      # 
      # # Separate sections
      # dat_all$row[dat_all$section == "Block1"] <- dat_all$row[dat_all$section == "Block1"] + 100
      # dat_all$line[dat_all$section == "Block1"] <- dat_all$line[dat_all$section == "Block1"] + 100
      # 
      # dat_all$row[dat_all$section == "North"] <- dat_all$row[dat_all$section == "North"] + 200
      # dat_all$line[dat_all$section == "North"] <- dat_all$line[dat_all$section == "North"] + 200
      # 
      # dat_all$row[dat_all$section == "North annex"] <- dat_all$row[dat_all$section == "North annex"] + 300
      # dat_all$line[dat_all$section == "North annex"] <- dat_all$line[dat_all$section == "North annex"] + 300
      
      # plot(dat_all$row, dat_all$line, type = "n", xlim = c(0, 50), ylim = c(0, 75))
      # text(dat_all$row, dat_all$line, 
      #      labels = as.character(dat_all$locality), cex = 0.25)
      # 
      # 
      # 

# Formatting height and RGR -----------------------------------------------

      ### Setting max height to main height if secondary stem is larger
      for(i in 1:nrow(dat_all)){
        if(!is.na(dat_all$height_max_2017[i])){
          if(dat_all$height_max_2017[i] > dat_all$height_2017[i]){
            dat_all$height_2017[i] <- dat_all$height_max_2017[i]
          }
        }
      }
      
    ## Calculate relative growth rates: over 3 years
      dat_all$rgr <- (log(dat_all$height_2017) - log(dat_all$height_2014)) / 3
      

      summary(dat_all$rgr)
      

# Reading in GBS data -----------------------------------------------------

     
  ## Reading in 012 genotype matrix created with vcftools, processing and cleaning in R
      
  ## MOMS IN GARDEN      
    ## Using iterative filtering
    snp_pos <- readr::read_tsv("./data/GBS_data/iterative filtering/cleaned_data/gbs451_FIL-4_noLD.012.pos", col_names = c("chrom", "pos"))
      
    snp_pos
    
    ## Read in genotype cols
    
    ## provide column names since they are not included in vcftools output
    ## First column is sequence of numbers 0-450, which we will remove later
    col_names <- c("NUM", paste0("snp_", snp_pos$chrom, "_", snp_pos$pos))

    ## Set missing values (-1) to NA

    # Iterative filtering
    genotypes <- readr::read_tsv("./data/GBS_data/iterative filtering/cleaned_data/gbs451_FIL-4_noLD.012",
                                 col_names = col_names, na = "-1")
    
    ## Remove first column
    genotypes <- dplyr::select(genotypes, -NUM)
    dim(genotypes) # Changes depends on filtering options
    genotypes
    
    ## Read in individual

    # Iterative filtering
    indv <- readr::read_tsv("./data/GBS_data/iterative filtering/cleaned_data/gbs451_FIL-4_noLD.012.indv",
                            col_names = "ID")
    dim(indv)
    indv 
    
    ## Read in file that has gbs_name, locality, and sork lab sample number, in the other of GBS samples
      sample_key <- readr::read_tsv("./data/GBS_data/gbs451.locality_sample_key.txt")
      sample_key
      
      ## Link up with accession number
      accession_key <- read_csv("./data/cleaned_data/Maternal tree locations and sample IDs - Locations and sample IDs.csv")
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
      
      dim(gen_dat)
      
      
      
  ## ALL TREES WITH GBS SAMPLES  
      
      ## Set missing values (-1) to NA
      # Iterative filtering
      genotypes_all <- readr::read_tsv("./data/GBS_data/iterative filtering/cleaned_data/gbs451_FIL-4_noLD.012",
                                       col_names = col_names, na = "-1")
      
      ## Remove first column
      genotypes_all <- dplyr::select(genotypes_all, -NUM)
      dim(genotypes_all) # Changes depends on filtering options
      genotypes_all
      
      ## Read in individual
      # Iterative filtering
      indv_all <- readr::read_tsv("./data/GBS_data/iterative filtering/cleaned_data/gbs451_FIL-4_noLD.012.indv",
                                  col_names = "ID")
      dim(indv_all)
      indv_all
    
      ## Join individual names and genotype data and subset to those that are in the common garden
      gen_dat_all <- dplyr::left_join(dplyr::select(sample_key, gbs_name, accession),
                                  dplyr::bind_cols(indv_all, genotypes_all),
                                  by = c("gbs_name" = "ID")) 
      
      dim(gen_dat_all)
      gen_dat_all

    
# Filtering data ----------------------------------------------------------
      
    ## Remove inviduals without an accession or progeny number
      dat_all <- dplyr::filter(dat_all, !is.na(accession))
    
    ## Remove individuals with mechanical damage
      mech <- dplyr::filter(dat_all, grepl('mechanical', comments_2015) |
                              grepl('mechanical', comments_2017))
      mech$comments_2017
      mech$comments_2015
      dat_all <- dplyr::filter(dat_all, !(accession_progeny %in% mech$accession_progeny))
      
      # Weed whacker
      mech <- dplyr::filter(dat_all, grepl('whacker', comments_2015) |
                              grepl('whacker', comments_2017))
      mech$comments_2015
      mech$comments_2017
      dat_all <- dplyr::filter(dat_all, !(accession_progeny %in% mech$accession_progeny))
      
      # Gopher
      mech <- dplyr::filter(dat_all, grepl('gopher', comments_2015) |
                              grepl('gopher', comments_2017))
      mech$comments_2017
      mech$comments_2015
      dat_all <- dplyr::filter(dat_all, !(accession_progeny %in% mech$accession_progeny))

    
    ## Filtering based on RGR
      if(response_variable == "rgr"){

          ## Filter out individuals without an estimated RGR
           dat_all <- dplyr::filter(dat_all, !is.na(rgr))
    
          # Exclude individuals with a negative growth rate
          dat_all <- dplyr::filter(dat_all, rgr >= 0)

      }
      
    ## Remove border trees  
    dat_all <- dplyr::filter(dat_all, border_tree == 0)
    
    ## Remove individuals without a locality
    dat_all <- dplyr::filter(dat_all, !is.na(locality))
      

    ## Filtering based on garden
     # dat_all <- dat_all[dat_all$site == "Chico", ]
    
    dim(dat_all)
     

