
# Finding the closest genes ----------------------------------------------


## Read in qlobata BLAST hits to arabidopsis
  blasted <- read_tsv("./output/finding_closest_genes/genemodels.TAIRpep", col_names = FALSE)
  
  dim(blasted)
  head(blasted)
  
  blasted <- blasted %>%
    dplyr::select(X1, X2, X11)
  
  colnames(blasted) <- c("qlob_gene_name", "athal_gene_name", "athal_description")
  head(blasted)
  
  # Remove duplicate genes
  blasted <- blasted %>%
    dplyr::distinct(qlob_gene_name, athal_gene_name, .keep_all = TRUE)
  
  ## Remove .ti from qlobata gene names
  blasted$qlob_gene_name <-unlist(lapply(strsplit(blasted$qlob_gene_name, 
                                                  split = "\\."), function(x) x[[1]]))
  
  blasted$athal_gene_name <- unlist(lapply(strsplit(blasted$athal_gene_name, 
                                                    split = "\\."), function(x) x[[1]]))
  
  head(blasted)
  
  ## Extract out gene name from description
  blasted$athal_gene_description <- unlist(lapply(strsplit(x = blasted$athal_description,
                                                           split = " \\| "),
                                                           function(x) x[[3]]))
  head(blasted)
  dim(blasted)

  
## Read in arabidopsis go terms
  ara_go <- read_tsv("./data/GBS_data/arabidopsis_go/ATH_GO_GOSLIM.txt",
                     col_names = c("locus_name", "TAIR_accession", "object_name",
                                   "relationship_type", "GO_term", "GO_ID",
                                   "TAIR_keyword_ID", "Aspect", "GOslim",
                                   "Evidence_code", "Evidence_with", "Reference",
                                   "Annotator", "Date_annotated"),
                     guess_max = 50000)
  
  head(ara_go)
  dim(ara_go)

# Subset down to just columns that are useful
  ara_go <- ara_go %>%
    dplyr::select(locus_name, relationship_type, GO_term, GOslim, 
                  GO_ID, Aspect, Evidence_with)
  
  head(ara_go)
  
# Subset down to just biological processes
  ara_go <- ara_go %>%
    dplyr::filter(Aspect == "P")
  
  dim(ara_go)
  

  
  
# Step 1 - output outlier genotypes to BED

# Write file with one column for chrom and one for position and 3rd for position + 1

# Top snps
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(unique(top_snps_long$snp), "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(unique(top_snps_long$snp), "_"),
                                               function(x) x[[3]])),
                           pos2 = unlist(lapply(strsplit(unique(top_snps_long$snp), "_"),
                                                function(x) x[[3]]))),
            "./output/finding_closest_genes/top_snps.bed",
            col_names = FALSE)

  
# Bottom snps  
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(unique(bottom_snps_long$snp), "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(unique(bottom_snps_long$snp), "_"),
                                               function(x) x[[3]])),
                           pos2 =unlist(lapply(strsplit(unique(bottom_snps_long$snp), "_"),
                                               function(x) x[[3]]))),
            "./output/finding_closest_genes/bottom_snps.bed",
            col_names = FALSE)
  
  
# All snps  
  write_tsv(x = tibble(chrom = snp_pos$chrom,
                           pos = snp_pos$pos,
                           pos2 = snp_pos$pos),
            "./output/finding_closest_genes/all_snps.bed",
            col_names = FALSE)  
  

 
#### Use bedtools to find closest genes - run in terminal
  bedtools sort -i ./output/finding_closest_genes/top_snps.bed > ./output/finding_closest_genes/top_snps_sorted.bed
  bedtools closest -a ./output/finding_closest_genes/top_snps_sorted.bed -b ./output/finding_closest_genes/qlobatagenes.bed -D ref > ./output/finding_closest_genes/top_closest_genes.txt
  
  bedtools sort -i ./output/finding_closest_genes/bottom_snps.bed > ./output/finding_closest_genes/bottom_snps_sorted.bed
  bedtools closest -a ./output/finding_closest_genes/bottom_snps_sorted.bed -b ./output/finding_closest_genes/qlobatagenes.bed -D ref > ./output/finding_closest_genes/bottom_closest_genes.txt
  
  bedtools sort -i ./output/finding_closest_genes/all_snps.bed > ./output/finding_closest_genes/all_snps_sorted.bed
  bedtools closest -a ./output/finding_closest_genes/all_snps_sorted.bed -b ./output/finding_closest_genes/qlobatagenes.bed -D ref > ./output/finding_closest_genes/all_closest_genes.txt
  
 
# Read back in as dataframes  
  col_names <- c("chrom_a", "pos1_a", "pos2_a", "chrom_b", "pos1_b", "pos2_b",
                 "qlob_gene_name",
                 "unk2", "unk3", "species", "gene_type", "unk4", "qlob_info", "distance")
 
  # Top genes
    top_closest_genes_raw <- read_tsv("./output/finding_closest_genes/top_closest_genes.txt", 
                                  col_names = col_names)
    head(top_closest_genes_raw)
    dim(top_closest_genes_raw)
    
    # Filter down and join with blast hits and go terms
    top_closest_genes <- top_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted) %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name")) %>%
      dplyr::mutate(outlier_type = "Warm advantageous")
    
    # Add in genotype, q value, and height prediction
     top_closest_genes <-  top_closest_genes %>%
        mutate(snp = paste("snp", chrom_a, pos1_a, sep = "_")) %>%
        dplyr::left_join(., top_snps_long, by = "snp")
     
     table(top_closest_genes$genotype)
     summary(top_closest_genes$q_val)
     summary(top_closest_genes$height_change_warmer_base0)
    
    head(top_closest_genes)
    dim(top_closest_genes)
    
    top_snps_long$snp %in% paste0("snp_", top_closest_genes$chrom_a, "_",
                                  top_closest_genes$pos1_a) # Should be all true
  
 # bottom genes
    bottom_closest_genes_raw <- read_tsv("./output/finding_closest_genes/bottom_closest_genes.txt", 
                                      col_names = col_names)
    head(bottom_closest_genes_raw)
    dim(bottom_closest_genes_raw)
    
    # Filter down
    bottom_closest_genes <- bottom_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted)  %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)  %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name")) %>%
      dplyr::mutate(outlier_type = "Warm disadvantageous")
    
    # Add in genotype, q value, and height prediction
    bottom_closest_genes <-  bottom_closest_genes %>%
      mutate(snp = paste("snp", chrom_a, pos1_a, sep = "_")) %>%
      dplyr::left_join(., bottom_snps_long, by = "snp")
    
    table(bottom_closest_genes$genotype)
    summary(bottom_closest_genes$q_val)
    summary(bottom_closest_genes$height_change_warmer_base0)
    
    head(bottom_closest_genes)
    dim(bottom_closest_genes)
    
    all(bottom_snps_long$snp %in% paste0("snp_", bottom_closest_genes$chrom_a, "_", bottom_closest_genes$pos1_a)) # Should be all true
    
    
 # All genes
    all_closest_genes_raw <- read_tsv("./output/finding_closest_genes/all_closest_genes.txt", 
                                         col_names = col_names)
    head(all_closest_genes_raw)
    dim(all_closest_genes_raw)
    
    # Filter down
    all_closest_genes <- all_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted)  %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)  %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name"))  %>%
      dplyr::mutate(outlier_type = "All")
    
    dim(all_closest_genes)
    

    
 ## Output to table
    
    outlier_genes_raw <- bind_rows(top_closest_genes, bottom_closest_genes)
    
  # Select columns we need  
    outlier_genes <- outlier_genes_raw %>%
      dplyr::select(outlier_type, chrom_a, pos1_a, genotype, q_val, height_change_warmer_base0, qlob_gene_name, distance,
                     athal_gene_name, athal_gene_description, GOslim, GO_term,
                    GO_ID)
    
    head(outlier_genes)
    
 outlier_genes_formatted = outlier_genes %>%
            dplyr::distinct(outlier_type, qlob_gene_name, 
                            athal_gene_name, athal_gene_description, GO_term, 
                            .keep_all = TRUE) %>%
            dplyr::group_by(outlier_type, chrom_a, pos1_a, genotype, q_val, height_change_warmer_base0, qlob_gene_name,
                            distance,  athal_gene_description, 
                          athal_gene_name) %>%
            dplyr::summarise(GO_terms = toString(GO_term)) %>%
            dplyr::ungroup() %>% 
            dplyr::rename(`Outlier type` = outlier_type,
                         Chromosome = chrom_a,
                         Position = pos1_a,
                         Genotype = genotype,
                         `Q value` = q_val,
                         `Change in RGR` = height_change_warmer_base0,
                         `Distance to gene (bp)` = distance,
                         `Q. lobata gene name` = qlob_gene_name,
                         `A. thaliana gene name` = athal_gene_name,
                         `Gene description` = athal_gene_description,
                         `GO terms` = GO_terms)
 
 # Add information about climate outliers
 outlier_genes_formatted <- outlier_genes_formatted %>%
   mutate(`LGM Climate outlier` = FALSE) %>%
   dplyr::select(`Outlier type`, Chromosome, Position, 
                 Genotype, `Q value`, `Change in RGR`, `LGM Climate outlier`,
                 everything())
 
outlier_genes_formatted$`LGM Climate outlier`[paste("snp", outlier_genes_formatted$Chromosome, 
      outlier_genes_formatted$Position, sep = "_") %in% double_outliers] <- TRUE

table(outlier_genes_formatted$`LGM Climate outlier`)
 
 head(outlier_genes_formatted)
 dim(outlier_genes_formatted)
 
 
 # write_tsv(outlier_genes_formatted,
 #           paste0("./figs_tables/Table S3 - formatted gene list with go terms ", Sys.Date(), ".txt"))
 
 
 # Quantify distances to genes
 
 # How many located within a gene?
 sum(outlier_genes_formatted$`Distance to gene (bp)` == 0) / 
   length(outlier_genes_formatted$`Distance to gene (bp)`)
 
 # Average distance to gene?
 summary(abs(outlier_genes_formatted$`Distance to gene (bp)`))
 
 
 ## Gene distances of all SNPs
 
all_genes_sub = all_closest_genes %>%
   dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)

summary(abs(all_genes_sub$distance))

sum(all_genes_sub$distance == 0) / length(all_genes_sub$distance)

wilcox.test(abs(outlier_genes_formatted$`Distance to gene (bp)`),
            abs(all_genes_sub$distance))

table(outlier_genes_formatted$`Distance to gene (bp)` <= 5000)
table(all_genes_sub$distance <= 5000)

fisher.test(rbind(table(outlier_genes_formatted$`Distance to gene (bp)` <= 5000), 
      table(all_genes_sub$distance <= 5000) ))

 
 ## Barplot of GOslim terms
 
 ## Add in total go terms across all blasted genes
 qlob_genes <- left_join(blasted, ara_go, 
                         by = c("athal_gene_name" = "locus_name")) %>%
               dplyr::distinct(qlob_gene_name, GOslim, .keep_all = TRUE) %>%
               dplyr::mutate(outlier_type = "Overall gene set")
 dim(qlob_genes)
 
 outlier_genes_barchart <- outlier_genes_raw %>%
   dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
   dplyr::select(outlier_type, GOslim) %>%
   dplyr::bind_rows(., dplyr::select(qlob_genes, outlier_type, GOslim)) 
 
 ## Filter down to terms that are in both beneficial and detrimental genotypes
 outlier_GOslims <- unique(c(outlier_genes_barchart$GOslim[outlier_genes_barchart$outlier_type == "Warm advantageous"], 
        outlier_genes_barchart$GOslim[outlier_genes_barchart$outlier_type == "Warm disadvantageous"]))
 
 outlier_genes_barchart <- outlier_genes_barchart %>%
   dplyr::filter(GOslim %in% outlier_GOslims)
 
 length(table(outlier_genes_barchart$GOslim))
 
 
 table(outlier_genes_barchart$outlier_type)
 
 outlier_genes_barchart$GOslim[is.na(outlier_genes_barchart$GOslim)] <- "unclassified"

 ggplot(outlier_genes_barchart, aes(x = GOslim, y = ..prop.., 
                           group = outlier_type, fill = outlier_type)) + 
   geom_bar(position = "dodge") +
   coord_flip() +
   xlab("GO term") +
   ylab("Percentage") +
   scale_fill_manual(name = "", values = c("grey70", "#6699CC", "#FF6633")) + 
   theme_bw(10) + 
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
 # Save to file
 # ggsave(filename = paste0("./figs_tables/Figure S8 - goslim terms barchart ",
 #                          Sys.Date(), ".pdf"),
 #        units = "cm",
 #        height = 12, width = 20,
 #        useDingbats = FALSE )


 
 #  Output gene names in arabidopsis format to use for enrichment tests on pantherdb
outlier_genes_raw %>%
   dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
   dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
   dplyr::select(outlier_type, athal_gene_name) %>%
   dplyr::filter(outlier_type == "Warm advantageous") %>%
   dplyr::select(athal_gene_name) %>%
 write_tsv(. ,
           "./output/finding_closest_genes/top_closest_genes_arath_names.txt", col_names = FALSE)

outlier_genes_raw %>%
  dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
  dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
  dplyr::select(outlier_type, athal_gene_name) %>%
  dplyr::filter(outlier_type == "Warm disadvantageous") %>%
  dplyr::select(athal_gene_name) %>%
  write_tsv(. ,
            "./output/finding_closest_genes/bottom_closest_genes_arath_names.txt", col_names = FALSE)


outlier_genes_raw %>%
  dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
  dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
  dplyr::select(outlier_type, athal_gene_name) %>%
  dplyr::select(athal_gene_name) %>%
  write_tsv(. ,
            "./output/finding_closest_genes/top_and_bottom_closest_genes_arath_names.txt", col_names = FALSE)

# Write list of genes within 5kb of snps
all_closest_genes %>%
  dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
  dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
  dplyr::select(athal_gene_name) %>%
  write_tsv(. ,
            "./output/finding_closest_genes/all_snps_arath_names.txt", col_names = FALSE)



 
 
 
 

# Lagniappe ---------------------------------------------------------------

 ## Make reference list based on qlobata genes
 arath_genes <- read_tsv("./data/GBS_data/genemodels.TAIRpep", 
                         col_names = c("qlob_gene_name", "arath_gene_name",
                                       "unk1", "unk2", "unk3", "unk4", "unk5", "unk6",
                                       "unk7", "unk8","unk9"))
 
 # Remove suffix
 stripped_gene_names <-  unlist(lapply(strsplit(arath_genes$arath_gene_name, 
                                                split = "\\."), function(x) x[[1]]))
 length(stripped_gene_names)
 
 # Check for duplicates
 sort(table(stripped_gene_names), decreasing = T)[1:10]
 stripped_gene_names <- unique(stripped_gene_names)

 write_tsv(as.data.frame(stripped_gene_names), 
           "./output/finding_closest_genes/arath_genes_reference.txt", col_names = FALSE) 
 
 
 
 
 
 # Look at distances between genes    
 bind_rows(top_closest_genes, mid_closest_genes, bottom_closest_genes) %>%
   ggplot(., aes(distance, fill = outlier_type)) + 
   geom_density(alpha = 0.5) +
   theme_bw(15)
 
    
  