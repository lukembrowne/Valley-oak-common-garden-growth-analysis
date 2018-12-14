


# Finding the closest genes ----------------------------------------------


## Read in qlobata BLAST hits to arabidopsis
  blasted <- read_tsv("./output/temp/genemodels.TAIRpep", col_names = FALSE)
  
  dim(blasted)
  head(blasted)
  
  blasted <- blasted %>%
    dplyr::select(X1, X2, X11)
  
  colnames(blasted) <- c("qlob_gene_name", "athal_gene_name", "athal_description")
  
  # Remove duplicate genes
  blasted <- blasted %>%
    dplyr::distinct(qlob_gene_name, athal_gene_name, .keep_all = TRUE)
  
  ## Remove .ti from qlobata gene names
  blasted$qlob_gene_name <-unlist(lapply(strsplit(blasted$qlob_gene_name, 
                                                  split = "\\."), function(x) x[[1]]))
  
  blasted$athal_gene_name <- unlist(lapply(strsplit(blasted$athal_gene_name, 
                                                    split = "\\."), function(x) x[[1]]))
  
  head(blasted)
  
  ## Extract out gene name
  blasted$athal_gene_description <- unlist(lapply(strsplit(x = blasted$athal_description,
                                                           split = " \\| "),
                                                           function(x) x[[3]]))
  
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
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(top_snps_long$snp, "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(top_snps_long$snp, "_"),
                                               function(x) x[[3]])),
                           pos2 = as.numeric(unlist(lapply(strsplit(top_snps_long$snp, "_"),
                                                           function(x) x[[3]])))),
            "./output/temp/top_snps.bed",
            col_names = FALSE)
  
# Mid snps  
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(mid_snps_long$snp, "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(mid_snps_long$snp, "_"),
                                               function(x) x[[3]])),
                           pos2 = as.numeric(unlist(lapply(strsplit(mid_snps_long$snp, "_"),
                                                           function(x) x[[3]])))),
            "./output/temp/mid_snps.bed",
            col_names = FALSE)
  
# Bottom snps  
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(bottom_snps_long$snp, "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(bottom_snps_long$snp, "_"),
                                               function(x) x[[3]])),
                           pos2 = as.numeric(unlist(lapply(strsplit(bottom_snps_long$snp, "_"),
                                                           function(x) x[[3]])))),
            "./output/temp/bottom_snps.bed",
            col_names = FALSE)
  
  
# All snps  
  write_tsv(x = tibble(chrom = snp_pos$chrom,
                           pos = snp_pos$pos,
                           pos2 = snp_pos$pos),
            "./output/temp/all_snps.bed",
            col_names = FALSE)  
  

 
#### Use bedtools to find closest genes - run in terminal
  bedtools sort -i ./output/temp/top_snps.bed > ./output/temp/top_snps_sorted.bed
  bedtools closest -a ./output/temp/top_snps_sorted.bed -b ./output/temp/genes.bed -D ref > ./output/temp/top_closest_genes.txt

  bedtools sort -i ./output/temp/mid_snps.bed > ./output/temp/mid_snps_sorted.bed
  bedtools closest -a ./output/temp/mid_snps_sorted.bed -b ./output/temp/genes.bed -D ref > ./output/temp/mid_closest_genes.txt
  
  bedtools sort -i ./output/temp/bottom_snps.bed > ./output/temp/bottom_snps_sorted.bed
  bedtools closest -a ./output/temp/bottom_snps_sorted.bed -b ./output/temp/genes.bed -D ref > ./output/temp/bottom_closest_genes.txt
  
  bedtools sort -i ./output/temp/all_snps.bed > ./output/temp/all_snps_sorted.bed
  bedtools closest -a ./output/temp/all_snps_sorted.bed -b ./output/temp/genes.bed -D ref > ./output/temp/all_closest_genes.txt
  
 
# Read back in as dataframes  
  col_names <- c("chrom_a", "pos1_a", "pos2_a", "chrom_b", "pos1_b", "pos2_b",
                 "qlob_gene_name",
                 "unk2", "unk3", "species", "gene_type", "unk4", "qlob_info", "distance")
 
  # Top genes
    top_closest_genes_raw <- read_tsv("./output/temp/top_closest_genes.txt", 
                                  col_names = col_names)
    head(top_closest_genes_raw)
    dim(top_closest_genes_raw)
    
    # Filter down and join with blast hits and go terms
    top_closest_genes <- top_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted) %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name")) %>%
      dplyr::mutate(outlier_type = "Beneficial genotype")
    
    dim(top_closest_genes)
    
    top_snps_long$snp %in% paste0("snp_", top_closest_genes$chrom_a, "_",
                                  top_closest_genes$pos1_a) # Should be all true
  
# Mid genes
    mid_closest_genes_raw <- read_tsv("./output/temp/mid_closest_genes.txt", 
                                      col_names = col_names)
    head(mid_closest_genes_raw)
    dim(mid_closest_genes_raw)
    
    # Filter down
    mid_closest_genes <- mid_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted)  %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)  %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name")) %>%
      dplyr::mutate(outlier_type = "mid")
    
    dim(mid_closest_genes)
    
    all(mid_snps_long$snp %in% paste0("snp_", mid_closest_genes$chrom_a, "_", mid_closest_genes$pos1_a)) # Should be all true
    which(!mid_snps_long$snp %in% paste0("snp_", mid_closest_genes$chrom_a, "_", mid_closest_genes$pos1_a))
    
    mid_snps_long$snp[13]
   
 # bottom genes
    bottom_closest_genes_raw <- read_tsv("./output/temp/bottom_closest_genes.txt", 
                                      col_names = col_names)
    head(bottom_closest_genes_raw)
    dim(bottom_closest_genes_raw)
    
    # Filter down
    bottom_closest_genes <- bottom_closest_genes_raw %>%
      dplyr::filter(gene_type == "gene") %>%
      dplyr::left_join(., blasted)  %>%
      dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)  %>%
      dplyr::left_join(., ara_go, by = c("athal_gene_name" = "locus_name")) %>%
      dplyr::mutate(outlier_type = "Detrimental genotype")
    
    dim(bottom_closest_genes)
    
    all(bottom_snps_long$snp %in% paste0("snp_", bottom_closest_genes$chrom_a, "_", bottom_closest_genes$pos1_a)) # Should be all true
    
    
 # All genes
    all_closest_genes_raw <- read_tsv("./output/temp/all_closest_genes.txt", 
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
    
    all_genes_raw <- bind_rows(top_closest_genes, bottom_closest_genes)
    
  # Select columns we need  
    all_genes <- all_genes_raw %>%
      dplyr::select(outlier_type, chrom_a, pos1_a, qlob_gene_name,
                    distance, athal_gene_name, athal_gene_description, GOslim, GO_term,
                    GO_ID)
    
    head(all_genes)
    

 all_genes_formatted =  all_genes %>%
            dplyr::distinct(outlier_type, qlob_gene_name, 
                            athal_gene_name, athal_gene_description, GO_term, 
                            .keep_all = TRUE) %>%
            dplyr::group_by(outlier_type, chrom_a, pos1_a, qlob_gene_name,
                            athal_gene_description, distance, 
                          athal_gene_name) %>%
            dplyr::summarise(GO_terms = toString(GO_term)) %>%
            dplyr::ungroup() %>% 
            dplyr::rename(`Outlier type` = outlier_type,
                         Chromosome = chrom_a,
                         Position = pos1_a,
                         `Q. lobata gene name` = qlob_gene_name,
                         `Distance to gene` = distance,
                         `A. thaliana gene name` = athal_gene_name,
                         `Gene description` = athal_gene_description,
                         `GO terms` = GO_terms)
 
 write_tsv(all_genes_formatted, 
           paste0("./figs_tables/Table S3 - formatted gene list with go terms ", Sys.Date(), ".txt"))
 
 
 ## Quantify distances to genes
 
 # How many located within a gene?
 sum(all_genes_formatted$`Distance to gene` == 0) / length(all_genes_formatted$`Distance to gene`)
 
 # Average distance to gene?
 summary(abs(all_genes_formatted$`Distance to gene`))
 
 
 ## Gene distances of all SNPs
 
all_genes_sub = all_closest_genes %>%
   dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE)

summary(abs(all_genes_sub$distance))

sum(all_genes_sub$distance == 0) / length(all_genes_sub$distance)

wilcox.test(abs(all_genes_formatted$`Distance to gene`),
            abs(all_genes_sub$distance))
            

 
 ## Barplot of GOslim terms
 
 ## Add in total go terms across all blasted genes
 
 qlob_genes <- left_join(blasted, ara_go, 
                         by = c("athal_gene_name" = "locus_name")) %>%
               dplyr::distinct(qlob_gene_name, GOslim, .keep_all = TRUE) %>%
               dplyr::mutate(outlier_type = "Overall gene set")
 dim(qlob_genes)
 
 all_genes_barchart <- all_genes_raw %>%
   dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
   dplyr::select(outlier_type, GO_term) %>%
   dplyr::bind_rows(., dplyr::select(qlob_genes, outlier_type, GO_term)) 
 
 ## Filter down to terms that are in both beneficial and detrimental genotypes
 outlier_go_terms <- unique(c(all_genes_barchart$GO_term[all_genes_barchart$outlier_type == "Beneficial genotype"], 
        all_genes_barchart$GO_term[all_genes_barchart$outlier_type == "Detrimental genotype"]))
 
 all_genes_barchart <- all_genes_barchart %>%
   dplyr::filter(GO_term %in% outlier_go_terms)
 
 length(table(all_genes_barchart$GO_term))
 
 
 table(all_genes_barchart$outlier_type)
 
 all_genes_barchart$GO_term[is.na(all_genes_barchart$GO_term)] <- "unclassified"

 ggplot(all_genes_barchart, aes(x = GO_term, y = ..prop.., 
                           group = outlier_type, fill = outlier_type)) + 
   geom_bar(position = "dodge") +
   coord_flip() +
   xlab("GO term") +
   ylab("Percentage") +
   scale_fill_manual(name = "", values = c("#6699CC", "#FF6633", "grey70")) + 
   theme_bw(10) + 
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
 # Save to file
 ggsave(filename = paste0("./figs_tables/Figure S6 - goslim terms barchart ",
                          Sys.Date(), ".pdf"),
        units = "cm",
        height = 10, width = 30,
        useDingbats = FALSE )

 
 
 #  
# all_genes_raw %>%
#    dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
#    dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
#    dplyr::select(outlier_type, athal_gene_name) %>%
#    dplyr::filter(outlier_type == "Beneficial genotype") %>%
#    dplyr::select(athal_gene_name) %>%
#  write_tsv(. , 
#            "./output/temp/top_closest_genes_arath_names.txt", col_names = FALSE) 
#   
# all_genes_raw %>%
#   dplyr::filter(distance <= 5000 & distance >= -5000) %>% # Restrict to just 5kb
#   dplyr::distinct(chrom_a, pos1_a, .keep_all = TRUE) %>%
#   dplyr::select(outlier_type, athal_gene_name) %>%
#   dplyr::filter(outlier_type == "Beneficial genotype") %>%
#   dplyr::select(athal_gene_name) %>%
#   write_tsv(. , 
#             "./output/temp/bottom_closest_genes_arath_names.txt", col_names = FALSE) 
# 
#  
 
 
 
 
 
 
 
 

# Lagniappe ---------------------------------------------------------------

 ## Make reference list based on qlobata genes
 arath_genes <- read_tsv("./data/GBS_data/genemodels.TAIRpep", 
                         col_names = c("qlob_gene_name", "arath_gene_name",
                                       "unk1", "unk2", "unk3", "unk4", "unk5", "unk6",
                                       "unk7", "unk8","unk9"))
 
 # Remove suffix
 stripped_gene_names <-  unlist(lapply(strsplit(arath_genes$arath_gene_name, 
                                                split = "\\."), function(x) x[[1]]))
 
 
 write_tsv(as.data.frame(unique(stripped_gene_names)), 
           "./output/temp/arath_genes_reference.txt", col_names = FALSE) 
 
 
 
 
 
 # Look at distances between genes    
 bind_rows(top_closest_genes, mid_closest_genes, bottom_closest_genes) %>%
   ggplot(., aes(distance, fill = outlier_type)) + 
   geom_density(alpha = 0.5) +
   theme_bw(15)
 
    
  