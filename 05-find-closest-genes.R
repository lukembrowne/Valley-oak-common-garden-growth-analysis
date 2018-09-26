## Find closest genes



# Step 1 - output top SNPs to BED

# Write file with one column for chrom and one for position and 3rd for position + 1
  write_tsv(x = data.frame(chrom = unlist(lapply(strsplit(top_snps_long$snp, "_"), 
                                                 function(x) x[[2]])),
                           pos = unlist(lapply(strsplit(top_snps_long$snp, "_"),
                                               function(x) x[[3]])),
                           pos2 = as.numeric(unlist(lapply(strsplit(top_snps_long$snp, "_"),
                                                           function(x) x[[3]]))) + 1),
            "./output/temp/top_snps.bed",
            col_names = FALSE)
 
  #### 
  bedtools sort -i ./output/temp/top_snps.bed > ./output/temp/top_snps_sorted.bed
  
  bedtools closest -a ./output/temp/top_snps_sorted.bed -b ./output/temp/genes.bed -D ref > ./output/temp/closest_genes.txt
  
  
  closest_genes <- read_tsv("./output/temp/closest_genes.txt", col_names = FALSE)
  
  head(closest_genes)
  dim(closest_genes)
  
  
  ## Read in BLAST hits to arabidopsis
  
  blasted <- read_tsv("./output/temp/genemodels.TAIRpep", col_names = FALSE)
  
  head(blasted)
  
  blasted <- blasted %>%
    dplyr::select(X1, X2, X11)
  
  colnames(blasted) <- c("qlobata", "arabidopsis", "description")  
  
  ## Join with closest genes
  
  closest_genes <- left_join(closest_genes, blasted, by = c("X7" = "qlobata"))
  
  View(closest_genes)
  
  
  closest_genes2 <- closest_genes %>%
    dplyr::filter(!is.na(description))
  
  View(closest_genes2)
  
  
  ## Testing out fisher test
  
  bedtools fisher -a ./output/temp/top_snps_sorted.bed -b ./output/temp/genes.bed -g ./output/temp/chrom_lengths.txt
  