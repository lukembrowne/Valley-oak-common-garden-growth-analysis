
## Load packages
library(vcfR)

## Read in GBS data

setwd("E:/Dropbox/Projects/2017 - Common garden/analyses/R scripts")


vcf <- read.vcfR("../data/GBS data/gbs451.GDP4.AN50.biallelic.QD10.filtered.recode.vcf")

head(vcf)


strwrap(vcf@meta)
queryMETA(vcf) ## Find 

queryMETA(vcf, element = "AD") ## Find definitions of symbols

getFIX(vcf) ## Get first columns

vcf@gt[1:6, 1:4] ## Get genotypes

head(vcf)

chrom <- create.chromR(vcf)

plot(chrom)

chromoqc(chrom)

### 

    tbl=read.table("../data/GBS data/gbs451.GDP4.AN50.biallelic.QD10.filtered.recode.pruned_ld.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(3),
          xlab="Individual #", ylab="Ancestry", border=NA)

