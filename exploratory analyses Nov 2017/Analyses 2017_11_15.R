
## Load required libraries
library(tidyverse)
library(ggplot2)
#library(devtools)
#install_github("ujjwalkarn/xda")
library(xda)
library(maptools)
library(scatterpie)

library(lme4)
library(MuMIn)
library(car)
library(qtlcharts)
library(sjPlot)
library(knitr)
library(reshape2)


library(raster)
library(rasterVis)
library(ggfortify)
library(factoextra)

#devtools::install_github('royfrancis/pophelper')
library(pophelper) ## For admixture results



setwd("E:/Dropbox/Projects/2017 - Common garden/analyses/exploratory analyses Nov 2017")

setwd("~/Dropbox/Projects/2017 - Common garden/analyses/exploratory analyses Nov 2017")

## Load in california outline
cali_outline <- readShapePoly("../data/california_outline/california_outline.shp",
                              proj4string = CRS("+proj=longlat +datum=WGS84"))
lobata_range <- readShapePoly("../data/valley oak range/querloba.shp",
                              proj4string = CRS("+proj=longlat +datum=WGS84"))


### Load in garden data 
  dat_all <- read_csv("../data/cleaned data/Growth data  2016_2017 census 2017_11_15.csv")
  dat_all$lat <- NULL; dat_all$lon <- NULL; dat_all$elev <- NULL ## Clean up before merge
  
  ## Read in climate data of individuals and merge with 
  climate_ind <- read_csv("../data/cleaned data/Climate data for all individuals 2017_09_27.csv")
  
  
  dat_all <- left_join(dat_all, climate_ind, by = c("prov", "mom", "progeny"))
  

## Load in climate of common garden locations  
garden_climate <- read_csv("../data/cleaned data/CommonGardenClimate 2017_09_28.csv")

## Average across years
garden_climate <- garden_climate %>%
  group_by(site) %>%
  summarise_all(mean)



### Converting row and column numbers to XY locations

  ## Check for duplicate locations within each section - should all be 0
   sum(table(dat_all$row[dat_all$section == "block1"], dat_all$column[dat_all$section == "block1"]) > 1)
   sum(table(dat_all$row[dat_all$section == "north"], dat_all$column[dat_all$section == "north"]) > 1 )
   sum(table(dat_all$row[dat_all$section == "north_annex"], dat_all$column[dat_all$section == "north_annex"]) > 1)
   sum(table(dat_all$row[dat_all$section == "placerville"], dat_all$column[dat_all$section == "placerville"]) > 1)
  
dat_all$row_orig <- dat_all$row
dat_all$row <- as.numeric(factor(dat_all$row, 
                                 levels = rev(c( "A","B","C","D","E","F","G","H","I","J","K","L","M",
                                             "N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
                                             "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK",
                                             "AL","AM","AN","AO","AP")))) ## With proper ordering, reversed so plotting looks OK

## Adding a row of 3 tree sizes in placerville for the road

dat_all$row[dat_all$row_orig == "AD"]

dat_all$row[dat_all$row <= 13 & dat_all$section == "placerville"] <- dat_all$row[dat_all$row <= 13 & dat_all$section == "placerville"] - 3


table(dat_all$section)

test <- subset(dat_all, dat_all$section == "north")

plot(test$column, test$row, pch = 22 )
text(x = test$column, y = test$row, test$row_orig, cex = .71)
table(test$row_orig)
table(test$row)


########
# Add set values so that xy coordiates don't overlap across sections

dat_all$row[dat_all$section == "north"] <- dat_all$row[dat_all$section == "north"] + 50
dat_all$row[dat_all$section == "north_annex"] <- dat_all$row[dat_all$section == "north_annex"] + 100
dat_all$row[dat_all$section == "placerville"] <- dat_all$row[dat_all$section == "placerville"] + 250

plot(dat_all$column, dat_all$row, pch = 22)


### Replace max heights to normal heights if higher
for(i in 1:nrow(dat_all)){
  if(!is.na(dat_all$height_max_2017[i])){
    if(dat_all$height_max_2017[i] > dat_all$height_2017[i]){
      dat_all$height_2017[i] <- dat_all$height_max_2017[i]
    }
  }
  if(!is.na(dat_all$height_max_2016[i])){
    if(dat_all$height_max_2016[i] > dat_all$height_2016[i]){
      dat_all$height_2016[i] <- dat_all$height_max_2016[i]
    }
  }
}


#### Round height values to the nearest 5
  dat_all$height_2016 <- round(dat_all$height_2016 / 5) * 5
  dat_all$height_2017 <- round(dat_all$height_2017 / 5) * 5


## Calculate relative growth rates
  
  dat_all$rgr <- (log(dat_all$height_2017) - log(dat_all$height_2016))/1
  
  dat_all$rgr <- (dat_all$height_2017 - dat_all$height_2016) / dat_all$height_2016
  
  
  summary(dat_all$rgr)
  
  plot(jitter(dat_all$height_2016), jitter(dat_all$height_2017), pch = 19, cex = 0.5, las = 1)
  abline(a = 0, b = 1, lwd = 2, col = "steelblue")
  
  
  plot(dat_all$height_2016, dat_all$rgr, pch = 19, cex = 0.5)
  cor.test(dat_all$height_2016, dat_all$rgr)
  
  # View(dat_all[dat_all$rgr < -1 & !is.na(dat_all$rgr), ])
  # 
  # View(dat_all[dat_all$rgr > 1 & !is.na(dat_all$rgr), ])
  


## Exclude trees that died of things like gophers or are shaded by trees
 dat_all <- dat_all[-which(dat_all$exclude == "y"), ]

 ## Change to NA outlier values in relative growth rate
 dat_all$rgr[dat_all$rgr <= quantile(dat_all$rgr, 0.01, na.rm = TRUE) |
               dat_all$rgr >= quantile(dat_all$rgr, 0.99, na.rm = TRUE)] <- NA

 


#########
### Reading in GBS clustering and PCA results data


## Admixture results
## Read in Admixture ancestry estimates

  admix_results <- read.table("../data/GBS data/gbs451.GDP4.AN50.biallelic.QD10.filtered.final.recode.pruned_ld.3.Q")

  ## Read in family names output from PLINK/Admixture
  admix_moms <- read_delim("../data/GBS data/gbs451.GDP4.AN50.biallelic.QD10.filtered.final.recode.pruned_ld.fam", delim = " ", col_names = FALSE)
  
  admix_results <- cbind(admix_moms$X1, admix_results)
  colnames(admix_results) <- c("mom", paste("cluster",
                                            1:(ncol(admix_results)-1), sep = ""))
  admix_results$mom <- as.character(admix_results$mom)
  head(admix_results)
  
  ## Read in datasheet that has provenance names of each individual with GBS data
  ## This datasheet was made my ham
  admix_moms_provs <- read_delim("../data/GBS data/gbs451 samples provenance list.txt", delim = "\t")
  
  ## Join with admix_results
  
  admix_results <- left_join(admix_results, admix_moms_provs)
  sum(is.na(admix_results$prov)) ## Should equal to 0 - might have changed since NOV prov format gets messed up by excel and changed to date
  # View(admix_results[is.na(admix_results$prov), ])
  
  
  ### Average cluster assignment at the provenance level
  ## So that we are able to use information from samples across the provenance instead of just samples that have maternal trees in the common garden
  
  admix_results_avg = admix_results %>% 
    group_by(prov) %>%
    dplyr::select(prov, contains("cluster")) %>%
    summarise_all(mean, na.rm = TRUE)
  
    ## Rename columns
  colnames(admix_results_avg)[grep("cluster", colnames(admix_results_avg))] <- paste(colnames(admix_results_avg)[grep("cluster", colnames(admix_results_avg))], "_avg", sep = "" ) ## Rename columns
  
  # barplot(t(as.matrix(admix_results[, c(2,3,4)])), col=rainbow(3),
  #         xlab="Individual #", ylab="Ancestry", names.arg = admix_results$prov)
  # 
  
  
# ### Read in GBS codes with lat long to match up with ADMIXTURE results
#   
#   gbs_moms_info <- read_csv("../data/GBS data/GBS Trees info 2017_Master.csv")
#   gbs_moms_info$mom <- paste(gbs_moms_info$`Site Abbreviation`, gbs_moms_info$`Mother Plant Field ID`,
#                              sep = "")
#   
#   
#   ## Join with lat long info 
#   admix_results <- left_join(admix_results, gbs_moms_info)
#   
#   ## Count number of samples for each prov
#   
#   admix_results <- admix_results %>%
#     dplyr::filter(!is.na(cluster1)) %>%
#     group_by(`Site Abbreviation`) %>%
#     count() %>%
#     rename(gbs_sample_count = n) %>%
#     left_join(admix_results, .)
#   

  #View(admix_results)
  
  ## Look at individuals with GBS data that didn't match up
  
  #View(admix_results[is.na(admix_results$Latitude), ])
  
  ## Joine back with main dataset
  
  dat_all <- left_join(dat_all, admix_results_avg, by = "prov")
  

  ## Assign cluster based on majority assignment
  
  dat_all$cluster_assigned <- NA
  dat_all$cluster_assigned[dat_all$cluster1_avg > 0.50 & !is.na(dat_all$cluster1_avg)] <- "1"
  dat_all$cluster_assigned[dat_all$cluster2_avg > 0.50 & !is.na(dat_all$cluster2_avg)] <- "2"
  dat_all$cluster_assigned[dat_all$cluster3_avg > 0.50 & !is.na(dat_all$cluster3_avg)] <- "3"
  dat_all$cluster_assigned[dat_all$cluster4_avg > 0.50 & !is.na(dat_all$cluster4_avg)] <- "4"
  dat_all$cluster_assigned[dat_all$cluster5_avg > 0.50 & !is.na(dat_all$cluster5_avg)] <- "5"
  

  table(dat_all$cluster_assigned)
  
  dat_all$cluster_assigned <- factor(dat_all$cluster_assigned)
  
  boxplot(rgr ~ cluster_assigned, dat_all)
  
  ## How many seedlings not assigned?
  sum(is.na(dat_all$cluster_assigned)) # 1360
  
  
  ### Plot circle graph onto map of population assignment
  
  cali_fort <- fortify(cali_outline) ## Fortify to print california outline
  
  dat_all_prov_avg <- dat_all %>% ## Average across provenances
    mutate(cluster_assigned = as.numeric(cluster_assigned)) %>%
    group_by(prov) %>%
    summarise_all(mean, na.rm = T) %>%
    mutate(cluster_assigned = factor(cluster_assigned))
  
  ## Make plot
  ggplot() + geom_scatterpie(aes(x = lon, y = lat, group = prov, r = .2),
                             
                             data = dat_all_prov_avg,
                             cols = colnames(dplyr::select(dplyr::select(dat_all_prov_avg, 
                                                                contains("cluster")),
                                                  contains("_avg"))))+
                             coord_equal() + 
    geom_path(data = cali_fort, aes(x = long, y = lat, group = group )) +
    theme_bw()
  
  ggplot(dat_all_prov_avg[dat_all_prov_avg$cluster_assigned != "NaN",]) + 
    geom_point(aes(x = lon, y = lat, col = cluster_assigned), size = 3) +
    coord_equal() + 
    geom_path(data = cali_fort, aes(x = long, y = lat, group = group )) +
    theme_bw()
  
  

  
#########  
#### Read in PCA of SNPs
  
  gen_pca <- read_table("../data/GBS data/plink.eigenval", col_names = FALSE)
  
  gen_pca / sum(gen_pca)
  
  
  gen_pca_vals <- read_tsv("../data/GBS data/plink.eigenvec")
  
  plot(gen_pca_vals$PC1, gen_pca_vals$PC2, pch = 19)
  text(gen_pca_vals$PC1, gen_pca_vals$PC2, gen_pca_vals$FID)
  
  ## Join with provenance names and average across provenances - same as what's done with admixture
  
  pca_results <- left_join(gen_pca_vals, admix_moms_provs, by = c("FID" =  "mom"))
  sum(is.na(pca_results$prov)) ## Should equal to 0 - might have changed since NOV prov format gets messed up by excel and changed to date

  ### Average PCA Values for each provenance
  ## So that we are able to use information from samples across the provenance instead of just samples that have maternal trees in the common garden
  
  pca_results_avg = pca_results %>% 
    group_by(prov) %>%
    dplyr::select(prov, PC1:PC20) %>%
    summarise_all(mean, na.rm = TRUE) 
  colnames(pca_results_avg)[grep("PC", colnames(pca_results_avg))] <- paste(colnames(pca_results_avg)[grep("PC", colnames(pca_results_avg))], "_avg", sep = "" ) ## Rename columns
  
  ## Join back with main dataset
 dat_all <- left_join(dat_all, pca_results_avg)

 
 
## Calculate season variables

## Seedlings
dat_all <- dat_all %>%
  mutate(Tmax_seas = apply( dplyr::select(., Tmax01:Tmax12), 1, sd))

dat_all <- dat_all %>%
  mutate(Tave_seas = apply( dplyr::select(., Tave01:Tave12), 1, sd))

dat_all <- dat_all %>%
  mutate(Tmin_seas = apply( dplyr::select(., Tmin01:Tmin12), 1, sd))

## Garden sites
garden_climate <- garden_climate %>%
  mutate(Tmax_seas = apply( dplyr::select(., Tmax01:Tmax12), 1, sd))

garden_climate <- garden_climate %>%
  mutate(Tave_seas = apply( dplyr::select(., Tave01:Tave12), 1, sd))

garden_climate <- garden_climate %>%
  mutate(Tmin_seas = apply( dplyr::select(., Tmin01:Tmin12), 1, sd))


### Calculate max and min temperatature of the warmest month
dat_all <- dat_all %>%
  mutate(Tmax = apply( dplyr::select(., Tmax01:Tmax12), 1, max))

dat_all <- dat_all %>%
  mutate(Tmin = apply( dplyr::select(., Tmin01:Tmin12), 1, min))

garden_climate <- garden_climate %>%
  mutate(Tmax = apply( dplyr::select(., Tmax01:Tmax12), 1, max))

garden_climate <- garden_climate %>%
  mutate(Tmin = apply( dplyr::select(., Tmin01:Tmin12), 1, min))

## calculate difference between planting site and source provenance in climate variables

## Removing highly correlated climate variables
## Calculate VIF and remove most correlated variable, repeat

# climate_vars <- c(#"lat", "elev",
#                     "MAP", "MSP", "AHM", "SHM", "PAS", "CMD", "RH", "PPT_wt", "PPT_sp",
#                     "MAT",  "MWMT", "MCMT",
#                   "EMT", "EXT", 
#                   "DD5", "NFFD",
#                   "Tmin_wt", "Tmax_sm", 
#                   "Tmax_seas", "Tave_seas", "Tmin_seas")


## All climate variables - both averaged across year and by season, ecluding radiaton  
climate_vars = colnames(
  #  dplyr::select(dat_all, Tmax_wt:RH,  ## Seasonal and annual variables
  dplyr::select(dat_all, MAT:RH, ## Just annual variables
                #  dplyr::select(dat_all, Tmax_wt:RH_at, ## Seasonal variables
                -contains("Rad"), - contains("MAR"), ## Getting rid of solar radiation
                -contains("DD_0_sm"), -contains("PAS_sm"))) ## Do not vary at all..

### Climate vars from previous sork papers

#climate_vars = c("CMD", "MWMT", "MCMT", "DD5")

climate_vars_dif <- paste(climate_vars, "_dif", sep = "")

x = 1

for(var in climate_vars){
  
  dat_all[dat_all$site == "chico", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "chico", var]) - dat_all[dat_all$site == "chico", var]
  
  dat_all[dat_all$site == "placerville", climate_vars_dif[x]] <- as.numeric(garden_climate[garden_climate$site == "placerville", var]) - dat_all[dat_all$site == "placerville", var]
  
  x = x + 1
}


## Convert categorical variables to factors
dat_all$site <- factor(dat_all$site)
table(dat_all$site)

dat_all$prov <- factor(dat_all$prov)
table(dat_all$prov)

dat_all$mom <- factor(dat_all$mom)
table(dat_all$mom)

dat_all$section <- factor(dat_all$section)
table(dat_all$section)

dat_all$block <- factor(dat_all$block)
table(dat_all$block)

## Make sure values look OK
# 
# options(max.print = 10000)
# 
# numSummary(dat_all)
# 
# charSummary(dat_all)

## Correlations between variables

# grpd <- dat_all %>%
#   group_by(mom) %>% summarise_all(mean)
# 
# iplotCorr(dat_all[, climate_vars], reorder=FALSE)
# 
# 

##### Calculating relative height  

# 
# 
# ## Add in 'relative height' and branches - average height of plants grown in the same block
# dat_all$height_relative <- NA
# dat_all$branches_relative <- NA
# 
# for(block in levels(dat_all$block)){
#   
#   # cat("Working on block:", block, "\n")
#   
#   ## Relative to block
#   dat_all$height_relative[dat_all$block == block] <-  dat_all$height[dat_all$block == block] / mean(dat_all$height[dat_all$block == block], na.rm = TRUE)
#   
#   dat_all$branches_relative[dat_all$block == block] <- dat_all$branches[dat_all$block == block] / mean(dat_all$branches[dat_all$block == block], na.rm = TRUE)
#   
# }
# 
# # ggplot(dat_all, aes(x = block, y = height_relative)) + geom_boxplot() + theme_bw()
# # 
# #  ggplot(dat_all, aes(x = site, y = height_relative)) + geom_boxplot() + theme_bw()
# # 
# # 
# # ggplot(dat_all, aes(x = block, y = branches_relative)) + geom_boxplot()
# 
# 
# ### Setting elevation groups based on arbitrary cutoffs
# 
# low = 175; mid = 500
# 
# dat_all$cluster[dat_all$elev <= low] <- "low elev <175m"
# dat_all$cluster[dat_all$elev > low & dat_all$elev < mid] <- "mid elev 175-500m"
# dat_all$cluster[dat_all$elev >= mid] <- "high elev >500m"
# 
# 
# dat_all$cluster <- factor(dat_all$cluster, levels = c("low elev <175m",
#                                                       "mid elev 175-500m",
#                                                       "high elev >500m"))
# 
# table(dat_all$cluster)
# 
# 
# ggplot(dat_all, aes(x = cluster, y = elev)) + geom_boxplot() + geom_jitter()



# Scaling climate data prior to PCA - grouping by mom since all progeny will have the same climate

dat_all_scaled_mom <- dat_all %>%
  group_by(mom) %>%
  summarise_all(mean, na.rm = TRUE)

x = 1

climate_var_means <- NA
climate_var_sds <- NA

for(var in c(climate_vars_dif, climate_vars)){
  
  climate_var_means[x] <- mean(dat_all_scaled_mom[[var]])
  climate_var_sds[x] <- sd(dat_all_scaled_mom[[var]])
  
  dat_all_scaled_mom[var] <- (dat_all_scaled_mom[var] - climate_var_means[x]) / climate_var_sds[x]
  
  x = x+1
}

names(climate_var_means) <- c(climate_vars_dif, climate_vars)
names(climate_var_sds) <- c(climate_vars_dif, climate_vars)


## Run PCA on climate variables
clim_pca <- prcomp(dplyr::select(dat_all_scaled_mom, contains("_dif")), scale = FALSE) 

summary(clim_pca)


## PCA on home variables
clim_pca_home <- prcomp(dat_all_scaled_mom[, climate_vars], scale = FALSE) 

summary(clim_pca_home)



## Format as dataframe
clim_pca_vals <- as.data.frame(clim_pca$x)
clim_pca_home_vals <- as.data.frame(clim_pca_home$x)
colnames(clim_pca_home_vals) <- paste(colnames(clim_pca_home_vals), "_home", sep = "")

## Scree plot of PCA
fviz_screeplot(clim_pca)

## Biplot of PCA
autoplot(clim_pca, data = dat_all_scaled_mom, colour = "prov", loadings = TRUE, 
         loadings.label = TRUE) + theme_bw() + theme(legend.position="none") 

## Contributions of each axis
fviz_contrib(clim_pca, choice = "var", axes = c(1))
fviz_contrib(clim_pca, choice = "var", axes = c(2))
fviz_contrib(clim_pca, choice = "var", axes = c(3))
fviz_contrib(clim_pca, choice = "var", axes = c(4))
fviz_contrib(clim_pca, choice = "var", axes = c(5))

## Arrows plot
fviz_pca_var(clim_pca, axes = c(1,2))

dat_all_pca <- left_join(dat_all, bind_cols(dat_all_scaled_mom[, "mom"], clim_pca_vals))
# dat_all_pca <- left_join(dat_all_pca, bind_cols(dat_all_scaled_mom[, "mom"],
#                                             clim_pca_home_vals))


### SHOULD THIS BE SCALED BASED ON MOMS??
elev_sd <- sd(dat_all_pca$elev)
elev_mean <- mean(dat_all_pca$elev)

dat_all_pca$elev <- (dat_all_pca$elev - elev_mean) / elev_sd

dat_all_pca$height_2016_scaled <- (dat_all_pca$height_2016 - mean(dat_all_pca$height_2016, na.rm = T)) / sd(dat_all_pca$height_2016, na.rm = T)

### PC model
fit_nocluster = lmer(rgr  ~ PC1 + I(PC1^2) +
             + PC2 + I(PC2^2)  
           + PC3 + I(PC3^2) 
           + height_2016_scaled
           + (1 | block) # + (1 | site) 
           + (1 | prov) + (1 | mom), 
           data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])

### PC model
fit = lmer(rgr  ~ PC1*cluster_assigned + I(PC1^2)*cluster_assigned 
             + PC2*cluster_assigned + I(PC2^2)*cluster_assigned  
           + PC3*cluster_assigned + I(PC3^2) *cluster_assigned
           + cluster_assigned
           + height_2016_scaled
           + (1 | block)
           + (1 | prov) + (1 | mom), 
           data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])


fit_pca = lmer(rgr  ~ PC1*PC2_avg + I(PC1^2)*PC2_avg 
           + PC2*PC2_avg + I(PC2^2)*PC2_avg  
           + PC3*PC2_avg + I(PC3^2) *PC2_avg
           + height_2016_scaled
           + (1 | block)
           + (1 | prov) + (1 | mom), 
           data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])


anova(fit, fit_nocluster, fit_pca)

ranef(fit)

## Model summary

summary(fit)
print(r.squaredGLMM(fit))
sjp.lmer(fit, type = "fe", p.kr = FALSE)
sjp.lmer(fit, type = "eff")

sjp.int(fit)



### Scaling climate variables

for( var in climate_vars_dif){
  dat_all_pca[var] = (dat_all_pca[var] - climate_var_means[var]) / climate_var_sds[var]
}


fit = lmer(rgr  ~ CMD_dif*cluster_assigned  + I(CMD_dif^2) +
             + MWMT_dif*cluster_assigned  + I(MWMT_dif^2)  
             + MCMT_dif*cluster_assigned  + I(MCMT_dif^2) 
             + DD5_dif*cluster_assigned   + I(DD5_dif^2)
             + cluster_assigned 
             + height_2016_scaled   
             + (1 | block) # + (1 | site) 
             + (1 | prov) + (1 | mom), 
             data = dat_all_pca)


summary(fit)
print(r.squaredGLMM(fit))
sjp.lmer(fit, type = "fe", p.kr = FALSE)
sjp.lmer(fit, type = "eff")

sjp.int(fit)

#sjp.glmer(fit, type = "ma")

#sjp.lmer(fit, type = "fe.slope")

#sjp.int(fit, type = "eff", swap.pred = TRUE, show.ci = FALSE)



### Testing model selection

# options(na.action = "na.fail") 
# 
# dd = dredge(fit, subset = c(dc(MCMT_dif, I(MCMT_dif^2)), 
#                             dc(CMD_dif, I(CMD_dif^2)),
#                             dc(MWMT_dif, I(MWMT_dif^2)),
#                             dc(DD5_dif, I(DD5_dif^2))))
# 
# dd = dredge(fit, subset = c(dc(PC1, I(PC1^2)), 
#                             dc(PC2, I(PC2^2)),
#                             dc(PC3, I(PC3^2)),
#                             dc(PC4, I(PC4^2)),
#                             dc(elev, I(elev^2))))
# 
# #save(dd, file = "dd_out.Rdata")
# 
# dd
# 
# fit <- model.avg(dd, subset = delta < 6, fit = TRUE)
# #fit <- model.avg(dd, cumsum(weight) <= .95, fit = TRUE)



### Interaction plot with lowest and highest elevation

pred_dat <- expand.grid(
  elev = c(quantile(dat_all_pca$elev, .10),
           quantile(dat_all_pca$elev, .5),
           quantile(dat_all_pca$elev, .90)),
  PC1 = c(0, seq(min(dat_all_pca$PC1),
                 max(dat_all_pca$PC1),
                 length.out = 15)),
  PC2 = c(0, seq(min(dat_all_pca$PC2),
                 max(dat_all_pca$PC2),
                 length.out = 15)),
  PC3 = c(0, seq(min(dat_all_pca$PC3),
                 max(dat_all_pca$PC3),
                 length.out = 15)),
  PC4 = c(0, seq(min(dat_all_pca$PC4),
                 max(dat_all_pca$PC4),
                 length.out = 15)))

# pred_dat <- expand.grid( elev = c(quantile(dat_all_pca$elev, .10),
#                      quantile(dat_all_pca$elev, .5),
#                      quantile(dat_all_pca$elev, .90)),
#             CMD_dif = c(0, seq(min(dat_all_pca$CMD_dif),
#                                   max(dat_all_pca$CMD_dif),
#                           length.out = 15)),
#              MCMT_dif = c(0, seq(min(dat_all_pca$MCMT_dif),
#                                   max(dat_all_pca$MCMT_dif),
#                           length.out = 15)),
#              MWMT_dif = c(0, seq(min(dat_all_pca$MWMT_dif),
#                                   max(dat_all_pca$MWMT_dif),
#                           length.out = 15)),
#              DD5_dif = c(0, seq(min(dat_all_pca$DD5_dif),
#                                   max(dat_all_pca$DD5_dif),
#                           length.out = 15)))

grpd <- dat_all_pca %>%
  #    filter(site == "chico") %>%
  group_by(block, cluster, prov) %>%
  summarise_all(mean)


elev_labs = as.character(round(unique(pred_dat$elev)*elev_sd + elev_mean, 0))


pred_dat$height_relative <- predict(fit, pred_dat, re.form = NA)
pred_dat$elev <- factor(pred_dat$elev)


## Begin plot

vars <- c("PC1", "PC2", "PC3", "PC4")
# vars <- c("CMD_dif", "MCMT_dif", "MWMT_dif", "DD5_dif")

for(var in vars){
  
  pred_dat_temp <- pred_dat
  
  if(var != vars[1]){
    pred_dat_temp <- pred_dat_temp[pred_dat_temp[vars[1]] == 0, ]
  }
  
  if(var != vars[2]){
    pred_dat_temp <- pred_dat_temp[pred_dat_temp[vars[2]] == 0, ]
  }
  
  if(var != vars[3]){
    pred_dat_temp <- pred_dat_temp[pred_dat_temp[vars[3]] == 0, ]
  }
  
  if(var != vars[4]){
    pred_dat_temp <- pred_dat_temp[pred_dat_temp[vars[4]] == 0, ]
  }
  
  ## Backtransforming x axis for plot
  x_lo = min(dat_all_pca[[var]]) * climate_var_sds[var] +
    climate_var_means[var]
  
  x_hi = max(dat_all_pca[[var]]) * climate_var_sds[var] +
    climate_var_means[var]
  
  ## Locations in scaled range
  x_axis_at = c(min(dat_all_pca[[var]]),
                ((0 -  climate_var_means[var]) / climate_var_sds[var]), ## tick mark at 0
                max(dat_all_pca[[var]]))
  
  x_axis_labels = round(c(x_lo, 0, x_hi))
  
  
  
  p = ggplot(pred_dat_temp, aes(x = pred_dat_temp[[var]],
                                y = pred_dat_temp$height_relative)) + 
    #   col = "black")) +
    scale_color_discrete(name = "Elevation", labels = c(elev_labs)) +
    geom_point(data = grpd, aes(x = grpd[[var]],
                                y = grpd$height_relative,
                                fill = grpd$prov),
               pch = 21, size = 2, alpha = 0.5,
               inherit.aes=FALSE) +
    geom_line(size = 2)  + ggtitle(paste("", var)) +
    geom_vline(aes(xintercept = (0 - climate_var_means[var]) /
                     climate_var_sds[var]), lty = 2) +
    xlab(var) + ylab("Relative height") +
    scale_x_continuous(breaks = x_axis_at, labels = x_axis_labels) +
    labs(fill = "Elevation group", col = "Elevation") +
    geom_hline(yintercept = 1, lty = 2) +
    theme_bw()   + theme(legend.position="none", 
                         plot.margin = margin(1, 1, 1, 1, "cm"))
  
  print(p)
  
}




##########################
### RASTER PROCESSING ####
##########################

## Read in elevation dem and create a hillshade
dem <- raster("../data/dem/ClimateNA_DEM_cropped.tif")

slope = terrain(dem, opt='slope')
aspect = terrain(dem, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)


## Get general file names without extension
general_files <- list.files(path="../data/climate data/NA_NORM_8110_Bioclim_netCDF/",
                            pattern = "_cropped.tif$", full.names = FALSE)
general_files <- gsub(pattern = ".tif", replacement = "", general_files)

var_names <- paste(gsub(pattern = "_cropped", "", general_files), "_dif", sep = "")

general_files <- paste("/", general_files, sep = "") ## Need this for grep to be specific

## Read full files names of current climates
current_files <- list.files(path="../data/climate data/NA_NORM_8110_Bioclim_netCDF/",
                            pattern = "_cropped.tif$", full.names = TRUE)

## Read full file names of future climate
future_files <- list.files(path="../data/climate data/NA_ENSEMBLE_rcp85_2080s_Bioclim_netCDF/",
                           pattern = "_cropped.tif$", full.names = TRUE)


## Loop through climate variables and calculate difference between future and current climate

## Initialize stack of rasters  
dif_rast <- stack()
cur_rast <- stack()
x = 1 # For looping through var names in plot title
for(var in general_files){
  
  cat("\n ******* \n")
  cat("Calculating difference in:", var, "... \n")
  
  current_path <- current_files[grep(var, current_files)]
  cat("Reading in as current file:", current_path, "... \n")
  
  ## Load in current raster
  current_rast <- raster(current_path) 
  
  future_path <- future_files[grep(var, future_files)]
  cat("Reading in as future file:", future_path, "... \n")
  
  ## Load in future raster 
  future_rast <-  raster(future_path) 
  
  ### Calculate difference as future minus current
  
  difference_rast <- future_rast - current_rast
  
  ## Plot output
  # plot(difference_rast, main = paste("Difference in:", var_names[x]))
  
  ## Add to stack of rasters
  dif_rast <- raster::stack(dif_rast, difference_rast)
  cur_rast <- raster::stack(cur_rast, current_rast)
  x = x+1
}  

names(dif_rast) <- var_names
names(cur_rast) <- gsub("_dif", "", var_names)

dif_rast_raw  <- dif_rast

#dif_rast_raw <- stack(dif_rast, cur_rast)

## Cropping based on lobata range
dif_rast = dif_rast_raw
# dif_rast <- raster::mask(dif_rast_raw, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))


## Extract values from raster into dataframe
extracted_vals = extract(dif_rast, 1:ncell(dif_rast))
extracted_vals <- as.data.frame(extracted_vals)


### Scale extracted values to match with PCA

extracted_vals_scaled <- extracted_vals

for(var in climate_vars_dif){
  extracted_vals_scaled[, var] <- (extracted_vals[, var] - climate_var_means[var]) / 
    climate_var_sds[var]
}

summary(extracted_vals_scaled)

## Make histogram of predicted climate differences vs. those from the gardens
# 
#   extracted_vals_nona <- extracted_vals[!is.na(extracted_vals[, 1]), ] # Remove NAs
# 
#   dif_df <- tidyr::gather(extracted_vals_nona, key = "climate_var_dif",
#                         value = "dif")
#   dif_df$type <-  "Predicted climate"
# 
#   observed <- tidyr::gather(dplyr::select(dat_all, climate_vars_dif),
#                             key = "climate_var_dif",
#                         value = "dif")
#   observed$type = "Provenance differences"
# 
#   dif_df <- rbind(dif_df, observed)
# 
#   for(var in levels(factor(dif_df$climate_var_dif))){
# 
#     p = ggplot(dif_df[dif_df$climate_var_dif == var, ], aes(x = dif, fill = type)) +
#         geom_density(col = "black", alpha = 0.5) + ggtitle(var) + xlab(var) +
#         theme_bw() + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
# 
#     print(p)
# 
#   }
# 
# 
# ## Plots of predicted changes
# 
# ## CMD
# levelplot(dif_rast[["CMD_dif"]], contour = FALSE, margin = FALSE, main = "CMD_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 
# ## MCMT
# levelplot(dif_rast[["MCMT_dif"]], contour = FALSE, margin = FALSE, main = "MCMT_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 
# ## MWMT
# levelplot(dif_rast[["MWMT_dif"]], contour = FALSE, margin = FALSE, main = "MWMT_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 
# ## DD5
# levelplot(dif_rast[["DD5_dif"]], contour = FALSE, margin = FALSE, main = "DD5_dif",
#           par.settings = viridisTheme()) +  layer(sp.polygons(cali_outline, lwd=2))  +
#   layer(sp.polygons(lobata_range, fill = "black", alpha = 0.35))
# 


### PCA model  
### Make new rasters of each PC by calculating their loadings
predictions <- predict(clim_pca, newdata = extracted_vals_scaled)
predictions <- as.data.frame(predictions)

summary(predictions)

# pc1_rast <- dif_rast[[1]]
# 
# values(pc1_rast) <- predictions$PC1
# 
# levelplot(pc1_rast, contour = TRUE, margin = FALSE)
# 
# 


#### Climate vars model - making predictions

#  predictions <- extracted_vals_scaled[, climate_vars_dif]

## Extract and scale elevation data

predictions$elev <- extract(dem,  1:ncell(dem))
predictions$elev <- (predictions$elev - elev_mean) / elev_sd


## HEIGHT PREDICTIONS - Use PC scores to make predictions about height

height_predictions =  predict(fit, predictions, re.form = NA, type = "response")

summary(height_predictions)


height_rast <- dif_rast[[1]] ## Initialize raster

values(height_rast) <- height_predictions


## Plot of predicted relative height

hsTheme <- modifyList(GrTheme(), list(regions=list(alpha=.15)))
heightTheme <- modifyList(RdBuTheme(), list(regions=list(alpha=1)))

pixel_num = 1e5 ## Can make resolution better by making this 1e6 

# levelplot(height_rast, 
levelplot(raster::mask(height_rast, rgeos::gBuffer(lobata_range, width = 0, byid = TRUE)),
          contour = FALSE, margin = FALSE, par.settings = heightTheme, 
          maxpixels = pixel_num, colorkey = TRUE, main = "Relative height") +
  levelplot(hill,  margin = FALSE, par.settings=hsTheme, maxpixels = pixel_num) +
  layer(sp.polygons(cali_outline, lwd=2)) + layer(sp.polygons(lobata_range,
                                                              fill = "black", alpha = 0.15))


## Histogram of changes in relative height

## Across all of california

height_relative_range = unlist(raster::extract(height_rast, lobata_range))
height_relative_cali = unlist(raster::extract(height_rast, 1:ncell(height_rast)))

summary(height_relative_range)
summary(height_relative_cali)

## Growth within range
ggplot(data.frame(height_relative_range),
       aes(height_relative_range)) + 
  geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
  geom_vline(xintercept = 1, lty = 2) + ggtitle("Growth within range") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

sum(height_relative_range > 1, na.rm = TRUE) / length(na.omit(height_relative_range))

## Growth within Cali
ggplot(data.frame(height_relative_cali),
       aes(height_relative_cali)) + 
  geom_histogram(color = "black", fill = "steelblue2") + theme_bw() +
  geom_vline(xintercept = 1, lty = 2) + ggtitle("Growth across California") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

sum(height_relative_cali > 1, na.rm = TRUE) / length(na.omit(height_relative_cali))





### Calculating future temps within q lobata range

temp = raster("../data/climate data/NA_ENSEMBLE_rcp85_2080s_Bioclim_netCDF/Tave_sm_cropped.tif")

temp = raster("../data/climate data/NA_NORM_8110_Bioclim_netCDF/Tave_sm_cropped.tif")

temp = raster::mask(temp, rgeos::gBuffer(lobata_range,                                                      width = 0, byid = TRUE))

plot(temp)

summary(temp)








###############################
##### Random forest testing

library(randomForest)
library(forestFloor)

ran_for_data <- dat_all_pca %>%
                dplyr::filter(site == "chico") %>%
                dplyr::filter(!is.na(rgr), !is.na(cluster_assigned)) %>%
                dplyr::select(rgr, cluster_assigned, block, height_2016_scaled,
                              PC1_avg:PC5_avg, PC1:PC5,
                            #  climate_vars_dif
                              )


height_rf = randomForest(rgr ~ ., data = ran_for_data,
                         keep.inbag= T,
                           importance = TRUE)

height_rf

varImpPlot(height_rf)


ff = forestFloor(
  rf.fit = height_rf, # mandatory
  X = ran_for_data, # mandatory
  calc_np = FALSE, # TRUE or FALSE both works, makes no difference
  binary_reg = FALSE # takes no effect here when rfo$type="regression"
)

#plot partial functions of most important variables first
plot(ff, # forestFloor object
    # plot_seq = 1:7, # optional sequence of features to plot
     orderByImportance=TRUE # if TRUE index sequence by importance, else by X column
)



# library(ranger)
# 
# height_ranger <- ranger(rgr ~ ., data = ran_for_data, importance = "permutation")
# 
# sort(height_ranger$variable.importance, decreasing = T)
# 
# importance_pvalues(height_ranger)
# 
# height_ranger$r.squared
# 

library(gbm)
library(dismo)


gbm_data <- dat_all_pca %>%
            dplyr::filter(!is.na(rgr)) %>%
          #  dplyr::filter(!is.na(cluster_assigned)) %>%
            dplyr::select(rgr, 
                          #cluster1_avg:cluster3_avg,
                          height_2016,
                          PC1_avg:PC20_avg, 
                          #PC1:PC5, 
                          block, site,
                        #   prov, mom,
                            climate_vars_dif,
            )
names(gbm_data)[-1]


# height_gbm = gbm(rgr ~ ., data = gbm_data, n.trees = 1000, shrinkage = 0.01,
#                  interaction.depth = 5, bag.fraction = 0.5, train.fraction = 1.0,
#                  cv.folds = 5, keep.data = TRUE, distribution = "gaussian")

## Step function

response <-  which(colnames(gbm_data) == "rgr")
predictors <- which(colnames(gbm_data) %in% setdiff(names(gbm_data), c("rgr", "name")))

## Step function
height_gbm <- gbm.step(data = as.data.frame(gbm_data),
                        gbm.y = response, gbm.x = predictors,
                       tree.complexity = 5,
                       learning.rate = 0.001, bag.fraction = 0.5,
                       family = "gaussian")

summary(height_gbm)

## Simplify by dropping variables

height_gbm_simp <- gbm.simplify(height_gbm, n.drops = 5)


gbm.plot(height_gbm, n.plots = 12)

gbm.plot.fits(height_gbm, v = 1:6)

find.int <- gbm.interactions(height_gbm)

find.int$interactions

find.int$rank.list

plot(gbm_data$rgr ~ gbm_data$mom)



print(height_gbm)

summary(height_gbm)

gbm.perf(height_gbm, method = "OOB")
gbm.perf(height_gbm, method = "test")

best.iter <- gbm.perf(height_gbm, method="cv")
print(best.iter)

summary(height_gbm, best.iter)

plot(height_gbm, 15, best.iter, las = 1, continuous.resolution = 1000)

interact.gbm(height_gbm, gbm_data, i.var = c(5, 17), n.trees = best.iter)

test = partial(height_gbm, "height_2016_scaled", n.trees = best.iter, grid.resolution = 100,
        smooth = TRUE)

plotPartial(test, smooth = TRUE)


test = partial(height_gbm, c("PC5_avg", "height_2016_scaled"), n.trees = best.iter, grid.resolution = 100)
plotPartial(test)





###################
library(xgboost)
#install_github("AppliedDataSciencePartners/xgboostExplainer")
library(xgboostExplainer)
library(pdp)

previous_na_action <- options('na.action')
options(na.action='na.pass')
# Do your stuff...

boost_dat = sparse.model.matrix(rgr ~.-1, data = gbm_data, na.action = "na.pass")

dim(boost_dat)

options(na.action=previous_na_action$na.action)

xg <- xgb.cv(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear", 
             nfold = 3, nrounds = 2)

xg <- xgboost(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear",
              nrounds = 50, max_depth = 4, fill = TRUE)

xg

xg_import <- xgb.importance(feature_names = colnames(boost_dat), model = xg, )

xg_import

xgb.plot.importance(importance_matrix = xg_import)

xgb.plot.tree(model = xg)

explain <- buildExplainer(xg, trainingData = boost_dat, type = "regression")

partial(xg)

### 
# Tune an XGBoost model using 10-fold cross-validation
library(caret) # functions related to classification and regression training

boston.xgb <- train(x = data.matrix(boost_dat),
                    y = gbm_data$rgr, method = "xgbTree", metric = "Rsquared",
                    trControl = trainControl(method = "cv", number = 10),
                    tuneLength = 10)



##### USing h2o

library(h2o)


## Tutorial: https://github.com/h2oai/h2o-3/blob/master/h2o-docs/src/product/tutorials/gbm/gbmTuning.Rmd

## Initialize h2o cluster
h2o.shutdown()
h2o.init(nthreads = -1, enable_assertions = FALSE) 

## Write data to csv file first
gbm_data <- dat_all_pca %>%
  dplyr::filter(!is.na(rgr)) %>%
  dplyr::select(rgr, cluster1_avg:cluster_assigned, height_2016_scaled,
                PC1_avg:PC5_avg, 
              #  lat:RH, Tmax_seas:Tmin,
                PC1:PC5, 
                prov, mom, block, 
                climate_vars_dif)

## Change to h2o format
gbm_data_h2o <- as.h2o(x = gbm_data, key = gbm_data_h2o)

summary(gbm_data_h2o) ## Make sure factors are read correctly


#######################
### SPLIT FIT PREDICT 
#######################

  response <-  "rgr"
  predictors <- setdiff(names(gbm_data_h2o), c(response, "name"))
  
  ## Split data for training, testing, and validation
  ## May not be correct since data is not independent
  
    splits <- h2o.splitFrame(
      data = gbm_data_h2o, 
      ratios = c(0.7,0.15),   ## only need to specify 2 fractions, the 3rd is implied
      destination_frames = c("train.hex", "valid.hex", "test.hex"), seed = 1234
    )
    train <- splits[[1]]
    valid <- splits[[2]]
    test  <- splits[[3]]
    
    
    print(paste("Training data has", ncol(train), "columns and", nrow(train), "rows, valid has",
                nrow(valid), "rows, test has", nrow(test)))
    
  ### Run gbm model
    gbm <- h2o.gbm(x = predictors,
                   y = response,
                   training_frame = h2o.rbind(train, valid),
                  # training_frame =  train,
                  # validation_frame = valid,
                   nfolds = 5,
                  
                   model_id = "gbm",
                   
                   ## more trees is better if the learning rate is small enough 
                   ## here, use "more than enough" trees - we have early stopping
                   ntrees = 10000,   
                   
                   max_depth = 10,
                   
                   ## smaller learning rate is better (this is a good value for most datasets, but see below for annealing)
                   learn_rate=0.01,                                                         
                   
                   ## early stopping once the validation MSE doesn't improve by 
                   # at least 0.01% for 5 consecutive scoring events
                   stopping_rounds = 5, stopping_tolerance = 1e-4, stopping_metric = "MSE", 
                  
                   ## sample 80% of rows per tree
                   sample_rate = 0.8,                                                       
                   
                   ## sample 80% of columns per split
                   col_sample_rate = 0.8,                                                   
                   
                   ## fix a random number generator seed for reproducibility
                   seed = 1234,                                                             
                   
                   ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
                   score_tree_interval = 10 
       
                   )
    
    ###### TESTING
    
    ## Depth 10 is usually plenty of depth for most datasets, but you never know
    hyper_params = list( 
      ## restrict the search to the range of max_depth established above
      max_depth = seq(1,30,2),                                      
      
      ## search a large space of row sampling rates per tree
      sample_rate = seq(0.2,1,0.01),                                             
      
      ## search a large space of column sampling rates per split
      col_sample_rate = seq(0.2,1,0.01),                                         
      
      ## search a large space of column sampling rates per tree
      col_sample_rate_per_tree = seq(0.2,1,0.01),                                
      
      ## search a large space of how column sampling per split should change as a function of the depth of the split
      col_sample_rate_change_per_level = seq(0.9,1.1,0.01),                      
      
      ## search a large space of the number of min rows in a terminal node
      min_rows = 2^seq(0,log2(nrow(train))-1,1),                                 
      
      ## search a large space of the number of bins for split-finding for continuous and integer columns
      nbins = 2^seq(4,10,1),                                                     
      
      ## search a large space of the number of bins for split-finding for categorical columns
      nbins_cats = 2^seq(4,12,1),                                                
      
      ## search a few minimum required relative error improvement thresholds for a split to happen
      min_split_improvement = c(0,1e-8,1e-6,1e-4),                               
      
      ## try all histogram types (QuantilesGlobal and RoundRobin are good for numeric columns with outliers)
      histogram_type = c("UniformAdaptive","QuantilesGlobal","RoundRobin")       
    )
    
    search_criteria = list(
      ## Random grid search
      strategy = "RandomDiscrete",      
      
      ## limit the runtime to 60 minutes
      max_runtime_secs = 3600,         
      
      ## build no more than 100 models
      max_models = 100,                  
      
      ## random number generator seed to make sampling of parameter combinations reproducible
      seed = 1234,                        
      
      ## early stopping once the leaderboard of the top 5 models is converged to 0.1% relative difference
      stopping_rounds = 5,                
      stopping_metric = "MSE",
      stopping_tolerance = 1e-3
    )
    
    grid <- h2o.grid(
      ## hyper parameters
      hyper_params = hyper_params,
      
      ## full Cartesian hyper-parameter search
      search_criteria = search_criteria,
      
      ## which algorithm to run
      algorithm="gbm",
      
      ## identifier for the grid, to later retrieve it
      grid_id="final_grid",
      
      ## standard model parameters
      x = predictors, 
      y = response, 
      training_frame = train, 
      validation_frame = valid,
      
      ## more trees is better if the learning rate is small enough 
      ## here, use "more than enough" trees - we have early stopping
      ntrees = 10000,                                                            
      
      ## smaller learning rate is better
      ## since we have learning_rate_annealing, we can afford to start with a bigger learning rate
      learn_rate = 0.05,                                                         
      
      ## learning rate annealing: learning_rate shrinks by 1% after every tree 
      ## (use 1.00 to disable, but then lower the learning_rate)
      learn_rate_annealing = 0.99,                                               

      ## fix a random number generator seed for reproducibility
      seed = 1234,                                                             
      
      ## early stopping once the validation AUC doesn't improve by at least 0.01% for 5 consecutive scoring events
      stopping_rounds = 5,
      stopping_tolerance = 1e-4,
      stopping_metric = "MSE", 
      
      ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
      score_tree_interval = 10                                                
    )
    
    ## by default, display the grid search results sorted by increasing logloss (since this is a classification task)
    grid                                                                       
    
    ## sort the grid models by decreasing AUC
    sortedGrid <- h2o.getGrid("final_grid", sort_by="MSE", decreasing = FALSE)    
    sortedGrid
    
    gbm <- h2o.getModel(sortedGrid@model_ids[[1]])
    
    
    # Run DRF
    drf <- h2o.randomForest(x = predictors,
                            y = response,
                            training_frame = h2o.rbind(train, valid),
                            # training_frame =  train,
                            # validation_frame = valid,
                            nfolds = 5,
                            model_id = "rf",
                            ntrees            = 250,
                            max_depth         = 30)
    
    
    # 4- Score on holdout set & report
    train_rmse_gbm  <- h2o.rmse(gbm, train = TRUE)
    xval_rmse_gbm   <- h2o.rmse(gbm, xval = TRUE)
    test_perf_gbm <- h2o.performance(model = gbm, newdata = test)
    test_rmse_gbm   <- h2o.rmse(object = test_perf_gbm)
    print(paste0("GBM rmse TRAIN = ", train_rmse_gbm, ", rmse XVAL = ", xval_rmse_gbm, ", rmse TEST = ",
                 test_rmse_gbm))
    
    train_rmse_drf  <- h2o.rmse(drf, train = TRUE)
    xval_rmse_drf   <- h2o.rmse(drf, xval = TRUE)
    test_perf_drf <- h2o.performance(model = drf, newdata = test)
    test_rmse_drf   <- h2o.rmse(object = test_perf_drf)
    print(paste0("DRF rmse TRAIN = ", train_rmse_drf, ", rmse XVAL = ", xval_rmse_drf, ", rmse TEST = ",
                 test_rmse_drf))
   
  
  ## Show detailed model summary
  gbm
  
  drf
  
  
  ## Get the Mean Squared Error on the validation set
  h2o.mse(h2o.performance(gbm, newdata = valid)) ## Best is 0.009

  h2o.varimp(gbm)
  h2o.varimp_plot(gbm, ncol(train))
  
  h2o.varimp(drf)
  h2o.varimp_plot(drf, ncol(train))

  h2o.partialPlot(gbm, h2o.rbind(train, valid))
  
  h2o.saveModel(gbm, "h2o_model")
  
  
  
  
  ### Make predictions
  
  pred_gbm <- h2o.predict(object = gbm,
                          newdata = test)
  
  plot(as.data.frame(test)$rgr, as.data.frame(pred_gbm)$predict, pch = 19, las = 1,
       xlab = "Observed", ylab = "Predicted")
  abline(a = 0, b = 1, lwd = 3, col = "steelblue")
  
  
  ## Partial plots
  
  h2o.partialPlot(gbm, h2o.rbind(train, valid), c("height_2016_scaled", "DD_0_dif"))
  h2o.partialPlot(gbm, h2o.rbind(train, valid), "prov", nbins = 100)
  
  
  
  
  
  
  

#############################################
###### GRAVEYARD



## Spatial autocorrelation in residuals
# 
#     resids <- data.frame(y = dat_all_pca$row[!is.na(dat_all_pca$rgr)], 
#                          x = dat_all_pca$column[!is.na(dat_all_pca$rgr)], 
#                          resids = residuals(fit))
#     coordinates(resids) <- c("x", "y")
#     
#     bubble(resids, zcol = "resids")
#     
#     
#     ## Variogram
#     library(gstat)
#     
#     vario <- variogram(resids ~ 1, data = resids)
#     plot(vario)
#     
#     library(ncf)
#     
#     plot(spline.correlog(x = resids$x, y = resids$y, z = resids$resids, resamp = 0))
#     
#     
#     ## Moran's I
#     dists <- as.matrix(dist(cbind(resids$x, resids$y)))
#     
#     dists.inv <- 1/dists
#     diag(dists.inv) <- 0
#     
#     library(ape)
#     
#     Moran.I(resids$resids, dists.inv)



