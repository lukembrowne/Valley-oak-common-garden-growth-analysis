
# Source files ------------------------------------------------------------

source("./01_clean-process-explore-data.R")
source("./03_climate-pca.R") ## Adds PCA climate dif variables


# Load libraries ----------------------------------------------------------

library(gamm4)
library(sjPlot)


# GAMM models -------------------------------------------------------------

# Null model with only random effects
fit_null<- gamm4(rgr  ~ 1,
                 random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                 data = dat_all_scaled)

# Height in 2016 only
fit_height_only <- gamm4(rgr  ~ s(height_2016),
                         random = ~ (1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                         data = dat_all_scaled)

## Individual climate variables
fit_CMD <- gamm4(rgr ~ s(height_2016)
                         + s(CMD_dif),
                        random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                        data = dat_all_scaled)

fit_Tmax <- gamm4(rgr ~ s(height_2016)
                 + s(Tmax_dif),
                 random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                 data = dat_all_scaled)

fit_Tmin <- gamm4(rgr ~ s(height_2016)
                  + s(Tmin_dif),
                  random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                  data = dat_all_scaled)

fit_DD5 <- gamm4(rgr ~ s(height_2016)
                  + s(DD5_dif),
                  random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                  data = dat_all_scaled)

# Checking model ----------------------------------------------------------

check_model <- function(gamm4_out){
  
  plot(gamm4_out$gam, residuals = FALSE, shade = TRUE, 
       scale = 0, seWithMean=FALSE, pages = 0)
  
  print(summary(gamm4_out$gam))
  
  print(summary(gamm4_out$mer))
  
  print(anova(gamm4_out$gam))
  
  #gam.check(gamm4_out$gam)

}


plot_gam <- function(mod, var, var_index, xlab, plot_raw_pts = TRUE, PCA = FALSE){
  
  plt <- plot(mod, pages = 1) ## Extract info from mgcv plot.gam function
  
  ## Format into DF
  plt_df <- data.frame(x_scaled = plt[[var_index]]$x,
                       fit = plt[[var_index]]$fit, ### Doesn't account for intercept!!
                       se = plt[[var_index]]$se)
  
  ## Scale x axis back to original units
  if(PCA == FALSE) {
   plt_df$x_orig <- (plt_df$x_scaled * scaled_var_sds[[var]]) + scaled_var_means[[var]]
  } else {
    plt_df$x_orig <- plt_df$x_scaled
  }
  
  ## Add back in intercept to fitted values
  plt_df$fit <- plt_df$fit + mod$coefficients[["(Intercept)"]]
  
  ## GGplot
 p =  ggplot(plt_df, aes(x = x_orig, y = fit)) + geom_line(lwd = 1.5) +
    geom_ribbon(aes(ymin = fit - 2*se, ymax = fit + 2*se), alpha = 0.3) +
    geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2) + 
    ylab("Relative growth rate") + xlab(xlab) + 
    theme_bw() + theme(axis.text=element_text(size=12),
                       axis.title=element_text(size=14),
                       plot.margin = unit(c(1,1,1,1), "cm"))
 
 if(plot_raw_pts == TRUE){
  p = p + geom_point(data = dat_all, aes(x = dat_all[[var]], y = rgr), 
                     alpha = 0.2, cex = 0.75)  
 }
 
 print(p)
  
}



  ## Height only
  check_model(fit_height_only)
  plot_gam(fit_height_only$gam, var = "height_2016", var_index = 1, 
           xlab = "Height 2016", plot_raw_pts = T)
  ggsave(filename = "./figs_tables/Height only.pdf", scale = .75)
  
  ## CMD
  check_model(fit_CMD)
  plot_gam(fit_CMD$gam, var = "CMD_dif", var_index = 2, 
           xlab = "CMD_dif", plot_raw_pts = F)
  ggsave(filename = "./figs_tables/CMD_dif no pts.pdf", scale = .75)
  
  ## Tmax
  check_model(fit_Tmax)
  plot_gam(fit_Tmax$gam, var = "Tmax_dif", var_index = 2, 
           xlab = "Tmax_dif", plot_raw_pts = F)
  ggsave(filename = "./figs_tables/Tmax_dif no pts.pdf", scale = .75)
  
  ## Tmin
  check_model(fit_Tmin)
  plot_gam(fit_Tmin$gam, var = "Tmin_dif", var_index = 2,
           xlab = "Tmin_dif", plot_raw_pts = F)
  ggsave(filename = "./figs_tables/Tmin_dif pts.pdf", scale = .75)
  
  
  ## DD5
  check_model(fit_DD5)
  plot_gam(fit_DD5$gam, var = "DD5_dif", var_index = 2,
           xlab = "DD5_dif", plot_raw_pts = F)
  ggsave(filename = "./figs_tables/DD5_dif no pts.pdf", scale = .75)
  
  
  

# Climate models with genetic interactions --------------------------------

  ### Need to rerun 01_clean process script with filtering for GBS data
  
  ## Individual climate variables
  fit_CMD_clust <- gamm4(rgr ~ s(height_2016)
                   + s(CMD_dif, by = cluster_assigned),
                   random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                   data = dat_all_scaled)
  
  fit_Tmin_clust <- gamm4(rgr ~ s(height_2016)
                    + s(Tmin_dif, by = cluster_assigned),
                    random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                    data = dat_all_scaled)
  
  fit_Tmax_clust <- gamm4(rgr ~ s(height_2016)
                          + s(Tmax_dif, by = cluster_assigned)
                          random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                          data = dat_all_scaled)
  
  fit_DD5_clust <- gamm4(rgr ~ s(height_2016)
                   + s(DD5_dif, by = cluster_assigned),
                   random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                   data = dat_all_scaled)
  

  ## Compare models with and without genetic cluster interactions
    anova(fit_CMD_clust$mer, fit_CMD$mer)
    anova(fit_Tmin_clust$mer, fit_Tmin$mer)
    anova(fit_Tmax_clust$mer, fit_Tmax$mer)
    anova(fit_DD5_clust$mer, fit_DD5$mer)
    
    
    fit_Tmax_clust <- gamm4(rgr ~ s(height_2016)
                            + s(Tmax_dif, by = Tmax),
                            random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                            data = dat_all_scaled)
    
    
    

# Adding elevation to models ----------------------------------------------

    ## Individual climate variables
    fit_CMD_elev <- gamm4(rgr ~ s(height_2016)
                           +s(CMD_dif)
                           +s(elev)
                           + t2(CMD_dif, elev),
                           random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                           data = dat_all_scaled)
    
    fit_Tmin_elev <- gamm4(rgr ~ s(height_2016)
                            + t2(Tmin_dif, elev),
                            random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                            data = dat_all_scaled)
    
    fit_Tmax_elev <- gamm4(rgr ~ s(height_2016)
                            +s(Tmax_dif)
                            +s(elev),
                          #  + t2(Tmax_dif, elev),
                            random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                            data = dat_all_scaled)
    
    ## DD5
    fit_DD5_elev <- gamm4(rgr ~ s(height_2016)
                           +s(DD5_dif)
                           +s(elev)
                           + t2(DD5_dif, elev),
                           random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                           data = dat_all_scaled)
    
    fit_DD5_elev_noint <- gamm4(rgr ~ s(height_2016)
                          +s(DD5_dif)
                          +s(elev),
                          random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                          data = dat_all_scaled)
    
    
    ## Compare models with and without genetic elever interactions
    anova(fit_CMD_elev$mer, fit_CMD$mer) ## No interaction
    anova(fit_Tmin_elev$mer, fit_Tmin$mer) ## 
    anova(fit_Tmax_elev$mer, fit_Tmax$mer) ## Elev but no interaction
    anova(fit_DD5_elev$mer, fit_DD5_elev_noint$mer, fit_DD5$mer) ##  Elev but no interaction
    
    
    anova(fit_Tmax$mer, fit_Tmax_old$mer)
    
    vis.gam(fit_Tmax_elev_int$gam, view = c("elev", "Tmax_dif"), theta = 45, phi = 25)
    
    vis.gam(fit_Tmax_elev_int$gam, view = c("elev", "Tmax_dif"), plot.type = "contour")
    
    ## Tmax
    check_model(fit_Tmax_elev)
    plot_gam(fit_Tmax_elev$gam, var = "elev", var_index = 3, 
             xlab = "elev", plot_raw_pts = F)
    
  

# Models with PC of climate variables -------------------------------------

    # Null model with only random effects
    fit_null<- gamm4(rgr  ~ 1,
                     random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                     data = dat_all_scaled)
    
    # Height in 2016 only
    fit_height_only <- gamm4(rgr  ~ s(height_2016),
                             random = ~ (1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                             data = dat_all_scaled)
    
    ## PC Climate variables
    fit_clim_pca <- gamm4(rgr ~ s(height_2016)
                     + PC1_clim_dif
                     + s(PC2_clim_dif)
                     + s(PC3_clim_dif),
                     random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                     data = dat_all_scaled)
    
    fit_clim_pca_clust <- gamm4(rgr ~ s(height_2016)
                          + s(PC1_clim_dif, by = cluster_assigned)
                          + s(PC2_clim_dif, by = cluster_assigned)
                          + s(PC3_clim_dif, by = cluster_assigned)
                          + s(PC4_clim_dif, by = cluster_assigned),
                          random = ~(1|block) + (1|prov) + (1|mom) + (1|row) + (1|column),
                          data = dat_all_scaled)
    
    
    
    ## Check PC climate variables
    check_model(fit_clim_pca)
    
    ## These don't work right if not a smoothing factor!!
    
    plot_gam(fit_clim_pca$gam, var = "PC1_clim_dif", var_index = 2, 
             xlab = "PC1_clim_dif", plot_raw_pts = F, PCA = TRUE)
    
    plot_gam(fit_clim_pca$gam, var = "PC2_clim_dif", var_index = 3, 
             xlab = "PC2_clim_dif", plot_raw_pts = F, PCA = TRUE)
    
    plot_gam(fit_clim_pca$gam, var = "PC3_clim_dif", var_index = 4, 
             xlab = "PC3_clim_dif", plot_raw_pts = F, PCA = TRUE)
    
    plot_gam(fit_clim_pca$gam, var = "PC4_clim_dif", var_index = 5, 
             xlab = "PC4_clim_dif", plot_raw_pts = F, PCA = TRUE)
    

    
      

# Graveyard ---------------------------------------------------------------



# Spatial autocorrelation -------------------------------------------------



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



# GLMM --------------------------------------------------------------------

  # ### PC model
  # fit_nocluster = lmer(rgr  ~ PC1 + I(PC1^2) +
  #                        + PC2 + I(PC2^2)  
  #                      + PC3 + I(PC3^2) 
  #                      + height_2016_scaled
  #                      + (1 | block) # + (1 | site) 
  #                      + (1 | prov) + (1 | mom), 
  #                      data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # ### PC model
  # fit = lmer(rgr  ~ PC1*cluster_assigned + I(PC1^2)*cluster_assigned 
  #            + PC2*cluster_assigned + I(PC2^2)*cluster_assigned  
  #            + PC3*cluster_assigned + I(PC3^2) *cluster_assigned
  #            + cluster_assigned
  #            + height_2016_scaled
  #            + (1 | block)
  #            + (1 | prov) + (1 | mom), 
  #            data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # 
  # fit_pca = lmer(rgr  ~ PC1*PC2_avg + I(PC1^2)*PC2_avg 
  #                + PC2*PC2_avg + I(PC2^2)*PC2_avg  
  #                + PC3*PC2_avg + I(PC3^2) *PC2_avg
  #                + height_2016_scaled
  #                + (1 | block)
  #                + (1 | prov) + (1 | mom), 
  #                data = dat_all_pca[!is.na(dat_all_pca$cluster_assigned),])
  # 
  # 
  # anova(fit, fit_nocluster, fit_pca)
  # 
  # ranef(fit)
  # 
  # ## Model summary
  # 
  # summary(fit)
  # print(r.squaredGLMM(fit))
  # sjp.lmer(fit, type = "fe", p.kr = FALSE)
  # sjp.lmer(fit, type = "eff")
  # 
  # sjp.int(fit)



# Random forest models ----------------------------------------------------

  # library(randomForest)
  # library(forestFloor)
  
  # ran_for_data <- dat_all_pca %>%
  #   dplyr::filter(site == "chico") %>%
  #   dplyr::filter(!is.na(rgr), !is.na(cluster_assigned)) %>%
  #   dplyr::select(rgr, cluster_assigned, block, height_2016_scaled,
  #                 PC1_avg:PC5_avg, PC1:PC5
  #                 #  climate_vars_dif
  #   )
  # 
  # 
  # height_rf = randomForest(rgr ~ ., data = ran_for_data,
  #                          keep.inbag= T,
  #                          importance = TRUE)
  # 
  # height_rf
  # 
  # varImpPlot(height_rf)
  # 
  # 
  # ff = forestFloor(
  #   rf.fit = height_rf, # mandatory
  #   X = ran_for_data, # mandatory
  #   calc_np = FALSE, # TRUE or FALSE both works, makes no difference
  #   binary_reg = FALSE # takes no effect here when rfo$type="regression"
  # )
  # 
  # #plot partial functions of most important variables first
  # plot(ff, # forestFloor object
  #      # plot_seq = 1:7, # optional sequence of features to plot
  #      orderByImportance=TRUE # if TRUE index sequence by importance, else by X column
  # )
  # 
  
  
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



# Gradient boosted machines -----------------------------------------------

  # library(gbm)
  # library(dismo)
  
  
  # gbm_data <- dat_all_pca %>%
  #   dplyr::filter(!is.na(rgr)) %>%
  #   #  dplyr::filter(!is.na(cluster_assigned)) %>%
  #   dplyr::select(rgr, 
  #                 #cluster1_avg:cluster3_avg,
  #                 height_2016,
  #                 PC1_avg:PC20_avg, 
  #                 #PC1:PC5, 
  #                 block, site,
  #                 #   prov, mom,
  #                 climate_vars_dif,
  #   )
  # names(gbm_data)[-1]
  # 
  # 
  # # height_gbm = gbm(rgr ~ ., data = gbm_data, n.trees = 1000, shrinkage = 0.01,
  # #                  interaction.depth = 5, bag.fraction = 0.5, train.fraction = 1.0,
  # #                  cv.folds = 5, keep.data = TRUE, distribution = "gaussian")
  # 
  # ## Step function
  # 
  # response <-  which(colnames(gbm_data) == "rgr")
  # predictors <- which(colnames(gbm_data) %in% setdiff(names(gbm_data), c("rgr", "name")))
  # 
  # ## Step function
  # height_gbm <- gbm.step(data = as.data.frame(gbm_data),
  #                        gbm.y = response, gbm.x = predictors,
  #                        tree.complexity = 5,
  #                        learning.rate = 0.001, bag.fraction = 0.5,
  #                        family = "gaussian")
  # 
  # summary(height_gbm)
  # 
  # ## Simplify by dropping variables
  # 
  # height_gbm_simp <- gbm.simplify(height_gbm, n.drops = 5)
  # 
  # 
  # gbm.plot(height_gbm, n.plots = 12)
  # 
  # gbm.plot.fits(height_gbm, v = 1:6)
  # 
  # find.int <- gbm.interactions(height_gbm)
  # 
  # find.int$interactions
  # 
  # find.int$rank.list
  # 
  # plot(gbm_data$rgr ~ gbm_data$mom)
  # 
  # 
  # 
  # print(height_gbm)
  # 
  # summary(height_gbm)
  # 
  # gbm.perf(height_gbm, method = "OOB")
  # gbm.perf(height_gbm, method = "test")
  # 
  # best.iter <- gbm.perf(height_gbm, method="cv")
  # print(best.iter)
  # 
  # summary(height_gbm, best.iter)
  # 
  # plot(height_gbm, 15, best.iter, las = 1, continuous.resolution = 1000)
  # 
  # interact.gbm(height_gbm, gbm_data, i.var = c(5, 17), n.trees = best.iter)
  # 
  # test = partial(height_gbm, "height_2016_scaled", n.trees = best.iter, grid.resolution = 100,
  #                smooth = TRUE)
  # 
  # plotPartial(test, smooth = TRUE)
  # 
  # 
  # test = partial(height_gbm, c("PC5_avg", "height_2016_scaled"), n.trees = best.iter, grid.resolution = 100)
  # plotPartial(test)
  


  # 
  # # library(xgboost)
  # # library(pdp)
  # 
  # previous_na_action <- options('na.action')
  # options(na.action='na.pass')
  # # Do your stuff...
  # 
  # boost_dat = sparse.model.matrix(rgr ~.-1, data = gbm_data, na.action = "na.pass")
  # 
  # dim(boost_dat)
  # 
  # options(na.action=previous_na_action$na.action)
  # 
  # xg <- xgb.cv(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear", 
  #              nfold = 3, nrounds = 2)
  # 
  # xg <- xgboost(data = boost_dat, label = gbm_data$rgr, objective = "reg:linear",
  #               nrounds = 50, max_depth = 4, fill = TRUE)
  # 
  # xg
  # 
  # xg_import <- xgb.importance(feature_names = colnames(boost_dat), model = xg, )
  # 
  # xg_import
  # 
  # xgb.plot.importance(importance_matrix = xg_import)
  # 
  # xgb.plot.tree(model = xg)
  # 
  # explain <- buildExplainer(xg, trainingData = boost_dat, type = "regression")
  # 
  # partial(xg)
  # 
  # ### 
  # # Tune an XGBoost model using 10-fold cross-validation
  # #library(caret) # functions related to classification and regression training
  # 
  # boston.xgb <- train(x = data.matrix(boost_dat),
  #                     y = gbm_data$rgr, method = "xgbTree", metric = "Rsquared",
  #                     trControl = trainControl(method = "cv", number = 10),
  #                     tuneLength = 10)



# Using h2o for RF and GBM ------------------------------------------------

  # 
  # #library(h2o)
  # 
  # 
  # ## Tutorial: https://github.com/h2oai/h2o-3/blob/master/h2o-docs/src/product/tutorials/gbm/gbmTuning.Rmd
  # 
  # ## Initialize h2o cluster
  # h2o.shutdown()
  # h2o.init(nthreads = -1, enable_assertions = FALSE) 
  # 
  # ## Write data to csv file first
  # gbm_data <- dat_all_pca %>%
  #   dplyr::filter(!is.na(rgr)) %>%
  #   dplyr::select(rgr, cluster1_avg:cluster_assigned, height_2016_scaled,
  #                 PC1_avg:PC5_avg, 
  #                 #  lat:RH, Tmax_seas:Tmin,
  #                 PC1:PC5, 
  #                 prov, mom, block, 
  #                 climate_vars_dif)
  # 
  # ## Change to h2o format
  # gbm_data_h2o <- as.h2o(x = gbm_data, key = gbm_data_h2o)
  # 
  # summary(gbm_data_h2o) ## Make sure factors are read correctly
  # 
  # 
  # 
  # ### SPLIT FIT PREDICT 
  #
  # 
  # response <-  "rgr"
  # predictors <- setdiff(names(gbm_data_h2o), c(response, "name"))
  # 
  # ## Split data for training, testing, and validation
  # ## May not be correct since data is not independent
  # 
  # splits <- h2o.splitFrame(
  #   data = gbm_data_h2o, 
  #   ratios = c(0.7,0.15),   ## only need to specify 2 fractions, the 3rd is implied
  #   destination_frames = c("train.hex", "valid.hex", "test.hex"), seed = 1234
  # )
  # train <- splits[[1]]
  # valid <- splits[[2]]
  # test  <- splits[[3]]
  # 
  # 
  # print(paste("Training data has", ncol(train), "columns and", nrow(train), "rows, valid has",
  #             nrow(valid), "rows, test has", nrow(test)))
  # 
  # ### Run gbm model
  # gbm <- h2o.gbm(x = predictors,
  #                y = response,
  #                training_frame = h2o.rbind(train, valid),
  #                # training_frame =  train,
  #                # validation_frame = valid,
  #                nfolds = 5,
  #                
  #                model_id = "gbm",
  #                
  #                ## more trees is better if the learning rate is small enough 
  #                ## here, use "more than enough" trees - we have early stopping
  #                ntrees = 10000,   
  #                
  #                max_depth = 10,
  #                
  #                ## smaller learning rate is better (this is a good value for most datasets, but see below for annealing)
  #                learn_rate=0.01,                                                         
  #                
  #                ## early stopping once the validation MSE doesn't improve by 
  #                # at least 0.01% for 5 consecutive scoring events
  #                stopping_rounds = 5, stopping_tolerance = 1e-4, stopping_metric = "MSE", 
  #                
  #                ## sample 80% of rows per tree
  #                sample_rate = 0.8,                                                       
  #                
  #                ## sample 80% of columns per split
  #                col_sample_rate = 0.8,                                                   
  #                
  #                ## fix a random number generator seed for reproducibility
  #                seed = 1234,                                                             
  #                
  #                ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
  #                score_tree_interval = 10 
  #                
  # )
  # 
  # ###### TESTING
  # 
  # ## Depth 10 is usually plenty of depth for most datasets, but you never know
  # hyper_params = list( 
  #   ## restrict the search to the range of max_depth established above
  #   max_depth = seq(1,30,2),                                      
  #   
  #   ## search a large space of row sampling rates per tree
  #   sample_rate = seq(0.2,1,0.01),                                             
  #   
  #   ## search a large space of column sampling rates per split
  #   col_sample_rate = seq(0.2,1,0.01),                                         
  #   
  #   ## search a large space of column sampling rates per tree
  #   col_sample_rate_per_tree = seq(0.2,1,0.01),                                
  #   
  #   ## search a large space of how column sampling per split should change as a function of the depth of the split
  #   col_sample_rate_change_per_level = seq(0.9,1.1,0.01),                      
  #   
  #   ## search a large space of the number of min rows in a terminal node
  #   min_rows = 2^seq(0,log2(nrow(train))-1,1),                                 
  #   
  #   ## search a large space of the number of bins for split-finding for continuous and integer columns
  #   nbins = 2^seq(4,10,1),                                                     
  #   
  #   ## search a large space of the number of bins for split-finding for categorical columns
  #   nbins_cats = 2^seq(4,12,1),                                                
  #   
  #   ## search a few minimum required relative error improvement thresholds for a split to happen
  #   min_split_improvement = c(0,1e-8,1e-6,1e-4),                               
  #   
  #   ## try all histogram types (QuantilesGlobal and RoundRobin are good for numeric columns with outliers)
  #   histogram_type = c("UniformAdaptive","QuantilesGlobal","RoundRobin")       
  # )
  # 
  # search_criteria = list(
  #   ## Random grid search
  #   strategy = "RandomDiscrete",      
  #   
  #   ## limit the runtime to 60 minutes
  #   max_runtime_secs = 3600,         
  #   
  #   ## build no more than 100 models
  #   max_models = 100,                  
  #   
  #   ## random number generator seed to make sampling of parameter combinations reproducible
  #   seed = 1234,                        
  #   
  #   ## early stopping once the leaderboard of the top 5 models is converged to 0.1% relative difference
  #   stopping_rounds = 5,                
  #   stopping_metric = "MSE",
  #   stopping_tolerance = 1e-3
  # )
  # 
  # grid <- h2o.grid(
  #   ## hyper parameters
  #   hyper_params = hyper_params,
  #   
  #   ## full Cartesian hyper-parameter search
  #   search_criteria = search_criteria,
  #   
  #   ## which algorithm to run
  #   algorithm="gbm",
  #   
  #   ## identifier for the grid, to later retrieve it
  #   grid_id="final_grid",
  #   
  #   ## standard model parameters
  #   x = predictors, 
  #   y = response, 
  #   training_frame = train, 
  #   validation_frame = valid,
  #   
  #   ## more trees is better if the learning rate is small enough 
  #   ## here, use "more than enough" trees - we have early stopping
  #   ntrees = 10000,                                                            
  #   
  #   ## smaller learning rate is better
  #   ## since we have learning_rate_annealing, we can afford to start with a bigger learning rate
  #   learn_rate = 0.05,                                                         
  #   
  #   ## learning rate annealing: learning_rate shrinks by 1% after every tree 
  #   ## (use 1.00 to disable, but then lower the learning_rate)
  #   learn_rate_annealing = 0.99,                                               
  #   
  #   ## fix a random number generator seed for reproducibility
  #   seed = 1234,                                                             
  #   
  #   ## early stopping once the validation AUC doesn't improve by at least 0.01% for 5 consecutive scoring events
  #   stopping_rounds = 5,
  #   stopping_tolerance = 1e-4,
  #   stopping_metric = "MSE", 
  #   
  #   ## score every 10 trees to make early stopping reproducible (it depends on the scoring interval)
  #   score_tree_interval = 10                                                
  # )
  # 
  # ## by default, display the grid search results sorted by increasing logloss (since this is a classification task)
  # grid                                                                       
  # 
  # ## sort the grid models by decreasing AUC
  # sortedGrid <- h2o.getGrid("final_grid", sort_by="MSE", decreasing = FALSE)    
  # sortedGrid
  # 
  # gbm <- h2o.getModel(sortedGrid@model_ids[[1]])
  # 
  # 
  # # Run DRF
  # drf <- h2o.randomForest(x = predictors,
  #                         y = response,
  #                         training_frame = h2o.rbind(train, valid),
  #                         # training_frame =  train,
  #                         # validation_frame = valid,
  #                         nfolds = 5,
  #                         model_id = "rf",
  #                         ntrees            = 250,
  #                         max_depth         = 30)
  # 
  # 
  # # 4- Score on holdout set & report
  # train_rmse_gbm  <- h2o.rmse(gbm, train = TRUE)
  # xval_rmse_gbm   <- h2o.rmse(gbm, xval = TRUE)
  # test_perf_gbm <- h2o.performance(model = gbm, newdata = test)
  # test_rmse_gbm   <- h2o.rmse(object = test_perf_gbm)
  # print(paste0("GBM rmse TRAIN = ", train_rmse_gbm, ", rmse XVAL = ", xval_rmse_gbm, ", rmse TEST = ",
  #              test_rmse_gbm))
  # 
  # train_rmse_drf  <- h2o.rmse(drf, train = TRUE)
  # xval_rmse_drf   <- h2o.rmse(drf, xval = TRUE)
  # test_perf_drf <- h2o.performance(model = drf, newdata = test)
  # test_rmse_drf   <- h2o.rmse(object = test_perf_drf)
  # print(paste0("DRF rmse TRAIN = ", train_rmse_drf, ", rmse XVAL = ", xval_rmse_drf, ", rmse TEST = ",
  #              test_rmse_drf))
  # 
  # 
  # ## Show detailed model summary
  # gbm
  # 
  # drf
  # 
  # 
  # ## Get the Mean Squared Error on the validation set
  # h2o.mse(h2o.performance(gbm, newdata = valid)) ## Best is 0.009
  # 
  # h2o.varimp(gbm)
  # h2o.varimp_plot(gbm, ncol(train))
  # 
  # h2o.varimp(drf)
  # h2o.varimp_plot(drf, ncol(train))
  # 
  # h2o.partialPlot(gbm, h2o.rbind(train, valid))
  # 
  # h2o.saveModel(gbm, "h2o_model")
  # 
  # 
  # 
  # 
  # ### Make predictions
  # 
  # pred_gbm <- h2o.predict(object = gbm,
  #                         newdata = test)
  # 
  # plot(as.data.frame(test)$rgr, as.data.frame(pred_gbm)$predict, pch = 19, las = 1,
  #      xlab = "Observed", ylab = "Predicted")
  # abline(a = 0, b = 1, lwd = 3, col = "steelblue")
  # 
  # 
  # ## Partial plots
  # 
  # h2o.partialPlot(gbm, h2o.rbind(train, valid), c("height_2016_scaled", "DD_0_dif"))
  # h2o.partialPlot(gbm, h2o.rbind(train, valid), "prov", nbins = 100)








