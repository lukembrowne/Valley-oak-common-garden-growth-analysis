

# Load libraries ----------------------------------------------------------

source("./01_clean-process-explore-data.R")
source("./03_adding-climate-data.R")

dim(dat_all_scaled)

# install_github("jdstorey/qvalue")
library(gamm4)
library(visreg)
library(MuMIn)
library(beepr)
library(patchwork)

# function to transform from raw to scaled variables
forward_transform <- function(x, var, means, sds){
  ( x - means[var]) / sds[var]
}

# Function to back transform variables
back_transform <- function(x, var, means, sds){
  (x * sds[var]) + means[var]
}


# Function for plotting predictions ---------------------------------------

plot_pred <- function(xvar,
                      xlab = "",
                      add_points = TRUE,
                      add_residuals = FALSE,
                      save = FALSE,
                      path = "./",
                      filename = "",
                      main = ""){
  
  # Test to see what type of variable it is - factor or numeric
  type <- ifelse(is.factor(pull(dat_all_scaled, xvar)), "factor", "numeric")
  
  # Make predictions
  v <- visreg(gam_all, xvar = xvar, 
              scale = "response", rug = FALSE, ylab = "5 yr height (cm)", plot = FALSE)
  
  dat_all_scaled_trans <- dat_all_scaled
  
  
  # If xvariable is not a factor, backtranform to original scale for plotting
  if(type == "numeric"){
    # Transform x axis + assign variables for lo and hi CI
    v$fit <- v$fit %>%
      # Back transform X axis
      mutate(!!xvar := back_transform(x = get(xvar),
                                      var = xvar,
                                      means = scaled_var_means_all,
                                      sds = scaled_var_sds_all))
    
    dat_all_scaled_trans <- dat_all_scaled %>%
      # Back transform X axis
      mutate(!!xvar := back_transform(x = get(xvar),
                                      var = xvar,
                                      means = scaled_var_means_all,
                                      sds = scaled_var_sds_all))
  } # End factor if
  
  
  # Add partial residuals
  dat_all_scaled_trans$resid <- v$res$visregRes
  
  ## Main GGPLOT
  gg <- 
    ggplot(dat_all_scaled_trans, aes(x = get(xvar), y = height_2017)) + 
    xlab(xvar) + ylab("5 yr height (cm)") +
    ggtitle(main) +
    xlab(xlab) +
    theme_bw()
  
  # Add raw points to graph?
  if(add_points == TRUE){
    gg <- gg + 
      geom_jitter(data = dat_all_scaled_trans, aes(x = get(xvar), y = height_2017),
                  size = 0.5, alpha = 0.5)
  }
  
  ## Add partial residuals to graph?
  if(add_residuals == TRUE){
    gg <- gg + 
      geom_point(data = dat_all_scaled_trans, aes(x = get(xvar), y = resid),
                 size = 0.5, alpha = 0.5)
  }
  
  
  # If numeric, add points and lines
  if(type == "numeric"){
    gg <- gg + 
      geom_ribbon(data = v$fit, aes(ymin = visregLwr, ymax = visregUpr), fill = "forestgreen", alpha = 0.75) +
      geom_line(data = v$fit, aes(x = get(xvar), y = visregFit), lwd = 2) 
  }
  
  # If factor, make point_range
  if(type == "factor"){
    gg <- gg + 
      geom_pointrange(data = v$fit, aes(x = get(xvar), y = visregFit, 
                                        ymin = visregLwr, ymax = visregUpr,
                                        fill = get(xvar)),
                      size = 1.5, shape = 21) +
      guides(fill=FALSE) # Remove legend
  }
  
  
  
  # Add verticle line at 0 if a climate transfer function
  if(grepl(pattern = "_dif", x = xvar)){
    gg <- gg + geom_vline(aes(xintercept = 0), lty = 2) 
  }
  
  # General theme and margins
  gg <- gg + theme(plot.margin = (margin(1.5,1.5,1.5,1.5, "cm")),
                   text = element_text(size = 15))
  
  # Printing to file
  if(save){
    ggsave(paste0(path, filename))
  }
  
  print(gg)
  
  return(gg)
  
} # End plot_pred function


# Numerical    

# Run FULL models with all seedlings (no genomic data) --------------------

# Convert factors
dat_all_scaled$section <- factor(dat_all_scaled$section)
dat_all_scaled$section_block <- factor(dat_all_scaled$section_block) 
dat_all_scaled$accession <- factor(dat_all_scaled$accession)
dat_all_scaled$site <- factor(dat_all_scaled$site)

var = "tmax_sum_dif"


# Set formula for gam
  form <- formula(paste0("height_2017 ~ section_block + s(height_2014, bs =\"cr\") + s(", var, " , bs = \"cr\") + s(accession, bs = \"re\")"))
  
  
  # With all 5,000+ seedlings
  gam_all <- bam(formula = form,
                 data = dat_all_scaled,
                 discrete = TRUE, 
                 nthreads = 8,
                 method = "fREML", family = "gaussian", 
                 control = list(trace = FALSE))
  
  summary(gam_all)

  # Plot overall model fit
  test_fit <- dat_all_scaled
  test_fit$pred <- gam_all$fitted.values
  

  ggplot(test_fit, aes(x = pred, y = height_2017, bg = section_block)) + geom_point(alpha = 0.75, pch = 21) + theme_bw(15) +
    geom_abline(slope = 1, intercept = 0, lwd = 1.5, col = "forestgreen")
  
  visreg(gam_all, partial = TRUE, ylab = "5 yr height (cm)")
  
  visreg(gam_all, partial = FALSE, ylab = "5 yr height (cm)")
  
  gam.check(gam_all)
  
  
  save_flag = FALSE
  path <- "./figs_tables/"
  
  
  p = plot_pred(xvar = var, xlab = "Maximum Temperature Summer Transfer Distance",
                add_points = FALSE, add_residuals = FALSE,
                save = FALSE, path = path, filename = "")
  
  
  # print(p)
  
  p =  plot_pred(xvar = var, add_points = FALSE, add_residuals = TRUE,
                 save = FALSE, path = path, filename = "")
  
  # print(p)
  
  dev.off()
  
