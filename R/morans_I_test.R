current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')
library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)

if(file.exists('C:/Users/Admin-Dell')){
  res_dir = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results"
}else{
  res_dir = "C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results"
}

setwd(res_dir)

load('WF_Poisson_NB_comparison_WFDGP_R100_04072024/sim_results_1(1).RData')

tt <- imputed_list[[1]]$df_miss

# Finding covariance within each time point.
N = length(unique(tt$facility))
Ttime = length(unqiue(tt$date))

tmp <- tt %>% filter(date == '2016-01-01')
var_y = var
y_exp = mean(tmp$y)
mean((tmp$y - y_exp)^2)

can use ape::moran.i to input adjacecny and data for each time point.

# I'm confused. honestly I'm hungry. But if I do this and average across all time points that could add variance because of the high-missing time points.
# If I compute with all points at once, then wouldn't that show a high I even if there's no spatial correlation? Like because of equal time points having similar trends. Yes. But the equal time points won't be used more than unequal time points. Boy there's gotta be a better way to say this.

