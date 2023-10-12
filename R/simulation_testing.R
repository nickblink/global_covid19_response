current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)

#### Combining all results (WF w/ MCAR, MAR, MNAR; freqGLM MCAR; CAR_331 MCAR, CAR_33025 MCAR) ####
files <- grep('2023_10_11',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

## Pull in MCAR results
{
files_MCAR = grep('mcar', files, value = T)

res_MCAR <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

res_MCAR$results$sim = 'WF_MCAR' 
res_MCAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 

res_MCAR_district <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)

res_MCAR_district$results$sim = 'WF_MCAR' 
res_MCAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 
}

## Pull in MNAR results
{
files_MNAR = grep('mnar', files, value = T)

res_MNAR <- combine_results_wrapper(files = files_MNAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

res_MNAR$results$sim = 'WF_MNAR'
res_MNAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MNAR' 

res_MNAR_district <- combine_results_wrapper(files = files_MNAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)

res_MNAR_district$results$sim = 'WF_MNAR'
res_MNAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MNAR' 
}

## Pull in CAR results
{
files <- grep('2023_10_06',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files_331 = grep('car331', files, value = T)
files_33025 = grep('car33025', files, value = T)

res_CAR331 <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res_CAR331$results$sim = 'CAR_331'
res_CAR331$results$sim_long = 'CAR_DGP(0.3,0.3,1,b0=6,b1=-0.25);MCAR' 
res_CAR331_district <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'), district_results = T)

res_CAR331_district$results$sim = 'CAR_331'
res_CAR331_district$results$sim_long = 'CAR_DGP(0.3,0.3,1,b0=6,b1=-0.25);MCAR' 

res_CAR33025 <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res_CAR33025$results$sim = 'CAR_33025'
res_CAR33025$results$sim_long = 'CAR_DGP(0.3,0.3,0.25,b0=6,b1=-0.25);MCAR' 
res_CAR33025_district <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'), district_results = T)

res_CAR33025_district$results$sim = 'CAR_33025'
res_CAR33025_district$results$sim_long = 'CAR_DGP(0.3,0.3,0.25,b0=6,b1=-0.25);MCAR' 
}

## Join these together
res_facility = rbind(res_MCAR$results, res_MNAR$results, res_CAR331$results, res_CAR33025$results)
res_district = rbind(res_MCAR_district$results, res_MNAR_district$results, res_CAR331_district$results, res_CAR33025_district$results)


#
#### Processing CAR results ####
files <- grep('2023_10_06',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files_331 = grep('car331', files, value = T)
files_33025 = grep('car33025', files, value = T)

res1 <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res2 <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

table(res1$params$tau2_DGP)
table(res2$params$tau2_DGP)
# Good

plot_all_methods(res = res1$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=1): MCAR')

plot_all_methods(res = res2$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)),  add_lines = list(0.95, F, F, F), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=0.25): MCAR')


res1_district <- combine_results_wrapper(files = files_331, district_results = T, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res2_district <- combine_results_wrapper(files = files_33025, district_results = T, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

plot_all_methods(res = res1_district$results, fix_axis = list(ylim(0,1), ylim(0,1), ylim(0.3, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=1): MCAR')

plot_all_methods(res = res2_district$results, fix_axis = list(ylim(0,1), ylim(0,1), ylim(0.3, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=.025): MCAR')

#save(res1, res2, res1_district, res2_district, file = 'C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/results_CARDGP_10102023.RData')

     
#
#### (one time run) Fix file naming ####
files <- grep('2023_10_06',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

res <- combine_results_wrapper(files = files, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

params = res$params

ind <- grep('\\(1\\)', params$file)
table(params$tau2_DGP[ind])
table(params$tau2_DGP[-ind])

# files_to_move = params$file[ind]
# for(f in files_to_move){
#   new_file = gsub('car331','car33025', f)
#   new_file = gsub('\\(1\\)','', new_file)
#   file.rename(f, new_file)
# }

#### Analyzing the missingness pattern ####
plot_missingness <- function(df_miss){
  df_spread = df_miss %>%
    select(date, facility,y) %>%
    tidyr::spread(., facility, y)
  
  tmp2 = df_spread[,-c(1)]
  for(col in colnames(tmp2)){
    tmp2[,col] = as.integer(is.na(tmp2[,col]))
  }
  
  tmp2 = as.matrix(tmp2)
  rownames(tmp2) = as.character(df_spread$date)
  heatmap(tmp2, keep.dendro = F, Rowv = NA, )
  
  gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = F, Colv = F, xlab = 'facilities', trace = 'none', key = F)
}

lst <- simulate_data(district_sizes = c(4, 6, 10), R = 1, end_date = '2019-12-01')

df <- lst$df_list[[1]] 

# simulation function!
df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0, alpha = 0, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.9, alpha = 0, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0, alpha = 0.9, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.9, alpha = 0.9, tau = 3)

plot_missingness(df_miss)


