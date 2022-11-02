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


## 7 different approaches. 
# WF full data
# WF CCA
# Bayes CCA
# FreqGLM_epi CCA
# Bayes + WF
# FreqGLM_epi + WF
# MICE + WF

#### Creating function for combining results ####
input_folder <- 'cluster_code/mnar_2022_11_01/'
results_file <- 'results/simulation_noST_MNARp2_R100_11012022.Rdata'

combine_results('cluster_code/mcar_2022_11_01/', 'results/TESTING_TESTING_123.RData') -> tmp

tt <- c()
ss <- c()
vv <- c()
for(i in 1:100){
  tt <- c(tt, mean(tmp[[i]]$y, na.rm = T))
  ss <- c(ss, mean(tmp[[i]]$y_true))
  vv <- c(vv, mean(tmp[[i]]$y_pred_CCA_WF))
}
length(unique(vv)) # good


#### 10/21/2022: Comparing WF MCAR, MAR, MNAR ####
imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

# MCAR
load('results/simulation_noST_MCARp2_R100_10162022.RData')
res_MCAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method %in% c('y_pred_baseline_WF', 'y_pred_CCA_WF'))

res_MCAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MCAR', res_MCAR$method)

# MAR
rm(imputed_list)
load('results/simulation_noST_MARp2_R80_10162022.RData')
res_MAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_CCA_WF')

res_MAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MAR', res_MAR$method)

# MNAR
rm(imputed_list)
load('results/simulation_noST_MNARp2_R100_10202022.RData')
res_MNAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_CCA_WF')

res_MNAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MNAR', res_MNAR$method)

# combine and plot
res <- rbind(res_MCAR, res_MNAR, res_MAR)

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF_MCAR",  "y_pred_CCA_WF_MAR",  "y_pred_CCA_WF_MNAR")
rename_vec = c('WF full data', 'WF CCA MCAR', 'WF CCA MAR', 'WF CCA MNAR')
color_vec = c('red', 'blue', 'green', 'orange')

plot_metrics_by_point(res = res, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

#### 10/20/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MNAR ####

R = 100

results_file <- 'results/simulation_noST_MNARp2_R100_10202022.RData'

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MNAR_sim(df, p = 0.2, direction = 'upper', gamma = 1, by_facility = T)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF baseline (complete data) imputation
    df_miss <- WF_baseline(df_miss, R_PI = 100)
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
    
    if(i %% 5 == 0){
      print(sprintf('saving results for i = %s',  i))
      save(imputed_list, i, file = results_file)
    }
  }
}) # ~15hrs

# save(imputed_list, file = results_file)

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_baseline_WF') %>%
  group_by(date) %>%
  summarize(num_missing = sum(num_missing))

# number missing points across simulations
plot(res$date, res$num_missing, ylab = 'number of times missing', xlab = 'date', main = 'missingness from R = 100 with 20 facilities')

for(i in 1:10){
  test <- imputed_list[[i]]  %>% filter(is.na(y))
  print(max(test$date))
}

#### 10/16/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MCAR ####

R = 100

results_file <- 'results/simulation_noST_MCARp2_R200_10162022.RData'

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MCAR_sim(df, p = 0.2, by_facility = T)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF baseline (complete data) imputation
    df_miss <- WF_baseline(df_miss, R_PI = 100)
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
    
    if(i %% 5 == 0){
      print(sprintf('saving results for i = %s',  i))
      save(imputed_list, i, file = results_file)
    }
  }
}) # 15 hrs R = 100

# save(imputed_list, file = results_file)

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, min_date = '2020-01-01', imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  arrange(coverage95) %>%
  filter(method == 'y_pred_baseline_WF')


#### 10/16/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MAR ####

R = 80

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.7, alpha = 0.7, tau = 3)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF baseline (complete data) imputation
    df_miss <- WF_baseline(df_miss, R_PI = 100)
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
  }
}) # Started at 8pm on 10/16. Should be 12 hours, so ready by the morn

# save(imputed_list, file = 'results/simulation_noST_MARp2_R200_10162022.RData')

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, min_date = '2020-01-01', imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  arrange(coverage95) %>%
  filter(method == 'y_pred_baseline_WF')


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


