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

#### 10/6/2023: Results from (incomplete) CAR DGP comparing different CAR methods ####
### Timing
files <- grep('2023_10_05',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

res <- combine_results_wrapper(files = files, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

plot_all_methods(res = res, fix_axis = list(F, F, F, F), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP: MCAR')


load(dir(files, full.names = T)[1])
tt = imputed_list[[5]]$df_miss
# ok it's not in there. Is this a memory issue? Ya probably

### Timing
timing = data.table::rbindlist(lapply(imputed_list, '[[', 4))
timing = timing[-c(5,8,9),]
colMeans(timing)
apply(timing, 2, median)
apply(timing, 2, sd)

### ESS of CARstan
res <- combine_results_wrapper(files = files, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'), give_method_name_err = F, ignore_size_err = T, return_unprocessed = T)

ESS = do.call('rbind',(res$CARstan_ESS))
min(colMeans(ESS))
min(ESS)
# Alright this is all good. 

### District level
res <- combine_results_wrapper(files = files, methods = c("y_pred_WF"), rename_vec = c('WF'), give_method_name_err = F, ignore_size_err = T, district_results = T)

#
#### 10/04/2023: Investigating freqGLMepi in CAR DGP ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files1 <- grep('car331', files, value = T)

res <- combine_results_wrapper(files = files1, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))
# from 6 to 26 object sizes less than half...why?

res2 <- combine_results_wrapper(files = files1, methods = c("y_pred_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF','CAR', 'freqGLM'))
# also from 6 to 26. Why did I think freqGLMepi was so bad?

res <- combine_results_wrapper(files = files2, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))

## The data I need to test this.
# load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_car331_beta6_n025_20230603/simulated_data.RData" )

#
#### 9/20/2023 (Kiko's bday!): Getting groups of results ####
### Get the file names
files_MNAR <- grep('nost_beta43_n025_2023_02_22',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files_MAR <- grep('nost_beta43_n025_2023_02_26',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)[1:6]

tmp <- grep('nost_2022_11',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)
files_MCAR = grep('mcar', tmp, value = T)[1:6]

### Get the results!
res_MCAR = combine_results_wrapper(files_MCAR, methods = c("y_pred_CCA_WF","y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), results_by_point = F)

res_MAR = combine_results_wrapper(files_MAR, methods = c("y_pred_CCA_WF","y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), results_by_point = F)

res_MNAR = combine_results_wrapper(files_MNAR, methods = c("y_pred_CCA_WF","y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), results_by_point = F)

#save(res_MAR, res_MCAR, res_MNAR, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/results_WF_3MGP_09202023.RData')

load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/results_WF_3MGP_09202023.RData')

p1 <- plot_all_methods(res = res_MCAR, fix_axis = list(F, F, F, F), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF: MCAR')

p2 <- plot_all_methods(res = res_MAR, fix_axis = list(F, F, F, F), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF: MAR')

p3 <- plot_all_methods(res = res_MNAR, fix_axis = list(F, F, F, F), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF: MNAR')

plot_grid(plotlist = list(p1, p2, p3), ncol = 1)

ggsave('figures/results_MCAR_MAR_MNAR.png', height = 7, width = 10)

#
#### 9/19/2023: Testing simulation-level results ####
files <- grep('beta6_n025_20230603',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

{
load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar01_car331_beta6_n025_20230603/sim_results_p0.1_mcar_1(50).RData")

# cutting to 9 (it's confusing to have a list of 20 with 20 facilities)
imputed_list = imputed_list[1:9]

df = imputed_list[[1]]$df_miss

imputed_list = lapply(imputed_list, '[[', 1)

#res = calculate_metrics_by_sim(df, methods = c('y_pred_WF','y_pred_CCA_CAR', 'y_pred_CCA_freqGLMepi'), imputed_only = F)

res = calculate_metrics(imputed_list, methods = c('y_pred_WF','y_pred_CCA_CAR'), date = NULL, min_date = '2020-01-01', results_by_point = T)

res2 = calculate_metrics_by_point(imputed_list, methods = c('y_pred_WF','y_pred_CCA_CAR'), min_date = '2020-01-01', imputed_only = F)

identical(res$date, res2$date)
identical(res$facility, res2$facility)
# good. Those match

identical(res$bias, res2$bias)
identical(res$outbreak_detection10, res2$outbreak_detection10)
identical(res$coverage95, res2$coverage95)
}
# ok it's all good. DELETE

res = combine_results_wrapper(files, methods = c("y_pred_WF","y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), results_by_point = F)

tt <- sapply(imputed_list, function(xx){
  return(any(grepl('freqGLMepi', names(xx$df_miss))))
})


p1 <- plot_all_methods(res = res, fix_axis = list(F, F, F, F), add_lines = list(0.95, F, F, F), metrics = c('coverage95', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), methods = c("y_pred_WF","y_pred_CCA_CAR", 'y_pred_CCA_freqGLMepi'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP')

p2 <- plot_all_methods(files, fix_axis = list(F, F, F, F), metrics = c('coverage95', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), methods = c("y_pred_WF","y_pred_CCA_CAR"), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = T)

#
#### 7/26/2023: Comparing the CAR n.sample and thinning performance ####
files <- grep('car_thin_comparison_0_20230725',dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_df <- process_CAR_params_wrapper(files, all_betas = F, CAR_names = c('CAR_summary', 'CAR_summary2', 'CAR_summary3', 'CAR_summary4'))


tmp = res_df %>% filter(param == 'betas')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'tau2')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'alpha')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'rho.S')
tmp[,c(1,8:11)]


#
#### 7/25/2023: Comparing the CAR param convergence for different sample sizes and plotting metrics ####
files <- grep('car_prior_comparison_0_20230720',dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_df <- process_CAR_params_wrapper(files, all_betas = F, CAR_names = c('CAR_summary', 'CAR_summary2', 'CAR_summary3', 'CAR_summary4', 'CAR_summary5'))
plot_CAR_params(res_df = res_df)

tmp = res_df %>% filter(param == 'tau2')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'betas')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'alpha')
tmp[,c(1,8:11)]

tmp = res_df %>% filter(param == 'rho.S')
tmp[,c(1,8:11)]

res = combine_results_wrapper(files, methods = c("y_pred_CAR_none", "y_pred_CAR_WF1", "y_pred_CAR_WF10", "y_pred_CAR_WF100", "y_pred_CAR_constant"), rename_vec = c('none', 'WF_1' ,'WF_10', 'WF_100', 'constant'))

plot_all_methods(res = res, fix_axis = F)

res2 = res %>% 
  filter(method != 'constant')

plot_all_methods(res = res2, fix_axis = F)

#
#### 7/09/2023: Comparing R values ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_full = combine_results_wrapper(files, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))

res100 <- lapply(1:10, function(ii){
  seq = (1 + (ii-1)*100):(ii*100)
  tmp = combine_results_wrapper(files, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), subset_results = seq)
  return(tmp)
})

res200 <- lapply(1:5, function(ii){
  seq = (1 + (ii-1)*200):(ii*200)
  tmp = combine_results_wrapper(files, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), subset_results = seq)
  return(tmp)
})

res_full$method = paste0(res_full$method, '_full')
res1 = res200[[1]]
res1$method = paste0(res1$method, '_200')
res2 = res200[[2]]
res2$method = paste0(res2$method, '_200_2')
res3 = res200[[3]]
res3$method = paste0(res3$method, '_200_3')

res = rbind(res_full, res1, res2, res3)

ind_WF = grep('WF', res$method)
res_WF = res[ind_WF,]
res_CAR = res[-ind_WF,]

plot_all_methods(res = res_WF, fix_axis = F)
plot_all_methods(res = res_CAR, fix_axis = F)


# How to subset??
# will need to adapt the combine_results and combine_results_wrapper function so that they can take subsets. 
# This will be inefficient in terms of re-combining the results each time, but it's easiest for now. Actually that *may* not be true


#### 6/09/2023: Getting CAR parameter convergence ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_df <- process_CAR_params_wrapper(files, all_betas = F)
plot_CAR_params(res_df = res_df)
#plot_CAR_params(files = files)

res_df <- process_CAR_params_wrapper(files, all_betas = T)
res2 = res_df %>%
  filter(!(param %in% c('tau2','rho.S','alpha','betas')))

table(res2$prop_missing, res2$mean_acceptance)

#
#### 6/08/2023: Plotting CAR DGP ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

f1 <- grep('simulated_data', dir(files[1], full.names = T), value = T)
load(f1)
df = lst$df_list[[1]]
plot_facility_fits(df, outbreak_points = c(3,5,10))

mean(df$y) # 248
median(df$y) # 107

#dist_df = lst$district_list[[1]]

load(dir(files[1], full.names = T)[1])
df2 = imputed_list[[1]]$df_miss
plot_facility_fits(df2, 
                   methods = c('y_pred_WF', 'y_pred_CCA_CAR'), 
                   imp_names = c('WF', 'CAR'),
                   outbreak_points = c(3,5,10))

### Now for some MCAR data
load('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/mcar02_nost_beta6_n025_2023_03_16/simulated_data.RData')
df = lst$df_list[[1]]

mean(df$y) # 185
median(df$y) # 105
df$y_var = df$y_exp

plot_facility_fits(df, outbreak_points = c(3,5,10))

#
#### 6/05/2023: Analyzing CAR DGP facility AND district results ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res = combine_results_wrapper(files, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))
res_dst = combine_results_wrapper(files, methods = c("y_pred_WF"), rename_vec = c('WF'), district_results = T)

tt <- plot_all_methods(res = res, methods = c("y_pred_WF", "y_pred_CCA_CAR"), fix_axis = F, rename_vec = c('WF','CAR'))

tt2 <- plot_all_methods(res = res_dist, methods = c("y_pred_WF"), rename_vec = c('WF'), fix_axis = F)





files1 = files[1]
lst_full <- combine_results(input_folder = files1, return_lst = T, results_file = NULL)

#
#### 5/30/2023: Analyzing CAR DGP expectation and variance ####

# after simulating the data
tt <- do.call('rbind', lapply(lst_WF$df_list, function(df){
  res = data.frame(mean_observed = mean(df$y),
                   mean_expected = mean(df$y_exp))
  res
}))
colMeans(tt)

vars <- do.call('rbind', lapply(lst_WF$df_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_var)
  return(var(z))
}))
plot(density(vars))
mean(vars)

## now for the district analysis
tt <- do.call('rbind', lapply(lst_WF$district_list, function(df){
  res = data.frame(mean_observed = mean(df$y),
                   mean_expected = mean(df$y_exp))
  res
}))
colMeans(tt)

vars <- do.call('rbind', lapply(lst_WF$district_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_var)
  return(var(z))
}))
plot(density(vars))
mean(vars)

# after simulating the data
tt <- do.call('rbind', lapply(lst_CAR$df_list, function(df){
  res = data.frame(mean_observed = mean(df$y),
                   mean_expected = mean(df$y_exp))
  res
}))
colMeans(tt)

vars <- do.call('rbind', lapply(lst_CAR$df_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_var)
  return(var(z))
}))
plot(density(vars))
mean(vars)

## now for the district analysis
tt <- do.call('rbind', lapply(lst_CAR$district_list, function(df){
  res = data.frame(mean_observed = mean(df$y),
                   mean_expected = mean(df$y_exp))
  res
}))
colMeans(tt)

vars <- do.call('rbind', lapply(lst_CAR$district_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_var)
  return(var(z))
}))
plot(density(vars))
mean(vars)

vars <- do.call('rbind', lapply(lst_CAR$district_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_var_ind)
  return(var(z))
}))
mean(vars)
# too small

vars <- do.call('rbind', lapply(lst_CAR$district_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_cov)
  return(var(z))
}))
mean(vars)
# too small

vars <- do.call('rbind', lapply(lst_CAR$district_list, function(df){
  z = (df$y - df$y_exp)/sqrt(df$y_cov/2 + df$y_var_ind)
  return(var(z))
}))
mean(vars)

#
#### 5/22/2023: CAR DGP results ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files1 <- grep('car331', files, value = T)
files2 <- grep('car731', files, value = T)

res <- combine_results_wrapper(files = files1, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))

tt <- plot_all_methods(files, methods = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), fix_axis = F)

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


