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

#### 7/18/2023: Assessing CAR thinning comparison ####
files <- grep('car_thin',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_full = combine_results_wrapper(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))

#### 7/09/2023: Comparing R values ####
files <- grep('20230603',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_full = combine_results_wrapper(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))

res100 <- lapply(1:10, function(ii){
  seq = (1 + (ii-1)*100):(ii*100)
  tmp = combine_results_wrapper(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), subset_results = seq)
  return(tmp)
})

res200 <- lapply(1:5, function(ii){
  seq = (1 + (ii-1)*200):(ii*200)
  tmp = combine_results_wrapper(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), subset_results = seq)
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
                   imp_vec = c('y_pred_WF', 'y_pred_CCA_CAR'), 
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

res = combine_results_wrapper(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'))
res_dst = combine_results_wrapper(files, imp_vec = c("y_pred_WF"), rename_vec = c('WF'), district_results = T)

tt <- plot_all_methods(res = res, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), fix_axis = F, rename_vec = c('WF','CAR'))

tt2 <- plot_all_methods(res = res_dist, imp_vec = c("y_pred_WF"), rename_vec = c('WF'), fix_axis = F)





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
files <- grep('20230517',dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

files1 <- grep('car331', files, value = T)
files2 <- grep('car731', files, value = T)

tt <- plot_all_methods(files, imp_vec = c("y_pred_WF", "y_pred_CCA_CAR"), rename_vec = c('WF','CAR'), fix_axis = F)

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


