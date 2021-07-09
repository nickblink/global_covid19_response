### Now doing this as an R script rather than Rmd because it's easier to work with.
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)


#### Comparing CARBayes models ####

R = 50

system.time({
  lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = R, rho = 0.3, alpha = 0.5, tau = 0.5)
  
  imp_vec = c("y_CB_fixed", "y_CB_intercept", "y_CB_facility")
  rename_vec = imp_vec
  color_vec = c('red','blue','green')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df = lst$df_list[[i]]
    
    # simulation function!
    df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
    
    # run the CARBayes imputation with constant coeffs across all facilities
    CAR_list1 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'fixed')
    df_miss = CAR_list1$facility_df
    colnames(df_miss) = gsub('CARBayes_ST', 'CB_fixed', colnames(df_miss))
    
    #  run the CARBayes imputation with different intercepts by facility
    CAR_list2 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_intercept')
    df_miss = CAR_list2$facility_df
    colnames(df_miss) = gsub('CARBayes_ST', 'CB_intercept', colnames(df_miss))
    
    #  run the CARBayes imputation with different coeffs by facility
    CAR_list3 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')
    df_miss = CAR_list3$facility_df
    colnames(df_miss) = gsub('CARBayes_ST', 'CB_facility', colnames(df_miss))
    
    imputed_list[[i]] = df_miss
  }
  
  #plot_metrics_by_point(imputed_list, imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CARBayes_ST"), min_missing = 0, rename_vec = c('glmFreq','glmFreq_epi','CARBayes'), color_vec = c('red','blue','green'))
  
  pfit = plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec, color_vec = color_vec)
  
  p1 <- plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, imputed_only = F, min_missing = 0)
  
  p2 <- p1 <- plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, imputed_only = T, min_missing = 0)
})

# interestingly, it seems like the intercept does basically as good as the diff. facility ones, but I simply don't buy it

#### MCAR p = 0.2 spatio-temporal ####

R = 5

system.time({
lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = R, rho = 0.3, alpha = 0.5, tau = 0.5)

imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CARBayes_ST")
rename_vec = c('glmFreq','glmFreq_epi','CARBayes')
color_vec = c('red','blue','green')

imputed_list = list()
res_full = res_imputed = NULL
for(i in 1:R){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
  
  # run the periodic imputation
  freqGLMepi_list = freqGLMepi_imputation(df_miss, prediction_intervals = 'bootstrap', R_PI = 100, verbose = F) 
  df_miss = freqGLMepi_list$df
  
  # run the periodic imputation
  periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
  df_miss = periodic_list$df
  
  # run the CARBayes imputation
  CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')
  df_miss = CAR_list$facility_df
  
  imputed_list[[i]] = df_miss
}

#plot_metrics_by_point(imputed_list, imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CARBayes_ST"), min_missing = 0, rename_vec = c('glmFreq','glmFreq_epi','CARBayes'), color_vec = c('red','blue','green'))

pfit = plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

p1 <- plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, imputed_only = F, min_missing = 10, rename_vec = rename_vec)

p2 <- p1 <- plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, imputed_only = T, min_missing = 10, rename_vec = rename_vec)
})

## 20 minutes for 50 iterations 
# save(imputed_list, p1, p2, file = 'results/simulation_MCARp2_R500_res.RData')

#
#### Assessing the simulated data ####

# lst <- simulate_data(district_sizes = c(2,3,4,5,7), n = 2)
# df = lst$df

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 2, rho = 0, alpha = 0.5, tau = 1)
df = lst$df_list[[2]]

lst2 <- simulate_data(district_sizes = c(4), R = 2)
df2 = lst2$df_list[[2]]

par(mfrow = c(2,2))
for(f in sample(unique(df$facility), 4)){
  tmp = df %>% filter(facility == f)
  tmp2 = df2 %>% filter(facility == f)
  
  plot(tmp$date, tmp$y, type = 'l', main = f, ylim = c(0, 1.2*max(tmp$y)))
  lines(tmp2$date, tmp2$y, col = 'red')
}

#### Assessing MAR simulation ####
lst <- simulate_data_spatiotemporal(district_sizes = c(10), R = 1, rho = 0.3, alpha = 0.5, tau = 0.5)
df = lst$df_list[[1]]


### no correlation
df_miss = MAR_spatiotemporal_sim(df, p = 0.3, rho = 0, alpha = 0, tau = 3, by_facility = T)
df_spread = df_miss %>%
  dplyr::select(date, district, facility, y) %>%
  tidyr::spread(facility, y)

tmp2 = df_spread[,-c(1,2)]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = 1 - as.matrix(tmp2)

gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = NA, Colv = T, xlab = 'facilities', trace = 'none', key = F)


### spatial correlation
df_miss = MAR_spatiotemporal_sim(df, p = 0.3, rho = 0.9, alpha = 0, tau = 3, by_facility = T)
df_spread = df_miss %>%
  dplyr::select(date, district, facility, y) %>%
  tidyr::spread(facility, y)

tmp2 = df_spread[,-c(1,2)]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = 1 - as.matrix(tmp2)

gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = NA, Colv = T, xlab = 'facilities', trace = 'none', key = F)


### temporal correlation
df_miss = MAR_spatiotemporal_sim(df, p = 0.3, rho = 0, alpha = 0.9, tau = 3, by_facility = T)
df_spread = df_miss %>%
  dplyr::select(date, district, facility, y) %>%
  tidyr::spread(facility, y)

tmp2 = df_spread[,-c(1,2)]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = 1 - as.matrix(tmp2)

gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = NA, Colv = T, xlab = 'facilities', trace = 'none', key = F)

### spatial and temporal correlation
df_miss = MAR_spatiotemporal_sim(df, p = 0.3, rho = 0.9, alpha = 0.9, tau = 3, by_facility = T)
df_spread = df_miss %>%
  dplyr::select(date, district, facility, y) %>%
  tidyr::spread(facility, y)

tmp2 = df_spread[,-c(1,2)]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = 1 - as.matrix(tmp2)

gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = NA, Colv = T, xlab = 'facilities', trace = 'none', key = F)


#
#### MAR p = 0.2 spatio-temporal ####

R = 5
system.time({
lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = R, rho = 0.5, alpha = 0.5, tau = 0.5)

imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')

imputed_list = list()
res_full = res_imputed = NULL

imputed_list = lapply(1:R, function(i){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.5, alpha = 0.5, tau = 3, by_facility = T)
  
  # run the periodic imputation
  periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
  df_miss = periodic_list$df
  
  # run the CARBayes imputation
  CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)
  df_miss = CAR_list$facility_df
  
  df_miss
})

p1 <- plot_metrics_by_point(imputed_list, rename_vec = c('glmFreq','CARBayes'), imputed_only = F, min_missing = 50)

p2 <- plot_metrics_by_point(imputed_list, rename_vec = c('glmFreq','CARBayes'), imputed_only = T, min_missing = 50)

}) # 38m for R = 500

# save(imputed_list, p1, p2, file = 'results/simulation_MARp2_R500_res.RData')

#### MNAR p = 0.2 spatio-temporal ####

R = 2000
system.time({
  lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = R, rho = 0.5, alpha = 0.5, tau = 0.5)
  
  imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  
  imputed_list = lapply(1:R, function(i){
    df = lst$df_list[[i]]
    
    # simulation function!
    df_miss = MNAR_sim(df, p = 0.2, direction = 'upper', gamma = 1, by_facility = T)
    
    # run the periodic imputation
    periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
    df_miss = periodic_list$df
    
    # run the CARBayes imputation
    CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)
    df_miss = CAR_list$facility_df
    
    df_miss
  })
  
  p1 <- plot_metrics_by_point(imputed_list, rename_vec = c('glmFreq','CARBayes'), imputed_only = F, min_missing = 50)
  
  p2 <- plot_metrics_by_point(imputed_list, rename_vec = c('glmFreq','CARBayes'), imputed_only = T, min_missing = 50)
}) # 151 minutes

# save(imputed_list, p1, p2, file = 'results/simulation_MNARp2_R2000_res.RData')

######## GRAVEYARD #########
#### MCAR p = 0.1 ####

df <- simulate_data(district_sizes = c(2,3,4,5,7))

# simulation function!
df_miss = MCAR_sim(df, p = 0.1, by_facility = T)

# run the periodic imputation
periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
df_miss = periodic_list$df

# run the CARBayes imputation
CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)
df_miss = CAR_list$facility_df

# get the county level imputation results
county_miss = merge(periodic_list$county_fit, CAR_list$county_df, by = 'date')

# make the plots for the results
imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')
imp_names = c('parametric','CAR')
color_vec = c('red','blue')
p1 <- plot_facility_fits(df_miss, imp_vec = imp_vec, imp_names = imp_names, color_vec = color_vec)
p2 <- plot_county_fits(county_miss, imp_vec, imp_names = imp_names, color_vec, title = 'p = 0.1')

# ggsave(plot = p1, filename = 'C:/Users/nickl/Documents/global_covid19_response/figures/MCAR_p01_simulation_fits_05172021.png', width = 15, height = 10)
# ggsave(plot = p2, filename = 'C:/Users/nickl/Documents/global_covid19_response/figures/MCAR_p01_simulation_countyFit_05172021.png', width = 10, height = 7)

# get the metrics to evaluate each method
res <- calculate_metrics(df_miss, imp_vec,imputed_only = F)
res2 <- calculate_metrics(df_miss, imp_vec,imputed_only = T)
res3 <- calculate_metrics(county_miss, imp_vec,imputed_only = F, median_estimate = T)


#### MCAR p = 0.5 ####

df <- simulate_data(district_sizes = c(2,3,4,5,7))

# simulation function!
df_miss = MCAR_sim(df, p = 0.5, by_facility = T)

# run the periodic imputation
periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
df_miss = periodic_list$df

# run the CARBayes imputation
CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)

# get the facility and county level imputation results
df_miss = CAR_list$facility_df
county_miss = merge(periodic_list$county_fit, CAR_list$county_df, by = 'date')

# make the plots for the results
imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')
color_vec = c('red','blue')
p3 <- plot_facility_fits(df_miss, imp_vec = imp_vec, color_vec = color_vec)
p4 <- plot_county_fits(county_miss, imp_vec, color_vec, title = 'p = 0.5')

# get the metrics to evaluate each method
res4 <- calculate_metrics(df_miss, imp_vec,imputed_only = F)
res5 <- calculate_metrics(df_miss, imp_vec,imputed_only = T)
res6 <- calculate_metrics(county_miss, imp_vec,imputed_only = F, median_estimate = T)

#### MNAR ####

# simulation function!
df_miss = MNAR_sim(df, p = 0.2, direction = 'upper', alpha = 1, by_facility = T)

# run the periodic imputation
periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
df_miss = periodic_list$df

# run the CARBayes imputation
CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)

# get the facility and county level imputation results
df_miss = CAR_list$facility_df
county_miss = merge(periodic_list$county_fit, CAR_list$county_df, by = 'date')

# make the plots for the results
imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')
color_vec = c('red','blue')
p5 <- plot_facility_fits(df_miss, imp_vec = imp_vec, color_vec = color_vec)
p6 <- plot_county_fits(county_miss, imp_vec, color_vec, title = 'p = 0.2')

# ggsave(plot = p5, filename = 'C:/Users/nickl/Documents/global_covid19_response/figures/MNAR_p02_simulation_fits_05172021.png', width = 15, height = 10)
# ggsave(plot = p6, filename = 'C:/Users/nickl/Documents/global_covid19_response/figures/MNAR_p02_simulation_countyFit_05172021.png', width = 10, height = 7)


res7 <- calculate_metrics(df_miss, imp_vec,imputed_only = F)
res8 <- calculate_metrics(df_miss, imp_vec,imputed_only = T)
res9 <- calculate_metrics(county_miss, imp_vec,imputed_only = F, median_estimate = T)
