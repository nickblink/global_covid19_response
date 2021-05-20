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


#### Assessing the simulated data ####

lst <- simulate_data(district_sizes = c(2,3,4,5,7), n = 2)
df = lst$df

lst <- simulate_data_spatiotemporal(district_sizes = c(4), n = 2, rho = 0.3, alpha = 0.5, tau = 0.5)
df = lst$df_list[[1]]

lst2 <- simulate_data(district_sizes = c(4), n = 2)
df2 = lst2$df_list[[1]]

par(mfrow = c(2,2))
for(f in sample(unique(df$facility), 4)){
  tmp = df %>% filter(facility == f)
  tmp2 = df2 %>% filter(facility == f)
  
  plot(tmp$date, tmp$y, type = 'l', main = f, ylim = c(0, 1.2*max(tmp$y)))
  lines(tmp2$date, tmp2$y, col = 'red')
}

#### MCAR p = 0.2 spatio-temporal
lst <- simulate_data_spatiotemporal(district_sizes = c(4), n = 3, rho = 0.3, alpha = 0.5, tau = 0.5)

imp_vec = c('y_pred_harmonic', 'y_CARBayes_ST')

for(i in 1:3){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
  
  # run the periodic imputation
  periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
  df_miss = periodic_list$df
  
  # run the CARBayes imputation
  CAR_list = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T)
  df_miss = CAR_list$facility_df
}

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
