library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)
library(cowplot)

source('R/imputation_functions.R')

if(file.exists('C:/Users/Admin-Dell')){
  res_dir = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}else{
  res_dir = "C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}

# most recent data on github.
D = readRDS(sprintf('%s/data/liberia_cleaned_01-06-2021.rds', res_dir)) %>%
  filter(county == 'Maryland') %>%
  select(date, district, facility, ari = indicator_count_ari_total, indicator_denom) %>%
  add_periodic_cov()

# get the dates
dates = unique(D$date)
eval_dates = dates[dates >= '2020-01-01']

# anonymize and make name matching fit the code.
uni_district = unique(D$district)
matid = match(D$district, uni_district)
D$district = toupper(letters[matid])

# rename the facilities
df = NULL
for(d in unique(D$district)){
  tmp = D %>% filter(district == d)
  uni_fac = unique(tmp$facility)
  matid = match(tmp$facility, uni_fac)
  tmp$facility = paste0(d, matid)
  df = rbind(df, tmp)
}

df$y <- df$ari

res_list <- list()

# run WF model
res_list[['WF']] <- WF_CCA(df, col = "y", family = 'poisson', R_PI = 200)
df <- res_list[['WF']]$df

# run freqGLM
res_list[['freqGLM']] <- freqGLMepi_CCA(df, R_PI = 10, verbose = F)
df <- res_list[['freqGLM']]$df

# run CAR
res_list[['CAR']] <- CARBayes_wrapper(df, burnin = 100, n.sample = 200, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
df <- res_list[['CAR']]$df
df$y_CARstan <- df$y_CARstan_0.5

# merge all the results together
district_df <- merge(res_list[['WF']]$district_df, res_list[['freqGLM']]$district_df, by= c('district', 'date')) %>%
  merge(res_list[['CAR']]$district_df, by = c('district', 'date'))

### Plot 1: Plotting the facility fits
tmp  = df %>% filter(district == 'A')
plot_facility_fits(tmp, methods = c('y_pred_WF','y_pred_freqGLMepi', 'y_CARstan'))

### Plot 2: Plotting the proportion outbreaks over time.



