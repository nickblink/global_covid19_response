library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)
library(cowplot)

# Rscript R/real_data_main.R R_PI=2:CARburnin=10:CARnsample=20:output_path=results/real_data_analysis_allnegbin_TEST.RData
# Rscript R/real_data_main.R R_PI=200:CARburnin=5000:CARnsample=10000:output_path=results/real_data_analysis_rolling_07232024.RData
# sbatch -J real_anal run_sim_real_anal.sh R_PI=200:CARburnin=5000:CARnsample=10000:output_path=results/real_data_analysis_rolling_07232024.RData

source('R/imputation_functions.R')
rstan_options(auto_write = TRUE)

# register the cores
registerDoParallel(cores = 20)

inputs <- c('R_PI=2:CARburnin=20:CARnsample=40:output_path=results/real_data_analysis_07052024.RData\r')
inputs <- commandArgs(trailingOnly = TRUE)

params <- list()
inputs <- gsub('\r', '', inputs)
for(str in strsplit(inputs,':')[[1]]){
  tmp = strsplit(str, '=')[[1]]
  nn = tmp[1]
  val = tolower(tmp[2])
  if(nn %in% c('R_PI','CARburnin','CARnsample')){
    val = as.numeric(val)
  }
  params[[nn]] = val
}

### Data prep
{
  data_file <- ifelse(file.exists('data/liberia_cleaned_01-06-2021.rds'), 'data/liberia_cleaned_01-06-2021.rds', 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/liberia_cleaned_01-06-2021.rds')
  Dfull <- readRDS(data_file)
  
  Dfull %>% group_by(county) %>%
    summarize(n = length(unique(facility)))
  
  # most recent data on github.
  D <- Dfull %>%
    filter(county == 'Maryland') %>%
    select(date, district, facility, ari = indicator_count_ari_total, indicator_denom) %>%
    add_periodic_cov()
  
  D %>% 
    filter(date >= '2020-01-01') %>%
    group_by(date) %>% 
    summarize(n_miss = length(unique(is.na(ari))))
  
  D %>% 
    filter(date < '2020-01-01') %>%
    summarize(prop_miss = mean(is.na(ari)))
  
  # get the dates
  dates = unique(D$date)
  eval_dates = dates[dates >= '2020-01-01']
  
  # anonymize and make name matching fit the code.
  dist_n <- D %>% group_by(district) %>%
    summarize(n = length(unique(facility))) %>% 
    arrange(n)
  uni_district = unique(dist_n$district)
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
}
res_list <- list()
R_PI = params[['R_PI']] # 200
burnin = params[['CARburnin']] # 1000
nsample = params[['CARnsample']] # 2000

### Rolling surveillance
# get dates
all_dates = sort(unique(df$date))
eval_dates = all_dates[all_dates >= '2020-01-01']

# full res_list and prediction
full_res_list <- list()
df_predict <- NULL

one_run <- function(i, models = c('WF','WF_NB','freqGLM','freqGLM_NB','CAR_nsample','CAR_phifit','CAR_phifit_negbin')){
  res_list <- list()
  # get the dates
  # dates <- all_dates[all_dates >= (eval %m-% months(48)) & all_dates <= eval]
  d = eval_dates[i]
  ind <- which(all_dates == d)
  dates <- all_dates[(ind - 48):ind]
  train_end <- all_dates[ind-1]
  
  # get the data frame
  df_roll <- df %>% 
    filter(date %in% dates)
  
  # run WF model
  if('WF' %in% models){
    res_list[['WF']] <- WF_CCA(df_roll, col = "y", family = 'poisson', R_PI = R_PI, train_end_date = train_end)
    df_roll <- res_list[['WF']]$df
  }
  
  # run WF NB model
  if('WF_NB' %in% models){
    res_list[['WF_NB']] <- WF_CCA(df_roll, col = "y", family = 'negbin', R_PI = R_PI, train_end_date = train_end)
    df_roll <- res_list[['WF_NB']]$df 
  }
  
  # run freqGLM
  if('freqGLM' %in% models){
    system.time({
      res_list[['freqGLM']] <- freqGLMepi_CCA(df_roll, R_PI = R_PI, verbose = F, train_end_date = train_end)
    }) # 20m
    df_roll <- res_list[['freqGLM']]$df
  }
  
  
  # run freqGLM NB
  if('freqGLM_NB' %in% models){
    system.time({
      res_list[['freqGLM_NB']] <- freqGLMepi_CCA(df_roll, R_PI = R_PI, verbose = F, family = 'negbin', train_end_date = train_end)
    }) # 20m
    df_roll <- res_list[['freqGLM_NB']]$df
  }
  
  # run CAR
  if('CAR' %in% models){
    system.time({
      res_list[['CAR']] <- CARBayes_wrapper(df_roll, burnin = burnin, n.sample = nsample, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', model_rename = 'y_CAR_sample', train_end_date = train_end)
    }) # 143s
    df_roll <- res_list[['CAR']]$df
    # df_roll$y_CARstan <- df_roll$y_CARstan_0.5
    # colnames(df_roll) <- gsub('y_CARstan', 'y_CAR_sample', colnames(df_roll))
    # colnames(res_list[['CAR']]$district_df) <- gsub('y_CARstan', 'y_CAR_sample', colnames(res_list[['CAR']]$district_df))
    
  }

  if('CAR_phifit' %in% models){
    system.time({
      res_list[['CAR_phifit']] <- CARBayes_wrapper(df_roll, burnin = burnin, n.sample = nsample, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', train_end_date = train_end, model_rename = 'y_CAR_phifit', use_fitted_phi = T)
    })
    df_roll <- res_list[['CAR_phifit']]$df
    # df_roll$y_CARstan <- df_roll$y_CARstan_0.5
    # colnames(df_roll) <- gsub('y_CARstan', 'y_CAR_phifit', colnames(df_roll))
    # colnames(res_list[['CAR_phifit']]$district_df) <- gsub('y_CARstan', 'y_CAR_phifit', colnames(res_list[['CAR_phifit']]$district_df))
    
  }
  
  if('CAR_phifit_negbin' %in% models){
    res_list[['CAR_phifit_negbin']] <- CARBayes_wrapper(df_roll, burnin = burnin, n.sample = nsample, prediction_sample = F, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', train_end_date = train_end, use_fitted_phi = T, model_rename = 'y_CAR_phifit_negbin', family = 'negbin')
    df_roll <- res_list[['CAR_phifit_negbin']]$df
  }
  
  # update results
  df_predict <- rbind(df_predict, 
                      df_roll[df_roll$date == d,])
  full_res_list[[d]] <- df_roll
  
  return(list(df = df_roll, res_list = res_list))
}

# print('starting one run')
# test <- one_run(1)
# print('done')

system.time({
  results_list <- foreach(i = 1:length(eval_dates)) %dorng% one_run(i, models = c('WF_NB', 'freqGLM_NB', 'CAR_phifit_negbin'))
})
names(results_list) <- eval_dates

# if(file.exists('C:/Users/nickl/') | file.exists('C:/Users/Admin-Dell/')){
#   results_list <- NULL
#   for(i in 1:length(eval_dates)){
#     print(i)
#     results_list[[i]] <- one_run(i, models = 'CAR_phifit_negbin')
#   }
# }

save(results_list, file = params[['output_path']])