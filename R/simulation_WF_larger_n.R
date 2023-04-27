library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
#library(cowplot)
library(doRNG)
library(doParallel)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')

# register the cores
registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST R #jobs name_output job_id
inputs <- c('0.1', '6', 'n0.1', 'mcar', 'noST','200','50','test','25')

# pull parameters into proper format
p <- as.numeric(inputs[[1]])
b0_mean <- as.numeric(strsplit(inputs[[2]], '/')[[1]])
b1_mean <- tryCatch({as.numeric(inputs[[3]])
}, warning = function(w){
  if(substr(inputs[[3]],1,1) == 'n'){
    return(-as.numeric(substr(inputs[[3]], 2, nchar(inputs[[3]]))))
  }
})
missingness <- tolower(inputs[[4]])
DGP <- tolower(inputs[[5]])
R <- as.integer(inputs[[6]])
num_jobs <- as.integer(inputs[[7]])
output_path <- tolower(inputs[[8]])
job_id <- as.integer(inputs[[9]])


# missingness parameters
params <- list()
params[['p']] <- p
params[['missingness']] <- missingness
params[['DGP']] <- DGP
params[['b0_mean']] <- b0_mean
params[['b1_mean']] <- b1_mean
if(missingness == 'mar'){
  params[['rho']] <- 0.7
  params[['alpha']] <- 0.7
  params[['tau']] <- 3
}else if(missingness =='mnar'){
  gamma = 1
  params[['gamma']] <- gamma
}

seq <- 1:R

# initialize the error catcher for each run
errors <- list(freqEpi = data.frame(i = NULL, error = NULL),
                 WF = data.frame(i = NULL, error = NULL),
                 CARBayes = data.frame(i = NULL, error = NULL))

# function to run all models for a specific dataset
one_run <- function(lst, i, models = c('freq', 'WF', 'CAR')){
  t0 <- Sys.time()
  print(sprintf('i = %i',i))
  df = lst$df_list[[i]]
  
   # add in missingness
  if(missingness == 'mcar'){
    df_miss = MCAR_sim(df, p = p, by_facility = T)
  }else if(missingness == 'mar'){
    df_miss = MAR_spatiotemporal_sim(df, p = p, rho = params[['rho']], alpha = params[['alpha']], tau = params[['tau']])
  }else{
    df_miss <- MNAR_sim(df, p = p, direction = 'upper', gamma = gamma, by_facility = T)
  }
  
  # initializing the return list
  return_list <- list(df_miss = df_miss, errors = errors, WF_betas = NULL)
  rm(df_miss)
  
  # run the freqGLM_epi complete case analysis
  if('freq' %in% models){
    print('running freqGLM_epi')
    return_list <- tryCatch({
      freqGLMepi_list = freqGLMepi_CCA(return_list[['df_miss']], R_PI = 200, verbose = F)
      return_list[['df_miss']] <- freqGLMepi_list$df
      return_list
    }, error = function(e){
      return_list[['errors']][['freqEpi']] <- rbind(return_list[['errors']][['freqEpi']], data.frame(i, error = e[[1]]))
      return_list
    })
  }
  
  
  # run the WF complete case analysis model
  if('WF' %in% models){
    print('running WF CCA')
    return_list <- tryCatch({
      res <- WF_CCA(return_list[['df_miss']], col = "y", family = 'poisson', R_PI = 200)
      return_list[['df_miss']] <- res$df
      return_list[['WF_betas']] <- res$betas
      return_list
    }, error = function(e){
      return_list[['errors']][['WF']] <- rbind(return_list[['errors']][['WF']], data.frame(i = i, error = e[[1]]))
      return_list
    })
  }
  
  
  # run the CAR complete case analysis model
  if('CAR' %in% models){
    print('running CARBayes')
    return_list <- tryCatch({
      return_list[['df_miss']] <- CARBayes_CCA(return_list[['df_miss']], burnin = 10000, n.sample = 20000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
      return_list
    }, error = function(e){
      print(e)
      return_list[['errors']][['CARBayes']] <- rbind(return_list[['errors']][['CARBayes']], data.frame(i, error = e[[1]]))
      return_list
    })
  }
  

  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  return(return_list)
}

#### Run this for 4 year, 8 years, and 12 years

## 4 years
# Simulate the data
if(DGP == 'nost'){
  lst4 <- simulate_data(district_sizes = c(4, 6, 10), R = R, start_date = '2016-01-01', end_date = '2020-12-01', b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = 0.1)
}else{
  stop('unrecognized data generating process')
}

system.time({
  imputed_list4 <- foreach(i=seq) %do% one_run(lst4, i, models = 'WF')
})

true_betas4 <- lst4$betas

## 8 years
# Simulate the data
if(DGP == 'nost'){
  lst8 <- simulate_data(district_sizes = c(4, 6, 10), R = R, start_date = '2012-01-01', end_date = '2020-12-01', b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = 0.1)
}else{
  stop('unrecognized data generating process')
}

system.time({
  imputed_list8 <- foreach(i=seq) %do% one_run(lst8, i, models = 'WF')
})

true_betas8 <- lst8$betas

## 12 years
# Simulate the data
if(DGP == 'nost'){
  lst12 <- simulate_data(district_sizes = c(4, 6, 10), R = R, start_date = '2008-01-01', end_date = '2020-12-01', b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = 0.1)
}else{
  stop('unrecognized data generating process')
}

system.time({
  imputed_list12 <- foreach(i=seq) %do% one_run(lst12, i, models = 'WF')
})

true_betas12 <- lst12$betas

save(imputed_list4, imputed_list8, imputed_list12, true_betas4, true_betas8, true_betas12, file = 'data/WF_n_comparisons_04162023.RData')

## 20 years
# Simulate the data
if(DGP == 'nost'){
  lst20 <- simulate_data(district_sizes = c(4, 6, 10), R = R, start_date = '2000-01-01', end_date = '2020-12-01', b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = 0.1)
}else{
  stop('unrecognized data generating process')
}

system.time({
  imputed_list20 <- foreach(i=seq) %do% one_run(lst20, i, models = 'WF')
})

true_betas20 <- lst20$betas

save(imputed_list4, imputed_list8, imputed_list12, imputed_list20, true_betas4, true_betas8, true_betas20, file = 'data/WF_n_comparisons_04272023.RData')
