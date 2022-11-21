library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
#library(cowplot)
library(doRNG)
library(doParallel)

source('R/imputation_functions.R')

# register the cores
registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
inputs <- c('0.7','mnar', 'noST','1000','250','7')
inputs <- commandArgs(trailingOnly = TRUE)
print(inputs)

# pull parameters into proper format
p <- as.numeric(inputs[[1]])
missingness <- tolower(inputs[[2]])
DGP <- tolower(inputs[[3]])
R <- as.integer(inputs[[4]])
num_jobs <- as.integer(inputs[[5]])
job_id <- as.integer(inputs[[6]])

# check that proper missingness is input
if(!(missingness %in% c('mcar','mar','mnar'))){
  stop('please input a proper missingness')
}else{
  print(sprintf('proceeding with %s missingness', missingness))
}

# missingness parameters
params <- list()
if(missingness == 'mar'){
  params[['rho']] <- 0.7
  params[['alpha']] <- 0.7
  params[['tau']] <- 3
}else if(missingness =='mnar'){
  gamma = 1
  params[['gamma']] <- gamma
}

# get sequence of simulation iterations to run
if(job_id < num_jobs){
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):(floor(R/num_jobs)*job_id)
}else{
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):R
}

# set up the output folder
date <- gsub('-','_', Sys.Date())
output_path <- sprintf('results/%s%s_%s_%s', missingness, gsub('\\.','',p), DGP, date)
if(!file.exists(output_path)){
  dir.create(output_path, recursive = T)
}
results_file <- sprintf('%s/sim_results_p%1.1f_%s_%i(%i).RData', output_path, p, missingness, job_id, num_jobs)

if(file.exists(results_file)){
  iter = 0
  while(file.exists(results_file)){
    iter = iter + 1
    results_file <- sprintf('%s/sim_results_p%1.1f_%s_%i(%i)(%i).RData', output_path, p, missingness, job_id, num_jobs, iter)
  }
}

#### MCAR p = 0.2 no ST ####
if(DGP == 'nost'){
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
}else{
  stop('unrecognized data generating process')
}

print('data made')

if(job_id == 1){
  save(lst, file = paste0(output_path, '/simulated_data.RData'))
}

one_run <- function(lst, i){
  set.seed(1)
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
  
  print('running freqGLM_epi')
  # run the freqGLM_epi complete case analysis
  freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 200, verbose = F)
  df_miss <- freqGLMepi_list$df
  
  # run the WF baseline (complete data) imputation
  # df_miss <- WF_baseline(df_miss, R_PI = 100)
  
  print('running CCA')
  # run the WF complete case analysis model
  df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 200)
  
  print('running CARBayes')
  # run the CAR complete case analysis model
  df_miss <- CARBayes_CCA(df_miss, burnin = 10000, n.sample = 20000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')

  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  return(df_miss)
}

for(i in seq){
  one_run(lst, i)
}

# %dorng% works on the cluster. %do% works at home
system.time({
  imputed_list <- foreach(i=seq) %dorng% one_run(lst, i)
})

save(imputed_list, seq, params, file = results_file)

