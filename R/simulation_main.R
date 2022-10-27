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
registerDoParallel(cores = 10)

# get the parameters (first line is for testing on my home computer)
params <- c('0.2','mcar','20','2','1')
params <- commandArgs(trailingOnly = TRUE)

# pull parameters into proper format
p <- as.numeric(params[[1]])
missingness <- tolower(params[[2]])
R <- as.integer(params[[3]])
num_jobs <- as.integer(params[[4]])
job_id <- as.integer(params[[5]])

# check that proper missingness is input
if(!(missingness %in% c('mcar','mar','mnar'))){
  stop('please input a proper missingness')
}else{
  print(sprintf('proceeding with %s missingness', missingness))
}

# missingness parameters
if(missingness == 'mar'){
  rho = 0.7
  alpha = 0.7
  tau = 3
}else if(missingness =='mnar'){
  gamma = 1
}

# get sequence of simulation iterations to run
if(job_id < num_jobs){
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):(floor(R/num_jobs)*job_id)
}else{
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):R
}

# set up the output folder
date <- gsub('-','_', Sys.Date())
output_path <- sprintf('results/%s_%s', missingness, date)
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
lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')

print('data made')

if(job_id == 1){
  save(lst, file = paste0(output_path, '/simulated_data.RData'))
}

one_run <- function(lst, i){
  print(length(lst))
  print(length(lst$df_list))
  print(i)
  df = lst$df_list[[i]]
  
  # add in missingness
  if(missingness == 'mcar'){
    df_miss = MCAR_sim(df, p = p, by_facility = T)
  }else if(missingness == 'mar'){
    df_miss = MAR_spatiotemporal_sim(df, p = p, rho = rho, alpha = alpha, tau = tau)
  }else{
    f_miss <- MNAR_sim(df, p = p, direction = 'upper', gamma = gamma, by_facility = T)
  }
  
  # run the WF baseline (complete data) imputation
  df_miss <- WF_baseline(df_miss, R_PI = 100)

  return(df_miss)
}



# %dorng% works on the cluster. %do% works at home
system.time({
  imputed_list <- foreach(i=seq) %dorng% one_run(lst, i)
})

if(missingness == 'mar'){
  save(imputed_list, rho, alpha, tau, file = results_file)
}else if(missingness == 'mnar'){
  save(imputed_list, gamma, file = results_file)
}else{
  save(imputed_list, file = results_file)
}



