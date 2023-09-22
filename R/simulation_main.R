library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)

source('R/imputation_functions.R')

# register the cores
registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id
inputs <- c('0.1', '6', 'n0.25', 'mcar', 'CAR', '0.3', '0.3', '1', '10','5','test','1')
inputs <- commandArgs(trailingOnly = TRUE)
print(inputs)

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
rho_DGP <- as.numeric(tolower(inputs[[6]]))
alpha_DGP <- as.numeric(tolower(inputs[[7]]))
tau2_DGP <- as.numeric(tolower(inputs[[8]]))
R <- as.integer(inputs[[9]])
num_jobs <- as.integer(inputs[[10]])
output_path <- tolower(inputs[[11]])
job_id <- as.integer(inputs[[12]])

# check that proper missingness is input
if(!(missingness %in% c('mcar','mar','mnar'))){
  stop('please input a proper missingness')
}else{
  print(sprintf('proceeding with %s missingness', missingness))
}

# storing all the parameters
params <- list()
params[['p']] <- p
params[['missingness']] <- missingness
params[['DGP']] <- DGP
if(DGP == 'car'){
  params[['rho_DGP']] <- rho_DGP
  params[['alpha_DGP']] <- alpha_DGP
  params[['tau2_DGP']] <- tau2_DGP
}
params[['b0_mean']] <- b0_mean
params[['b1_mean']] <- b1_mean
if(missingness == 'mar'){
  params[['rho_MAR']] <- 0.7
  params[['alpha_MAR']] <- 0.7
  params[['tau2_MAR']] <- 9
}else if(missingness =='mnar'){
  gamma = 1
  params[['gamma']] <- gamma
}
params[['R']] <- R
params[['num_jobs']] <- num_jobs
params[['job_id']] <- job_id

# get sequence of simulation iterations to run
if(job_id < num_jobs){
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):(floor(R/num_jobs)*job_id)
}else{
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):R
}

# set up the output folder
date <- gsub('-','_', Sys.Date())
if(output_path == 'na'){
  output_path <- sprintf('results/%s%s_%s_%s', missingness, gsub('\\.','',p), DGP, date)
}else{
  output_path <- sprintf('results/%s', output_path)
}

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

# Simulate the data
if(DGP == 'nRost'){
  lst <- simulate_data(district_sizes = c(4, 6, 10), 
                       R = R, 
                       end_date = '2020-12-01', 
                       b0_mean = b0_mean, 
                       b1_mean = b1_mean)
}else if(DGP == 'car'){
  lst <- simulate_data(district_sizes = c(4, 6, 10),
                       R = R, 
                       end_date = '2020-12-01',
                       b0_mean = b0_mean, 
                       b1_mean = b1_mean,
                       type = 'CAR',
                       rho = rho_DGP, 
                       alpha = alpha_DGP, 
                       tau2 = tau2_DGP)
}else{
  stop('unrecognized data generating process')
}

print('data made')

# save the data
if(job_id == 1){
  save(lst, file = paste0(output_path, '/simulated_data.RData'))
}

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
    df_miss = MAR_spatiotemporal_sim(df, p = p, rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau = params[['tau2_MAR']])
  }else{
    df_miss <- MNAR_sim(df, p = p, direction = 'upper', gamma = gamma, by_facility = T)
  }
  
  # initializing the return list
  return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, WF_betas = NULL, CAR_summary = NULL)
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
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
        res$district_df
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
      res <- CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
      #res <- CARBayes_wrapper(return_list[['df_miss']], burnin = 10000, n.sample = 20000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      return_list[['CAR_summary']] <- res$summary_stats
      return_list
    }, error = function(e){
      print(e)
      return_list[['errors']][['CARBayes']] <- rbind(return_list[['errors']][['CARBayes']], data.frame(i, error = e[[1]]))
      return_list
    })
  }

  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  # check I got the correct result names
  outcome_name_checker(return_list, models = models)
  
  return(return_list)
}

set.seed(1)

# run the models for each simulation dataset
system.time({
  imputed_list <- foreach(i=seq) %dorng% one_run(lst, i)
})

res <- one_run(lst, 1, models = 'CAR')

true_betas <- lst$betas

save(imputed_list, seq, params, true_betas, file = results_file)

