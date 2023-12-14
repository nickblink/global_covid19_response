library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)

source('R/imputation_functions.R')
rstan_options(auto_write = TRUE)

# register the cores
registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id

### Get parameters (for home and cluster)
get_data <- function(inputs){
  # pull parameters into proper format
  params <- list()
  
  inputs[[1]] <- gsub('\r', '', inputs[[1]])
  params[['job_id']] <- as.integer(inputs[[2]])
  for(str in strsplit(inputs[[1]],':')[[1]]){
    tmp = strsplit(str, '=')[[1]]
    nn = tmp[1]
    val = tolower(tmp[2])
    if(nn %in% c('p','rho_DGP','alpha_DGP','tau2_DGP','rho_MAR','alpha_MAR','tau2_MAR','gamma','theta')){
      val = as.numeric(val)
    }else if(nn == 'b0_mean'){
      val = as.numeric(strsplit(val, '/')[[1]])
    }else if(nn == 'b1_mean'){
      val = tryCatch({as.numeric(val)
      }, warning = function(w){
        if(substr(val,1,1) == 'n'){
          return(-as.numeric(substr(val, 2, nchar(val))))
        }
      })
    }else if(nn %in% c('R','num_jobs')){
      val = as.integer(val)
    }
    params[[nn]] = val
  }
  
  print(params)
  
  # check that proper missingness is input
  if(!(params[['missingness']] %in% c('mcar','mar','mnar'))){
    stop('please input a proper missingness')
  }else{
    print(sprintf('proceeding with %s missingness', params[['missingness']]))
  }
  
  # get sequence of simulation iterations to run
  # (deprecated - now I just simulate R_new # of data frames)
  if(params[['job_id']] < params[['num_jobs']]){
    seq <- (floor(params[['R']]/params[['num_jobs']])*(params[['job_id']] - 1) + 1):(floor(params[['R']]/params[['num_jobs']])*params[['job_id']])
  }else{
    seq <- (floor(params[['R']]/params[['num_jobs']])*(params[['job_id']] - 1) + 1):params[['R']]
  }
  
  R_new = length(seq)
  
  # input arguments to data simulation
  arguments = list(district_sizes = c(4, 6, 10), 
                   R = R_new, 
                   seed = params[['job_id']],
                   end_date = '2020-12-01', 
                   b0_mean = params[['b0_mean']], 
                   b1_mean = params[['b1_mean']])
  
  if(params[['DGP']] == 'car'){
    arguments = c(arguments, 
                     list(type = 'CAR',
                     rho = params[['rho_DGP']], 
                     alpha = params[['alpha_DGP']], 
                     tau2 = params[['tau2_DGP']]))
  }else if(params[['DGP']] == 'freqglm'){
    arguments = c(arguments,
                  list(type = 'freqGLM',
                       rho = params[['rho_DGP']],  
                       alpha = params[['alpha_DGP']]))
  }
  
  if(!is.null(params[['family']])){
    if(params[['family']] == 'quasipoisson'){
      arguments = c(arguments,
                    list(family = 'quasipoisson',
                         theta = params[['theta']]))
    }else if(params[['family']] != 'poisson'){
      stop('improper family for DGP')
    }
  }
  
  
  # Simulate the data
  lst <- do.call(simulate_data, arguments)
  print(length(lst[[1]]))
  
  print('data made')
  
  # initialize the error catcher for each run
  errors <- list(freqEpi = data.frame(i = NULL, error = NULL),
                 WF = data.frame(i = NULL, error = NULL),
                 CARBayesST = data.frame(i = NULL, error = NULL),
                 CARstan = data.frame(i = NULL, error = NULL))
  
  return(lst)
}

compute_variance <- function(lst, iters = NULL){
  if(is.null(iters)){
    iters = length(lst[[1]])
  }
  sd_df <- c()
  for(i in 1:iters){
    df = lst$df_list[[i]]
    z = (df$y - df$y_exp)/sqrt(df$y_var)
    sd_df <- c(sd_df, sd(z))
    print(sprintf('sd df: %s', sd(z)))
  }
  print(sprintf('mean(sd df): %s', mean(sd_df)))
  
  sd_district <- c()
  for(i in 1:iters){
    df = lst$district_list[[i]]
    z = (df$y - df$y_exp)/sqrt(df$y_var)
    sd_district <- c(sd_district, sd(z))
    print(sprintf('sd distrct: %s', sd(z)))
  }
  print(sprintf('mean(sd district): %s', mean(sd_district)))
}

### WF Poisson 
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=WF:R=100:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05','3')
lst <- get_data(inputs)
compute_variance(lst)
# good.

### WF QP-4
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=WF:R=100:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05:theta=4:family=quasipoisson\r','3')
lst <- get_data(inputs)
compute_variance(lst)
# good.

### WF QP-9
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=WF:R=100:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05:theta=9:family=quasipoisson\r','3')
lst <- get_data(inputs)
compute_variance(lst)
# good.

### freqGLM Poisson
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=freqglm:R=100:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:TEST','3')
lst <- get_data(inputs)
compute_variance(lst)
# good

### freqGLM QP 4
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=freqglm:R=100:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:TEST:theta=4:family=quasipoisson','3')
lst <- get_data(inputs)
compute_variance(lst)
# ok a little higher than 1, but I'll take it.

### CAR Poisson 
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=CAR:R=200:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST','3')
lst <- get_data(inputs)
compute_variance(lst)
# good enough

### CAR QP 4
inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=CAR:R=500:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST:family=quasipoisson:theta=4','3')
lst <- get_data(inputs)
compute_variance(lst)
# good enough. I'll take that
