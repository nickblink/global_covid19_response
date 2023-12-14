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

### functions
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

compare_outcomes <- function(input_lst, names){
  lst <- get_data(input_lst[[1]])
  outcomes <- lst$df_list[[1]] %>% group_by(date) %>%
    summarize(mean_y = mean(y),
              median_y = median(y))
  outcomes$DGP = names[1]
  
  full_outcomes <- outcomes
  
  for(i in 2:length(input_lst)){
    lst <- get_data(input_lst[[i]])
    outcomes <- lst$df_list[[1]] %>% group_by(date) %>%
      summarize(mean_y = mean(y),
                median_y = median(y))
    outcomes$DGP = names[i]
    
    full_outcomes <- rbind(full_outcomes, outcomes)
  }
  
  p1 <- ggplot(data = full_outcomes, aes(x = date, y = mean_y, color = DGP)) + 
    geom_line()
  
  p2 <- ggplot(data = full_outcomes, aes(x = date, y = median_y, color = DGP)) + 
    geom_line()
  
  cowplot::plot_grid(p1, p2, nrow = 2)
}

### Comparisons same betas
compare_outcomes(list(c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=CAR:R=20:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST','3'), 
                      c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=freqglm:R=20:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:output_path=TEST','3'),
                      inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=WF:R=20:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05','3')), 
                 names = c('CAR','freqGLM','WF'))

### Comparisons match b0 = 6, b1 = -0.25
compare_outcomes(list(c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=CAR:R=20:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST','3'), 
                      c('p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=freqglm:R=20:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:output_path=TEST','3'),
                      inputs <- c('p=0.1:b0_mean=6:b1_mean=n0.25:missingness=mcar:DGP=WF:R=20:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05','3')), 
                 names = c('CAR','freqGLM','WF'))
# ok, but freqGLM is actually a bit lower for EB0 = 5.

# EB0 = 5.5 is almost exactly right here

### Comparisons match b0 = 6, b1 = 0
compare_outcomes(list(c('p=0.1:b0_mean=6:b1_mean=0:missingness=mcar:DGP=CAR:R=20:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST','3'), 
                      c('p=0.1:b0_mean=5.5:b1_mean=0:missingness=mcar:DGP=freqglm:R=20:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:output_path=TEST','3'),
                      inputs <- c('p=0.1:b0_mean=6:b1_mean=0:missingness=mcar:DGP=WF:R=20:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05','3')), 
                 names = c('CAR','freqGLM','WF'))
# same. Got it.

### Comparisons match b0 = 2, b1 = 0
compare_outcomes(list(c('p=0.1:b0_mean=2:b1_mean=0:missingness=mcar:DGP=CAR:R=20:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST','3'), 
                      c('p=0.1:b0_mean=1.5:b1_mean=0:missingness=mcar:DGP=freqglm:R=20:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:output_path=TEST','3'),
                      inputs <- c('p=0.1:b0_mean=2:b1_mean=0:missingness=mcar:DGP=WF:R=20:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05','3')), 
                 names = c('CAR','freqGLM','WF'))
# Yep perfect. Got it.


### Do things change with QP? They shouldn't
compare_outcomes(list(c('p=0.1:b0_mean=6:b1_mean=0:missingness=mcar:DGP=CAR:R=20:num_jobs=20:rho_DGP=0.3:alpha_DGP=0.3:tau2_DGP=0.25:output_path=TEST:family=quasipoisson:theta=4','3'), 
                      c('p=0.1:b0_mean=5.5:b1_mean=0:missingness=mcar:DGP=freqglm:R=20:num_jobs=20:rho_DGP=0.2:alpha_DGP=0.2:output_path=TEST:family=quasipoisson:theta=4','3'),
                      inputs <- c('p=0.1:b0_mean=6:b1_mean=0:missingness=mcar:DGP=WF:R=20:num_jobs=20:output_path=mcar01_WF_QPtheta9_beta06_beta1n025_ID499135_2023_12_05:family=quasipoisson:theta=4','3')), 
                 names = c('CAR','freqGLM','WF'))
# ya same. Of course. Not sure why the median values are decreasing here, but could just be randomness.