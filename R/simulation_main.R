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

## Interactive job command
# Rscript R/simulation_main.R p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=freqGLM:rho_DGP=0.2:alpha_DGP=0.2:R=50:num_jobs=50:output_path=TEST_2024_07_05:family=poisson:empirical_betas=F:CARburnin=100:CARnsample=200:models=5,6:R_PI=10 3

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id
inputs <- c('p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=freqGLM:rho_DGP=0.2:alpha_DGP=0.2:R=1000:num_jobs=50:output_path=mcar01_WF_QPtheta9_beta055_beta1n025_ID499135_2023_12_05:family=quasipoisson:theta=5:empirical_betas=F:CARburnin=2000:CARnsample=4000:R_PI=200:models=1,2,3,4,5,6\r','3')
inputs <- commandArgs(trailingOnly = TRUE)
print(inputs)

### Get parameters (for home and cluster)
{
# pull parameters into proper format
params <- list()

inputs[[1]] <- gsub('\r', '', inputs[[1]])
params[['job_id']] <- as.integer(inputs[[2]])
for(str in strsplit(inputs[[1]],':')[[1]]){
  tmp = strsplit(str, '=')[[1]]
  nn = tmp[1]
  val = tolower(tmp[2])
  if(nn %in% c('p','rho_DGP','alpha_DGP','tau2_DGP','rho_MAR','alpha_MAR','tau2_MAR','gamma','theta', 'dispersion','CARburnin','CARnsample','R_PI')){
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
  }else if(val %in% c('t','f')){
    val = as.logical(ifelse(val == 't', T, F))
  }else if(nn == 'models'){
    val <- as.numeric(strsplit(val, ',')[[1]])
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

# setting the CAR burnin and nsample params if not provided.
if(is.null(params[['CARburnin']])){
  params[['CARburnin']] <- 1000
}
if(is.null(params[['CARnsample']])){
  params[['CARnsample']] <- 2000
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
    if('theta' %in% names(params)){
      arguments = c(arguments,
                    list(family = 'quasipoisson',
                         theta = params[['theta']]))
    }else{
      stop('need to put in a theta value for quasipoisson')
    }
  }else if(params[['family']] == 'negbin'){
    arguments = c(arguments,
                  list(family = 'negbin'))
    if('dispersion' %in% names(params)){
      arguments = c(arguments,
                    list(dispersion = params[['dispersion']]))
    }
  }else if(params[['family']] != c('poisson')){
    stop('improper family for DGP')
  }
}

if('empirical_betas' %in% names(params)){
  arguments <- c(arguments, 
                 list(empirical_betas = params[['empirical_betas']]))
}
}

# Simulate the data
lst <- do.call(simulate_data, arguments)

print('data made')

# initialize the error catcher for each run
# errors <- list(freqEpi = data.frame(i = NULL, error = NULL),
#                WF = data.frame(i = NULL, error = NULL),
#                CARBayesST = data.frame(i = NULL, error = NULL),
#                CARstan = data.frame(i = NULL, error = NULL))



### File saving (for cluster only)
{
# set up the output folder
date <- gsub('-','_', Sys.Date())
if(params[['output_path']] == 'na'){
  params[['output_path']] <- sprintf('results/%s%s_%s_%s', params[['missingness']], gsub('\\.','', params[['p']]), params[['DGP']], date)
}else{
  params[['output_path']] <- sprintf('results/%s', params[['output_path']])
}

if(!file.exists(params[['output_path']])){
  dir.create(params[['output_path']], recursive = T)
}
results_file <- sprintf('%s/sim_results_p%1.1f_%s_%i(%i).RData', params[['output_path']], params[['p']], params[['missingness']], params[['job_id']], params[['num_jobs']])

if(file.exists(results_file)){
  iter = 0
  while(file.exists(results_file)){
    iter = iter + 1
    results_file <- sprintf('%s/sim_results_p%1.1f_%s_%i(%i)(%i).RData', params[['output_path']], params[['p']], params[['missingness']], params[['job_id']], params[['num_jobs']], iter)
  }
}

# save the data
if(params[['job_id']] == 1){
  save(lst, file = paste0(params[['output_path']], '/simulated_data.RData'))
}

}

one_run <- function(lst, i, model_list = NULL){
  
  t0 <- Sys.time()
  print(sprintf('i = %i',i))
  df = lst$df_list[[i]]
  
  set.seed(i)
  
  # add in missingness to the data
  if(params[['missingness']] == 'mcar'){
    df_miss = MCAR_sim(df, p = params[['p']], by_facility = T)
  }else if(params[['missingness']] == 'mar'){
    df_miss = MAR_spatiotemporal_sim(df, p = params[['p']], rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau2 = params[['tau2_MAR']])
  }else{
    df_miss <- MNAR_sim(df, p = params[['p']], direction = 'upper', gamma = params[['gamma']], by_facility = T)
  }
  
  # initialize errors 
  errors <- list()
  for(j in 1:length(model_list)){
    errors[[model_list[[j]]$model]] <- data.frame(i = NULL, error = NULL)
  }
  
  # initializing the return list
  return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
  
  # cycle through models
  for(sub_lst in model_list){
    t1 = Sys.time()
    model = sub_lst$model
    print(sprintf('------(%s)%s-------', i, model))
    params = sub_lst$params
    if(model == 'WF'){
      fit_fxn <- function(df){
        tmp <- WF_CCA(df, col = "y", family = 'poisson', R_PI = params[['R_PI']])
        return(tmp)
      }
    }else if(model == 'WF_NB'){
      fit_fxn <- function(df){
        tmp <- WF_CCA(df, col = "y", family = 'negbin', R_PI = params[['R_PI']])
        return(tmp)
      }
    }else if(model == 'freqGLM'){
      fit_fxn <- function(df){
        tmp <- freqGLMepi_CCA(df, R_PI = params[['R_PI']], verbose = F) 
        return(tmp)
      }
    }else if(model == 'freqGLM_NB'){
      fit_fxn <- function(df){
        tmp <- freqGLMepi_CCA(df, R_PI = params[['R_PI']], verbose = F, family = 'negbin') 
        return(tmp)
      }
    }else if(model == 'CAR_sample'){
      fit_fxn <- function(df){
        tmp <- CARBayes_wrapper(df, burnin = params[['burnin']], n.sample = params[['n.sample']], prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', 
                                model_rename = 'y_CAR_sample')
        return(tmp)
      }
    }else if(model == 'CAR_phifit'){
      fit_fxn <- function(df){
        tmp <- CARBayes_wrapper(df, burnin = params[['burnin']], n.sample = params[['n.sample']], prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', use_fitted_phi = T, 
                                model_rename = 'y_CAR_phifit')
        return(tmp)
      }
    }
    
    # run the model
    res <- tryCatch({
      fit_fxn(return_list[['df_miss']])
    }, error = function(e){
      e[[1]]
    })
    
    # store the results
    if(class(res) == 'list'){
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      if(model == 'WF'){
        return_list[['WF_betas']] <- res$betas
      }else if(model == 'WF_NB'){
        return_list[['WF_NB_betas']] <- res$betas
        return_list[['WF_NB_overdisp']] <- res$overdisp
      }else if(model == 'freqGLM'){
        return_list[['freqGLM_params']] <- res$param_results
      }else if(model == 'freqGLM_NB'){
        return_list[['freqGLM_NB_params']] <- res$param_results
      }else if(model == 'CAR_sample'){
        return_list[['CAR_sample_summary']] <- res$CARstan_summary
      }else if(model == 'CAR_phifit'){
        return_list[['CAR_phifit_summary']] <- res$CARstan_summary
      }
    }else{
      return_list[['errors']][[model]] <- rbind(return_list[['errors']][[model]], data.frame(i = i, error = res[[1]]))
    }
    
    # keep track of timing
    return_list[['timing']][[model]] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  return(return_list)
}

model_list <- list(list(model = 'WF',
                        params = list(R_PI = params[['R_PI']])),
                   list(model = 'WF_NB',
                        params = list(R_PI = params[['R_PI']])),
                   list(model = 'freqGLM',
                        params = list(R_PI = params[['R_PI']])),
                   list(model = 'freqGLM_NB',
                        params = list(R_PI = params[['R_PI']])),
                   list(model = 'CAR_sample',
                        params = list(burnin = params[['CARburnin']], n.sample = params[['CARnsample']])),
                   list(model = 'CAR_phifit',
                        params = list(burnin = params[['CARburnin']], n.sample = params[['CARnsample']])))

model_list <- model_list[params[['models']]]

# run the models for each simulation dataset
system.time({
  imputed_list <- foreach(i=1:R_new) %dorng% one_run(lst, i, model_list)
})

true_betas <- lst$betas

save(imputed_list, seq, params, arguments, true_betas, file = results_file)