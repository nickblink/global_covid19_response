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
inputs <- c('p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=freqGLM:rho_DGP=0.2:alpha_DGP=0.2:R=1000:num_jobs=50:output_path=mcar01_WF_QPtheta9_beta055_beta1n025_ID499135_2023_12_05:family=quasipoisson:theta=5:empirical_betas=F:CARburnin=2000:CARnsample=4000\r','3')
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
  if(nn %in% c('p','rho_DGP','alpha_DGP','tau2_DGP','rho_MAR','alpha_MAR','tau2_MAR','gamma','theta', 'dispersion','CARburnin','CARnsample')){
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

print('Running on 24 facilities!!')

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

{
# # function to run all models for a specific dataset
# one_run <- function(lst, i, models = c('freq', 'WF', 'CARBayesST','CARstan'), WF_params = list(R_PI = 200, family = 'poisson'), freqGLM_params = list(R_PI = 200), MCMC_params = list(burnin.stan = 1000, n.sample.stan = 2000, burnin.CARBayesST = 5000, n.sample.CARBayesST = 10000)){
#   
#   
#   t0 <- Sys.time()
#   print(sprintf('i = %i',i))
#   df = lst$df_list[[i]]
#   
#   set.seed(i)
#   
#    # add in missingness to the data
#   if(params[['missingness']] == 'mcar'){
#     df_miss = MCAR_sim(df, p = params[['p']], by_facility = T)
#   }else if(params[['missingness']] == 'mar'){
#     df_miss = MAR_spatiotemporal_sim(df, p = params[['p']], rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau2 = params[['tau2_MAR']])
#   }else{
#     df_miss <- MNAR_sim(df, p = params[['p']], direction = 'upper', gamma = params[['gamma']], by_facility = T)
#   }
#   
#   # initializing the return list
#   return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
#   rm(df_miss)
#   
#   set.seed(i)
#   # run the freqGLM_epi complete case analysis
#   if('freq' %in% models){
#     t1 = Sys.time()
#     print('running freqGLM_epi')
#     return_list <- tryCatch({
#       freqGLMepi_list = freqGLMepi_CCA(return_list[['df_miss']], R_PI = freqGLM_params[['R_PI']], verbose = F)
#       return_list[['df_miss']] <- freqGLMepi_list$df
#       return_list[['district_df']] <- merge(return_list[['district_df']], freqGLMepi_list$district_df, by = c('district','date'))
#       return_list[['freqGLM_params']] <- freqGLMepi_list$param_results
#       return_list
#     }, error = function(e){
#       return_list[['errors']][['freqEpi']] <- rbind(return_list[['errors']][['freqEpi']], data.frame(i, error = e[[1]]))
#       return_list
#     })
#     return_list[['timing']][['freq']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
#   }
#   
#   set.seed(i)
#   # run the WF complete case analysis model
#   if('WF' %in% models){
#     t1 = Sys.time()
#     print('running WF CCA')
#     return_list <- tryCatch({
#       res <- WF_CCA(return_list[['df_miss']], col = "y", family = WF_params[['family']], R_PI = WF_params[['R_PI']])
#       return_list[['df_miss']] <- res$df
#       return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
#         res$district_df
#       return_list[['WF_betas']] <- res$betas
#       return_list
#     }, error = function(e){
#       return_list[['errors']][['WF']] <- rbind(return_list[['errors']][['WF']], data.frame(i = i, error = e[[1]]))
#       return_list
#     })
#     return_list[['timing']][['WF']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
#   }
#   
#   set.seed(i)
#   # run the CAR complete case analysis model
#   if('CARBayesST' %in% models){
#     t1 = Sys.time()
#     print('running CARBayes with CARBayesST')
#     return_list <- tryCatch({
#       res <- CARBayes_wrapper(return_list[['df_miss']], burnin = MCMC_params[['burnin.CARBayesST']], n.sample = MCMC_params[['n.sample.CARBayesST']], prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', MCMC_sampler = 'CARBayesST')
#       # rename the columns
#       colnames(res$df) <- gsub('y_pred_CAR', 'y_CARBayesST', colnames(res$df))
#       colnames(res$district_df) <- gsub('y_pred_CAR', 'y_CARBayesST', colnames(res$district_df))
#       
#       # update the results list
#       return_list[['df_miss']] <- res$df
#       return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
#       return_list[['CAR_summary']] <- res$CARBayesST_summary
#       return_list
#     }, error = function(e){
#       print(e)
#       return_list[['errors']][['CARBayesST']] <- rbind(return_list[['errors']][['CARBayesST']], data.frame(i, error = e[[1]]))
#       return_list
#     })
#     return_list[['timing']][['CARBayesST']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
#   }
#   
#   set.seed(i)
#   if('CARstan' %in% models){
#     t1 = Sys.time()
#     print('running CARBayes with stan')
#     return_list <- tryCatch({
#       res <- CARBayes_wrapper(return_list[['df_miss']], burnin = MCMC_params[['burnin.stan']], n.sample = MCMC_params[['n.sample.stan']], prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
#       
#       # rename the columns
#       colnames(res$df) <- gsub('y_pred_CAR', 'y_CARstan', colnames(res$df))
#       colnames(res$district_df) <- gsub('y_pred_CAR', 'y_CARstan', colnames(res$district_df))
#       # update the results list
#       return_list[['df_miss']] <- res$df
#       return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
#       return_list[['CARstan_summary']] <- res$CARstan_summary
#       return_list
#     }, error = function(e){
#       print(e)
#       return_list[['errors']][['CARstan']] <- rbind(return_list[['errors']][['CARstan']], data.frame(i, error = e[[1]]))
#       return_list
#     })
#     return_list[['timing']][['CARstan']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
#   }
# 
#   print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
#   
#   # check I got the correct result names
#   # outcome_name_checker(return_list, models = models)
#   
#   return(return_list)
# }
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
    }else if(model == 'CARstan'){
      fit_fxn <- function(df){
        tmp <- CARBayes_wrapper(df, burnin = params[['burnin']], n.sample = params[['n.sample']], prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
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
      }else if(model == 'CARstan'){
        return_list[['CARstan_summary']] <- res$CARstan_summary
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
                        params = list(R_PI = 200)),
                   list(model = 'WF_NB',
                        params = list(R_PI = 200)),
                   list(model = 'freqGLM',
                        params = list(R_PI = 200)),
                   list(model = 'CARstan',
                        params = list(burnin = params[['CARburnin']], n.sample = params[['CARnsample']])))
# Next (Maybe) CARBayesST

# run the models for each simulation dataset
system.time({
  imputed_list <- foreach(i=1:R_new) %dorng% one_run(lst, i, model_list)
})

true_betas <- lst$betas

save(imputed_list, seq, params, arguments, true_betas, file = results_file)