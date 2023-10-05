library(MASS)
library(CARBayesST)
library(Matrix)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)
library(dplyr)

source('R/imputation_functions.R')
rstan_options(auto_write = TRUE)

# register the cores
registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id
inputs <- c('0.1', '6', 'n0.25', 'mcar', 'CAR', '0.3', '0.3', '1', '10','5','test','1')

### Get parameters (for home and cluster)
{
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

load("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_car331_beta6_n025_20230603/simulated_data.RData" )

print('data made')

# initialize the error catcher for each run
errors <- list(freqEpi = data.frame(i = NULL, error = NULL),
               WF = data.frame(i = NULL, error = NULL),
               CARBayesST = data.frame(i = NULL, error = NULL),
               CARstan = data.frame(i = NULL, error = NULL))

}

# function to run all models for a specific dataset
one_run <- function(lst, i, models = c('freq', 'WF', 'CARBayesST','CARstan'), WF_params = list(R_PI = 200), freqGLM_params = list(R_PI = 200), MCMC_params = list(burnin.stan = 1000, n.sample.stan = 2000, burnin.CARBayesST = 5000, n.sample.CARBayesST = 10000)){
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
  return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
  rm(df_miss)
  
  # run the freqGLM_epi complete case analysis
  if('freq' %in% models){
    t1 = Sys.time()
    print('running freqGLM_epi')
    return_list <- tryCatch({
      freqGLMepi_list = freqGLMepi_CCA(return_list[['df_miss']], R_PI = freqGLM_params[['R_PI']], verbose = F)
      return_list[['df_miss']] <- freqGLMepi_list$df
      return_list
    }, error = function(e){
      return_list[['errors']][['freqEpi']] <- rbind(return_list[['errors']][['freqEpi']], data.frame(i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['freq']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  
  # run the WF complete case analysis model
  if('WF' %in% models){
    t1 = Sys.time()
    print('running WF CCA')
    return_list <- tryCatch({
      res <- WF_CCA(return_list[['df_miss']], col = "y", family = 'poisson', R_PI = WF_params[['R_PI']])
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
        res$district_df
      return_list[['WF_betas']] <- res$betas
      return_list
    }, error = function(e){
      return_list[['errors']][['WF']] <- rbind(return_list[['errors']][['WF']], data.frame(i = i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['WF']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  # run the CAR complete case analysis model
  if('CARBayesST' %in% models){
    t1 = Sys.time()
    print('running CARBayes with CARBayesST')
    return_list <- tryCatch({
      res <- CARBayes_wrapper(return_list[['df_miss']], burnin = MCMC_params[['burnin.CARBayesST']], n.sample = MCMC_params[['n.sample.CARBayesST']], prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', MCMC_sampler = 'CARBayesST')
      #res <- CARBayes_wrapper(return_list[['df_miss']], burnin = 10000, n.sample = 20000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
      colnames(res$df) <- gsub('y_pred_CAR', 'y_CARBayesST', colnames(res$df))
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      return_list
    }, error = function(e){
      print(e)
      return_list[['errors']][['CARBayesST']] <- rbind(return_list[['errors']][['CARBayesST']], data.frame(i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['CARBayesST']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }
  
  if('CARstan' %in% models){
    t1 = Sys.time()
    print('running CARBayes with stan')
    return_list <- tryCatch({
      res <- CARBayes_wrapper(return_list[['df_miss']], burnin = MCMC_params[['burnin.stan']], n.sample = MCMC_params[['n.sample.stan']], prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
      #res <- CARBayes_wrapper(return_list[['df_miss']], burnin = 10000, n.sample = 20000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
      colnames(res$df) <- gsub('y_pred_CAR', 'y_CARstan', colnames(res$df))
      return_list[['df_miss']] <- res$df
      return_list[['district_df']] <- merge(return_list[['district_df']], res$district_df, by = c('district','date'))
      return_list[['CARstan_ESS']] <- res$ESS
      return_list
    }, error = function(e){
      print(e)
      return_list[['errors']][['CARstan']] <- rbind(return_list[['errors']][['CARstan']], data.frame(i, error = e[[1]]))
      return_list
    })
    return_list[['timing']][['CARstan']] <- as.numeric(difftime(Sys.time(), t1, units = 'm'))
  }

  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  # check I got the correct result names
  # outcome_name_checker(return_list, models = models)
  
  return(return_list)
}

set.seed(1)
# lst_OG <- lst
# lst_new <- lst
i = 208
df = lst$df_list[[i]]
df$facility = factor(df$facility, levels = unique(df$facility))
if(missingness == 'mcar'){
  df_miss = MCAR_sim(df, p = p, by_facility = T)
}else if(missingness == 'mar'){
  df_miss = MAR_spatiotemporal_sim(df, p = p, rho = params[['rho_MAR']], alpha = params[['alpha_MAR']], tau = params[['tau2_MAR']])
}else{
  df_miss <- MNAR_sim(df, p = p, direction = 'upper', gamma = gamma, by_facility = T)
}

# initializing the return list
return_list <- list(df_miss = df_miss, district_df = lst$district_list[[i]], errors = errors, timing = list())
rm(df_miss)

# run the freqGLM_epi complete case analysis
freqGLMepi_list = freqGLMepi_CCA(return_list[['df_miss']], R_PI = 10, verbose = F)
tt = freqGLMepi_list$df

# res <- one_run(lst, i, models = c('freq'), freqGLM_params = list(R_PI = 20), MCMC_params = list(burnin.stan = 1000, n.sample.stan = 2000, burnin.CARBayesST = 5000, n.sample.CARBayesST = 10000))

true_betas <- lst$betas

load("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_car331_beta6_n025_20230603/sim_results_p0.0_mcar_11(50).RData")

df = imputed_list[[8]]$df_miss
freqGLMepi_list = freqGLMepi_CCA(df, R_PI = 10, verbose = F)
