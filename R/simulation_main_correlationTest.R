library(MASS)
# library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)

source('R/imputation_functions.R')
rstan_options(auto_write = TRUE)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id
inputs <- c('p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=WF:DGP_theta_shape=2.5:DGP_theta_rate=0.33:rho_DGP=0.2:alpha_DGP=0.2:tau2_DGP=1:R=20:num_jobs=2:output_path=NULL:family=poisson:theta=5:empirical_betas=F:CARburnin=1000:CARnsample=2000:R_PI=200:models=5,6,7\r','2')

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
  if(nn %in% c('p','rho_DGP','alpha_DGP','tau2_DGP','rho_MAR','alpha_MAR','tau2_MAR','gamma','theta', 'dispersion','CARburnin','CARnsample','R_PI','DGP_theta_shape','DGP_theta_rate')){
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
                 end_date = '2020-01-01', 
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
    if('DGP_theta_shape' %in% names(params)){
      arguments = c(arguments,
                    list(DGP_theta_shape = params[['DGP_theta_shape']]))
    }
    if('DGP_theta_rate' %in% names(params)){
      arguments = c(arguments,
                    list(DGP_theta_rate = params[['DGP_theta_rate']]))
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

arguments$district_sizes <- c(4,6,10)

# Simulate the data
lst <- do.call(simulate_data, arguments)

res1 <- lapply(lst$df_list, function(xx){
  get_morans_I(xx)
})

arguments2 <- arguments
arguments2$district_sizes <- c(4,4,6,6,10,10)

lst2 <- do.call(simulate_data, arguments)

res2 <- lapply(lst2$df_list, function(xx){
  get_morans_I(xx)
})

arguments3 <- arguments
arguments3$district_sizes <- c(rep(2,10), rep(3, 10))

lst3 <- do.call(simulate_data, arguments)

res3 <- lapply(lst3$df_list, function(xx){
  get_morans_I(xx)
})


#### Ok fine let's try CAR
inputs <- c('p=0.1:b0_mean=5.5:b1_mean=n0.25:missingness=mcar:DGP=CAR:DGP_theta_shape=2.5:DGP_theta_rate=0.33:rho_DGP=0.9:alpha_DGP=0.2:tau2_DGP=10:R=20:num_jobs=2:output_path=NULL:family=negbin:theta=5:empirical_betas=F:CARburnin=1000:CARnsample=2000:R_PI=200:models=5,6,7\r','2')

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
    if(nn %in% c('p','rho_DGP','alpha_DGP','tau2_DGP','rho_MAR','alpha_MAR','tau2_MAR','gamma','theta', 'dispersion','CARburnin','CARnsample','R_PI','DGP_theta_shape','DGP_theta_rate')){
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
                   end_date = '2020-01-01', 
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
      if('DGP_theta_shape' %in% names(params)){
        arguments = c(arguments,
                      list(DGP_theta_shape = params[['DGP_theta_shape']]))
      }
      if('DGP_theta_rate' %in% names(params)){
        arguments = c(arguments,
                      list(DGP_theta_rate = params[['DGP_theta_rate']]))
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

lst <- do.call(simulate_data, arguments)

res1 <- lapply(lst$df_list, function(xx){
  get_morans_I(xx)
})
