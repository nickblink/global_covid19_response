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
#registerDoParallel(cores = 20)

# get the parameters (first line is for testing on my home computer)
# p b0 b1 missingness ST rho alpha tau2 R #jobs name_output job_id
inputs <- c('0.1', '6', 'n0.25', 'mcar', 'CAR', '0.3', '0.3', '1', '10','5','test','1')
#inputs <- commandArgs(trailingOnly = TRUE)
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

# Simulate the data
if(DGP == 'nost'){
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

# initialize the error catcher for each run
errors <- list(freqEpi = data.frame(i = NULL, error = NULL),
                 WF = data.frame(i = NULL, error = NULL),
                 CARBayes = data.frame(i = NULL, error = NULL))


df = lst$df_list[[2]]
df_miss = MCAR_sim(df, p = p, by_facility = T)
# initializing the return list
return_list <- list(df_miss = df_miss, district_df = lst$district_list[[1]], errors = errors, WF_betas = NULL, CAR_summary = NULL)

#### Testing different CAR models ####
res1 = CARBayes_wrapper(return_list[['df_miss']], burnin = 4000, n.sample = 16000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 10, return_chain = T)

names(res1)

# test the ESS calculations
head(res1$summary_stats)

betas = res1$model_chain$samples$beta
ESS1 = coda::effectiveSize(betas)
diff1 = ESS1 - res1$summary_stats[1:160,6]
plot(density(abs(diff1/ESS1)))
# ok so the results are very close. Max 3% off
# I'm gonna stop here. I'm not gonna worry about mcmcse::ess. Because this is damn close enough

res2 = CARBayes_wrapper(return_list[['df_miss']], burnin = 4000, n.sample = 16000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 10, return_chain = T, S_glm_debug = T, return_raw_fit = T)

df = data.frame(CARBayesST = res1$summary_stats[1:160,6],
                S_glm = res2$model_chain$summary.results[,6])
colMeans(df)

# burn = 1k, n = 2k: S_glm is on average 5X better (avg ESS = 19 compared to 4). However, it has a lot of zeros. Meaning the space was not explored at all for those params. Should I run a bigger sample here to be sure? My comp will probably crash.

# burn = 2k, n =4k: S_glm is on average 10x better (avg ESS = 50 compared to 5). Still some zeroes this time, but not as much as last time. The zeroes were all on a few of the cos3's and sin3's. Could this be an inherent problem with the model? Some ESS are at 200 and some even a bit over, which doesn't make sense to me.

# burn = 4k, n = 8k: S_glm is almost 20x better (avg ESS = 110 compared to 5.7) (median 50:5). Still the cos3's and sin3's are getting zeroes. Why the f is that happening?

# burn = 4k, n = 16k: S_glm is almost 20x better (avg ESS = 146 compared to 9.5) (median 66:8). Now getting some 0s on the cos1's. It appears, maybe, that the betas are sampled in MCMC in groups of 5. Because there's always 5 zeroes in a row

#### Testing priors ####
res1 = CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 10, return_chain = T)

res2 = CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 10, return_chain = T, prior = 'WF', prior_var_scale = 1)

res3 = CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 10, return_chain = T, prior = 'constant', prior_mean = rep(2, 160), prior_var = rep(1, 160))

#### Test thinning and n.sample ####
# run the models
set.seed(10)
system.time({
res <- CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin  = 1, return_chain = T)
}) # 182s

set.seed(10)
system.time({
res2 <- CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 2000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin = 10, return_chain = T)
}) # 44s

set.seed(10)
system.time({
  res3 <- CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 5000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin = 50, return_chain = T)
}) # 56s

set.seed(10)
system.time({
  res4 <- CARBayes_wrapper(return_list[['df_miss']], burnin = 1000, n.sample = 5000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01', thin = 10, return_chain = T)
}) # 72s

names(res)
tt = res$model_chain$samples$beta
ss = res2$model_chain$samples$beta

# ok they are the same. Good.
sum(abs(ss[1,]) - abs(tt[10,]))
sum(abs(ss[10,]) - abs(tt[100,]))

diff = res$summary_stats - res2$summary_stats

View(res$summary_stats)
View(res2$summary_stats)
View(diff)



#### Analyzing the n.sample results  ####
load('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/CARBayes_debug_res2_07092023.RData')

# compare the convergence statistics
# Compare the beta convergence from each of these (will require averaging the convergence or calling the getting results function)
tt <- res$summary_stats
ind = which(!(rownames(tt) %in% c('tau2','rho.S','rho.T')))

apply(res$summary_stats[ind,], 2, mean)

apply(res2$summary_stats[ind,], 2, mean)

apply(res3$summary_stats[ind,], 2, mean)

apply(res4$summary_stats[ind,], 2, mean)

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

true_betas <- lst$betas

save(imputed_list, seq, params, true_betas, file = results_file)

