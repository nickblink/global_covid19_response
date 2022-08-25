### Now doing this as an R script rather than Rmd because it's easier to work with.
#setwd('C:/Users/nickl/Documents/global_covid19_response/')
setwd('C:/Users/nickl/Documents/github_projects/global_covid19_response')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)
#library(doRNG)
#library(doParallel)

params <- commangArgs(trailingOnly = TRUE)

R <- params[[1]]
num_jobs <- params[[2]]
job_id <- params[[3]]

if(job_id < num_jobs){
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):(floor(R/num_jobs)*job_id)
}else{
  seq <- (floor(R/num_jobs)*(job_id - 1) + 1):R
}
 
job_name <- 'test'

# set up the output folder
tmp_data_path <- 'data/tmp_data/'
output_path <- paste0(tmp_data_path, job_name)
if(!file.exists(output_path)){
  dir.create(output_path)
}

#### MCAR p = 0.2 freqGLM_epi - 20 years ####
lst <- res <- simulate_data_freqGLM_epi(district_sizes = 4, R = R, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

if(job_id == 1){
  save(lst, file = 'simulated_data.RData')
}

# registerDoParallel(cores = 4)

one_run <- function(lst, i){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
  
  # run the freqGLM_epi imputation
  freqGLMepi_list = freqGLMepi_imputation(df_miss, prediction_intervals = 'bootstrap', R_PI = 100, verbose = F)
  df_miss = freqGLMepi_list$df
  
  # run the periodic imputation
  periodic_list = periodic_imputation(df_miss, col = "y", family = 'poisson', group = 'facility', R_PI = 100)
  df_miss = periodic_list$df
  
  #  run the CARBayes imputation with different intercepts by facility
  CAR_list2 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_intercept')
  df_miss = CAR_list2$facility_df
  colnames(df_miss) = gsub('CARBayes_ST', 'CB_intercept', colnames(df_miss))

  #  run the CARBayes imputation with different coeffs by facility
  CAR_list3 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')
  df_miss = CAR_list3$facility_df
  colnames(df_miss) = gsub('CARBayes_ST', 'CB_facility', colnames(df_miss))

  return(df_miss)
}

# This might not work. %do% worked but not %dorng%. We'll see on the cluster
system.time({
imputed_list <- foreach(i=seq) %dorng% one_run(lst, i)
})

save(imputed_list, file = sprintf('%s/sim_results_%i.RData', output_path, job_id))

