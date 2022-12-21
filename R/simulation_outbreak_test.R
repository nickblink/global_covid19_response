library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
#library(cowplot)

setwd('C:/Users/Admin-Dell/Documents/github_projects/global_covid19_response')
source('R/imputation_functions.R')

R = 1000
simulate_data_testing <- function(district_sizes = 1, R = R, empirical_betas = F, seed = 10, ...){
  # set seed
  set.seed(seed)
  # set up data frame
  df = initialize_df(district_sizes, ...)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility names
  facilities = unique(df$facility)
  
  betas = matrix(c(3.108460605,	-0.334426667,	-0.253100157,	0.473832506,	-0.21421039, 0.351340989,	-0.123478103,	0.342161078), nrow = 1)
  
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
  rownames(betas) = facilities
  
  # initialize list of data frames
  df_lst = list()
  
  # make n sampled sets of data
  for(i in 1:R){
    
    # simulate values given the betas
    tmp_lst = lapply(facilities, function(xx){
      tmp = df %>% filter(facility == xx)
      
      # keep the 1 for intercepts
      X = tmp %>% 
        mutate(intercept = 1) %>%
        dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
      
      # error checking
      if(!identical(colnames(betas), colnames(X))){
        print(colnamaes(betas))
        print(colnames(X))
        stop(sprintf('colnames of betas not equal to X: %s, %s',paste(colnames(betas), collapse = ';'), paste(colnames(X), collapse = ';') ))
      }
      
      # make the 8x1 beta vector for this facility
      beta_f = t(betas[xx,,drop = F])
      
      # get mean prediction from linear model
      mu = as.matrix(X)%*%beta_f
      tmp$y_exp = exp(mu)
      
      # simluate random values
      tmp$y = rpois(length(mu), exp(mu))
      
      return(tmp)
    })
    
    # combine values into one data frame
    df_lst[[i]] = do.call('rbind', tmp_lst)
    
  }
  
  # make list of values to return
  res_lst = list(df_list = df_lst, betas = betas)
  return(res_lst)
}

system.time({
lst <- simulate_data_testing(district_sizes = 1, R = 1000, end_date = '2020-12-01')
})
# 10/11s

# function to run all models for a specific dataset
one_run <- function(lst, i){
  t0 <- Sys.time()
  df = lst$df_list[[i]]
  
  # run the WF complete case analysis model
  print('running CCA')
  res <- WF_CCA(df, col = "y", family = 'poisson', R_PI = 200)
    
  print(sprintf('i = %i; time = %f minutes', i, difftime(Sys.time(), t0, units = 'm')))
  
  return(res)
}

set.seed(1)
system.time({
  res_lst <- lapply(1:R, function(ii) {one_run(lst, ii)})
})
# 1000s

save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/results_fixed_betas_12212022.RData')


# 
# simulate_data_testing2 <- function(district_sizes = 1, R = 1, empirical_betas = F, seed = 10, ...){
#   # set seed
#   set.seed(seed)
#   # set up data frame
#   df = initialize_df(district_sizes, ...)
#   
#   # make periodic covariates
#   df = add_periodic_cov(df)
#   
#   # get all facility names
#   facilities = unique(df$facility)
#   
#   betas = matrix(c(3.108460605,	-0.334426667,	-0.253100157,	0.473832506,	-0.21421039, 0.351340989,	-0.123478103,	0.342161078), nrow = 1)
#   
#   colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
#   rownames(betas) = facilities
#   
#   # initialize list of data frames
#   df_lst = list()
#   
#   # make n sampled sets of data
#   # for(i in 1:R){
#   df_lst <- lapply(1:R, function(i){ 
#     
#     # simulate values given the betas
#     tmp_lst = lapply(facilities, function(xx){
#       tmp = df %>% filter(facility == xx)
#       
#       # keep the 1 for intercepts
#       X = tmp %>% 
#         mutate(intercept = 1) %>%
#         dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
#       
#       # error checking
#       if(!identical(colnames(betas), colnames(X))){
#         print(colnamaes(betas))
#         print(colnames(X))
#         stop(sprintf('colnames of betas not equal to X: %s, %s',paste(colnames(betas), collapse = ';'), paste(colnames(X), collapse = ';') ))
#       }
#       
#       # make the 8x1 beta vector for this facility
#       beta_f = t(betas[xx,,drop = F])
#       
#       # get mean prediction from linear model
#       mu = as.matrix(X)%*%beta_f
#       tmp$y_exp = exp(mu)
#       
#       # simluate random values
#       tmp$y = rpois(length(mu), exp(mu))
#       
#       return(tmp)
#     })
#     
#     # combine values into one data frame
#     tt <- do.call('rbind', tmp_lst)
#     
#     tt
#   })
#   # make list of values to return
#   res_lst = list(df_list = df_lst, betas = betas)
#   return(res_lst)
# }
# 
# 
# system.time({
#   lst <- simulate_data_testing2(district_sizes = 1, R = 1000, end_date = '2020-12-01')
# })
# # samesies
