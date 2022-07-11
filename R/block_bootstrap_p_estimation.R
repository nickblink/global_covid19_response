setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)
library(boot)

optimal_block_length <- function(D, m, block_vec = 1:10, metric = 'MSE_quantiles'){
  uni_dates <- unique(D$date) %>% sort()
  n = length(uni_dates)
  D$date_fac <- paste0(D$date, '_', D$facility)
  
  # quant_cols <- grep('pred_freqGLMepi_', colnames(tt), value = T)
  quant_cols <- c("y_pred_freqGLMepi_0.05", "y_pred_freqGLMepi_0.25", "y_pred_freqGLMepi_0.5", "y_pred_freqGLMepi_0.75", "y_pred_freqGLMepi_0.95")
  
  # warning('running optimal block length for only one facility')
  # tmp = D
  
  full_boot_list = list()
  boot_list <- list()
  res_list <- list()
  
  system.time({
  # cycle through probabilities
  for(blocksize in block_vec){
    print(blocksize)
    # run bootstrap on full D for this block size
    tt = freqGLMepi_imputation(D, prediction_intervals = 'stationary_bootstrap', R_PI = 100, scale_by_num_neighbors = T, blocksize = blocksize)$df
    
    # store results for this block size
    tmp <- list(tt[,c('date', 'facility', 'date_fac', quant_cols)])
    names(tmp) <- blocksize
    full_boot_list = append(full_boot_list, tmp)
    
    # cycle through subsets of data - no wraparound
    #for(i in 1:(n - m + 1)){
    
    # cycle through with wraparound
    for(i in 1:n){
      if(i + m - 1 <= n){
        ind <- i:(i + m - 1)
      }else{
        ind <- c(i:n, 1:(m - n + i -1))
        
      }
 
      D2 <- D %>% filter(date %in% uni_dates[ind])

      # run bootstrap on D2
      tt = freqGLMepi_imputation(D2, prediction_intervals = 'stationary_bootstrap', R_PI = 100, scale_by_num_neighbors = T, blocksize = blocksize)$df
      
      # store the results
      boot_list <- append(boot_list, list(tt[,c('date', 'facility', 'date_fac', quant_cols)]))
    }
       
    if(metric == 'MSE_quantiles'){
      # include quantiles at 5%, 25%, 50%, 75%, 95% for MSE estimation
      full <- full_boot_list[[as.character(blocksize)]]
      
      res_df <- sapply(boot_list, function(xx){
        # match data frames by date and facility
        xx$date_fac <- paste0(xx$date, '_', xx$facility)
        tmp <- full[match(xx$date_fac, full$date_fac),]
        
        # get the MSEs
        diff <- xx[,quant_cols] - tmp[,quant_cols]
        MSEs <- apply(diff, 2, function(yy) {mean(yy^2, na.rm = T)})
        return(MSEs)
      }) %>%
        t()
      
      tmp <- list(list(res_all = res_df, res_mean = colMeans(res_df)))
      names(tmp) <- blocksize
      res_list <- append(res_list, tmp)
      
    }else{
      stop('please input a proper metric')
    }
    
  }
  })
  
  #save(res_list, file = 'results/tmp_blocklength_test_m40.RData')
  #save(res_list, file = 'results/tmp_blocklength_test_m24.RData')
  
  # create table of results
  
  # pick the final p!
  
  # scale up the probability/block length by the diff of m and n from the equation in the paper
  
}


res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 2, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

tmp = res$df_list[[1]]
tmp$y_true = tmp$y
D = tmp

freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'stationary_bootstrap', R_PI = 100, scale_by_num_neighbors = T) 





