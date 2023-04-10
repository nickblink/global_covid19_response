library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')

run_sim <- function(b0, b1 = -0.25, seasonal = T, R = 10000){
  df = initialize_df(district_sizes = 1)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  if(seasonal){
    betas = matrix(data = c(b0, b1, rnorm(6, 0, 0.25)), ncol = 1, dimnames = list(c("intercept","year","cos1","sin1","cos2","sin2","cos3","sin3"), NULL))
  }else{
    betas = matrix(data = c(b0, b1, rep(0, 6)), ncol = 1, dimnames = list(c("intercept","year","cos1","sin1","cos2","sin2","cos3","sin3"), NULL))
  }
  
  print(betas)
  
  # keep the 1 for intercepts
  X = df %>% 
    mutate(intercept = 1) %>%
    dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3) %>%
    as.matrix()
  
  mu = X%*%betas
  
  sim_lst <- lapply(1:R, function(i){
    tmp = df
  
    # simluate random values
    tmp$y = rpois(length(mu), exp(mu))
    
    return(tmp)
  })

  model_form <- as.formula("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  
  beta_hats <- sapply(sim_lst, function(sim_data){
    mod_col <- glm(model_form, data = sim_data, family=poisson)
    beta_hat <- mod_col$coefficients
    return(beta_hat)
  })
  
  bias_res <- as.data.frame(t(sapply(1:R, function(ii){
    y <- sim_lst[[ii]]$y
    tmp <- matrix(c(sum(y == 0), sum(y), beta_hats[1, ii] - b0))
    return(tmp)
  })))
  colnames(bias_res) = c('zeros', 'sum_y', 'bias')
  
  mean_bias = rowMeans(beta_hats)[1] - b0
  median_bias = median(beta_hats[1,]) - b0
  
  lst = list(mean_bias = mean_bias, median_bias = median_bias, bias_res)
  return(lst)
}

b0_vec = 1:10
bias_lst = lapply(b0_vec, function(b0){
  run_sim(b0, b1 = 0, seasonal = F, R = 10000)
})

mean_bias = unlist(lapply(bias_lst, '[[', 1))
median_bias = unlist(lapply(bias_lst, '[[', 2))

plot(b0_vec, mean_bias, xlab = 'true b0 (intercept)', ylab = 'b0 estimate bias', main = 'black = mean; blue = median')
points(b0_vec, median_bias, col = 'blue')
abline(a = 0, b = 0, col= 'red')
