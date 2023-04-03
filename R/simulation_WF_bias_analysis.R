library(MASS)
library(Matrix)
library(dplyr)

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
  
  bias = rowMeans(beta_hats)[1] - b0
  
  return(bias)
}

b0_vec = 1:10
bias_vec = sapply(b0_vec, function(b0){
  run_sim(b0, b1 = 0, seasonal = F, R = 1000)
})

plot(b0_vec, bias_vec, xlab = 'true b0 (intercept)', ylab = 'b0 estimate bias')
abline(a = 0, b = 0, col= 'red')
