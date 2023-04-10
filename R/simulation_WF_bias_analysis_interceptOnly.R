library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')

# function to run the simulation for one set of beta values
run_sim <- function(b0, n = 50, R = 10000){
  
  # simulate y values and then estimate
  b0_estimates <- sapply(1:R, function(i){
    y = rpois(n, exp(b0))
    tmp = data.frame(y = rpois(n, exp(b0)))
    mod_col <- glm('y ~ 1', data = tmp, family=poisson)
    b0_hat <- mod_col$coefficients[1]
    return(b0_hat)
  })
  
  # get the bias
  mean_bias = mean(b0_estimates) - b0
  
  return(mean_bias)
}

# simulate for beta0 values ranging from 1 to 10
b0_vec = 1:10
bias_vec = sapply(b0_vec, function(b0){
  run_sim(b0, n = 1000, R = 1000)
})

# plot the results
plot(b0_vec, bias_vec, xlab = 'true b0', ylab = 'b0 bias')
abline(h = 0, col = 'red')
