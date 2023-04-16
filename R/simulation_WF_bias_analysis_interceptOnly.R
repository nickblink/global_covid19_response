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
  run_sim(b0, n = 50, R = 10000)
})

# plot the results
plot(b0_vec, bias_vec, xlab = 'true b0', ylab = 'b0 bias')
abline(h = 0, col = 'red')


#### Now with b0 and b1 ####
# function to run the simulation for one set of beta values
run_sim <- function(b0, b1, n = 50, R = 10000){
  
  # simulate y values and then estimate
  beta_estimates <- t(sapply(1:R, function(i){
    x1 = rnorm(n)
    tmp = data.frame(y = rpois(n, exp(b0 + b1*x1)), x = x1)
    mod_col <- glm('y ~ x', data = tmp, family=poisson)
    beta_hat <- mod_col$coefficients
    return(beta_hat)
  }))
  
  # get the bias
  mean_bias = colMeans(beta_estimates) - c(b0, b1)
  
  return(mean_bias)
}

# simulate for beta0 values ranging from 1 to 10
b0_vec = 1:10
b1 = 1
bias = t(sapply(b0_vec, function(b0){
  run_sim(b0, b1, n = 50, R = 10000)
}))

# plot the results
plot(b0_vec, bias[,1], xlab = 'true b0', ylab = 'b0 bias')
abline(h = 0, col = 'red')

plot(b0_vec, bias[,2], xlab = 'true b0', ylab = 'b1 bias')
abline(h = 0, col = 'red')


# simulate for beta0 values ranging from 1 to 10
b0 = 1
b1_vec = -5:5
bias2 = t(sapply(b1_vec, function(b1){
  run_sim(b0, b1, n = 50, R = 10000)
}))

plot(b1_vec, bias2[,1], xlab = 'true b1', ylab = 'b0 bias')
abline(h = 0, col = 'red')

plot(b1_vec, bias2[,2], xlab = 'true b1', ylab = 'b1 bias')
abline(h = 0, col = 'red')
