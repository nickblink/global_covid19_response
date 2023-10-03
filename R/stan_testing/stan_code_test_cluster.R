library(dplyr)
library(rstan)
rstan_options(auto_write = T)

#
#### Take 1 - intercept only ####
pois_sdata <- list(
  N = 100,  # number of observations
  y = rpois(100, lambda = 30)  # outcome variable (yellow card)
)

m3 <- stan(model_code = 'data {
  int<lower=0> N;  // number of observations
  int y[N];  // data array (counts);
}
parameters {
  real log_lambda;  // log of rate parameter
}
model {
  y ~ poisson_log(log_lambda);
  // prior
  log_lambda ~ normal(0, 5);
}
generated quantities {
  real lambda = exp(log_lambda);
  int yrep[N];
  for (i in 1:N) {
    yrep[i] = poisson_log_rng(log_lambda);
  }
}
', data = pois_sdata, 
           # Below are optional arguments
           warmup = 100,
           iter = 200,
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m3, pars = c("lambda", "log_lambda"))


shinystan::launch_shinystan(m3)
# that's cool. I'll explore that more later on. But that's very neat. Thank you to whoever made that.

#### Take 2: Let's do a simple regression ####

beta = 0.3

x = rnorm(100, 10, 5)
pois_sdata <- list(
  N = 100,  # number of observations
  x = x, 
  y = rpois(100, lambda = exp(x*beta))  # outcome variable (yellow card)
)

m3 <- stan(file = "poisson_model2.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m3, pars = c("beta"))

#### Take 3 - multivariate regression! ####

beta = c(0.3, 0.6)

X = data.frame(ind_1 = rnorm(100, 10, 5), ind_2 = rnorm(100, 5,3))

pois_sdata <- list(
  N = 100,  # number of observations
  p = ncol(X),
  X = X, 
  y = rpois(100, lambda = exp(X[,1]*0.3 + X[,2]*0.6))  # outcome variable 
)

m3 <- stan(file = "poisson_model3.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m3, pars = c("beta"))


shinystan::launch_shinystan(m3)
#
#### Take 4 - just passing a string ####

m3 <- stan(model_code = "data {
  int<lower=0> N;  // number of observations
  int y[N];  // data array (counts);
}
parameters {
  real log_lambda;  // log of rate parameter
}
model {
  y ~ poisson_log(log_lambda);
  // prior
  log_lambda ~ normal(0, 5);
}
generated quantities {
  real lambda = exp(log_lambda);
  int yrep[N];
  for (i in 1:N) {
    yrep[i] = poisson_log_rng(log_lambda);
  }
}
", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           chains = 4, 
           cores = min(parallel::detectCores(), 4))
