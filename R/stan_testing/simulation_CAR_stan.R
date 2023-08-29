library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)

source('../imputation_functions.R')

# data creation
{

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

}

#### Just one facility linear regression ####
df2 = df %>%
  filter(facility == 'A1')

formula = as.formula("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")

X = model.matrix(formula, data = df2)
y = df2$y

pois_sdata <- list(
  N = nrow(X), # number of observations
  p = ncol(X), # number of variables
  X = X, # design matrix
  y = y)  # outcome variable 

m3 <- stan(file = "single_fac_reg.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m3, pars = c("beta"))

params = extract(m3)


#
#### QR style ####
m4 <- stan(file = "QR_regression.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m4, pars = c("beta"))

print(m4, pars = c('y_exp', 'y_exp2'))
# it seems like the predictions the same. Noice. Thank you QR decomposition and the people who figured out that it improves modeling.

#
#### Multiple facilities regression ####
formula = as.formula("y ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3")

df3 = df %>% filter(district == 'B')
# works with A
# works with B on single_fac_reg but not QR_regression

X = model.matrix(formula, data = df3)
y = df3$y

pois_sdata <- list(
  N = nrow(X), # number of observations
  p = ncol(X), # number of variables
  X = X, # design matrix
  y = y)  # outcome variable 

m5 <- stan(file = "QR_regression.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 8000, 
           warmup = 4000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# Stan model 'anon_model' does not contain samples.
# What does that mean?
# this is because the mean of the model is too damn big

print(m5, pars = c("beta"))
tt = extract(m5)

max(tt$y_exp)
max(tt$y_pred)

y_pred = colMeans(tt$y_pred)
y_exp = colMeans(tt$y_exp)
plot(y, y_pred)
plot(y, y_exp)

#### Setting MVN priors ####
formula = as.formula("y ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3")

X = model.matrix(formula, data = df)
y = df$y

lm_fit <- glm(formula, family = 'poisson', data = df)
coef_mat <- summary(lm_fit)$coefficients
prior_mean_beta <- coef_mat[,1]
sigma_beta = 10*vcov(lm_fit)

pois_sdata <- list(
  N = nrow(X), # number of observations
  p = ncol(X), # number of variables
  X = X, # design matrix
  y = y,
  mu = prior_mean_beta,
  Sigma = sigma_beta)  # outcome variable 

m5 <- stan(file = "regression_priors.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 8000, 
           warmup = 4000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# most definitely slower. So be it.
check_treedepth(m5)
# apparently this is a problem if we have saturated maximum tree depth, because that indicates premature ending of the NUTS sampler.
