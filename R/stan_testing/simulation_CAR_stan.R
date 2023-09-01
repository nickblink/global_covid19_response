library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)
library(bayestestR)

source('../imputation_functions.R')

#### data creation ####
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

#
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

m1 <- stan(file = "single_fac_reg.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m1, pars = c("beta"))

params = extract(m1)


#
#### QR style ####
m2 <- stan(file = "QR_regression.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 4, 
           cores = min(parallel::detectCores(), 4))

print(m2, pars = c("beta"))

print(m2, pars = c('y_exp'))
# it seems like the predictions the same (when I am also computing y_exp2 as a test, that is). Noice. Thank you QR decomposition and the people who figured out that it improves modeling.

#
#### Multiple facilities regression ####
formula = as.formula("y ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3")

df3 = df %>% filter(district == 'C')
# works with A
# works with B on single_fac_reg but not QR_regression (maybe not true - it seems to work when I run it later).
# works with C

X = model.matrix(formula, data = df3)
y = df3$y

pois_sdata <- list(
  N = nrow(X), # number of observations
  p = ncol(X), # number of variables
  X = X, # design matrix
  y = y)  # outcome variable 

m3 <- stan(file = "QR_regression.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 8000, 
           warmup = 4000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# Stan model 'anon_model' does not contain samples.
# What does that mean?
# this is because the mean of the model is too damn big

print(m3, pars = c("beta"))
tt = extract(m3)

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

m4 <- stan(file = "regression_priors.stan", data = pois_sdata, 
           # Below are optional arguments
           iter = 8000, 
           warmup = 4000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# most definitely slower. So be it.
check_treedepth(m4)
# apparently this is a problem if we have saturated maximum tree depth, because that indicates premature ending of the NUTS sampler.

#### Running a simple CAR model ####
# I will run a model with a CAR random effect. This will not have a parameter for spatial correlation nor will it have any temporal correlation. 
# I will use the pairwise formulation of the likelihood this model so that I will input edges and nodes for the adjacency matrix and compute from that rather than computing the precision matrix directly. I can later test if this is faster or not.
# So I will have the same mean model with the CAR random effect for each location (not differing over time, remember). Hmmm I guess it will be constant over time? That seems wack. I will add in independent CAR random effects for each time point. How about that? Not perfectly ideal, but I'll have to figure this out at some point anyway


### Preparing the stan data
# This function takes in the data frame and prepares the data to be input to stan for a CAR model. 
# The df is ordered correctly and then the nodes are created so that phi can be a vector with the time points stacked on top of one another
prep_stan_data <- function(df, formula){
  N_T <- length(unique(df$date))
  N_F <- length(unique(df$facility))
  
  # order the data frame
  df <- df %>% arrange(date, facility)
  
  W = make_district_adjacency(df)
  
  base_df = data.frame(node1 = NULL, node2 = NULL)
  
  # Make the base df for one time point
  for(i in 1:(nrow(W)-1)){
    j = which(W[i,] == 1)
    j = j[j > i]
    base_df = rbind(base_df, 
                    data.frame(node1 = rep(i, length(j)),
                               node2 = j))
  }
  
  # Make the repeated nodes for all time points
  node_df = base_df
  for(t in 2:N_T){
    node_df = rbind(node_df, base_df + (t-1)*N_F)
  }
  
  # make the model matrix and outcome
  X = model.matrix(formula, data = df)
  y = df$y
  
  # compute the priors
  lm_fit <- glm(formula, family = 'poisson', data = df)
  coef_mat <- summary(lm_fit)$coefficients
  prior_mean_beta <- coef_mat[,1]
  sigma_beta = 10*vcov(lm_fit)
  
  # make the stan data frame
  stan_data <- list(
    N = nrow(X), # number of observations
    p = ncol(X), # number of variables
    N_F = N_F, # number of facilities
    N_T = N_T, # number of time points
    X = X, # design matrix
    y = y, # outcome variable 
    mu = prior_mean_beta, # prior mean
    Sigma = sigma_beta, # prior variance
    N_edges = nrow(node_df),
    node1 = node_df$node1, # pairs of nodes that are adjacent
    node2 = node_df$node2)  
  
  return(stan_data)
}

formula = as.formula("y ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3")

stan_data <- prep_stan_data(df, formula)

m5 <- stan(file = "regression_ICAR.stan", data = stan_data, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# well it didn't crash!


#### Running a Leroux CAR model - Q matrix style ####

prep_stan_data_leroux <- function(df, formula){
  N_T <- length(unique(df$date))
  N_F <- length(unique(df$facility))
  
  # order the data frame
  df <- df %>% arrange(date, facility)
  
  W_star = make_district_W2_matrix(df)
  
  # make the model matrix and outcome
  X = model.matrix(formula, data = df)
  y = df$y
  
  # compute the priors
  lm_fit <- glm(formula, family = 'poisson', data = df)
  coef_mat <- summary(lm_fit)$coefficients
  prior_mean_beta <- coef_mat[,1]
  sigma_beta = 10*vcov(lm_fit)
  
  # make the stan data frame
  stan_data <- list(
    N = nrow(X), # number of observations
    p = ncol(X), # number of variables
    N_F = N_F, # number of facilities
    N_T = N_T, # number of time points
    X = X, # design matrix
    y = y, # outcome variable 
    mu = prior_mean_beta, # prior mean
    Sigma = sigma_beta, # prior variance
    W_star = W_star, 
    I = diag(1.0, N_F))
  
  return(stan_data)
}

formula = as.formula("y ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3")

stan_data <- prep_stan_data_leroux(df, formula)

m6 <- stan(file = "regression_leroux.stan", data = stan_data, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

# it ran! But the ESS is too low? Not sure. Moving on anyway

shinystan::launch_shinystan(m6)
tt = extract(m6, permuted = F)
print(m6, pars = 'beta')

#### Running a Rushworth CAR model - Leggo baby! ####
# all I need to add is the autocorrelation component. I got dis.

# (run the previous data setup for leroux. It's the same here)
m7 <- stan(file = "regression_rushworth.stan", data = stan_data, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))


# ok so it ran in 6-7 minutes for 2000 iterations and 1 chain. Getting an NA for at least one of the rhats. No bueno. And there are apparently very low bulk ESS and tail ESS that are too low.

print(m7, pars = 'beta')
stan_est = extract(m7, pars = 'beta', permuted = F)
ESS <- effective_sample(m7)
ESS[grep('beta', ESS[,1]),]
plot(density(ESS[grep('beta', ESS[,1]),2]))
# not bad! I mean, compared to CARBayesST, not bad. I like that.
ESS[is.nan(ESS[,2]),]
# So the NaNs might all be the elements of the Q matrix that are zero. That makes sense


# That's all good! What now? Missing data? Improving speed? 
# ok 

#### Now with missing data ####
df = df_miss
prep_stan_data_rushworth <- function(df, formula){
  N = nrow(df)
  N_T <- length(unique(df$date))
  N_F <- length(unique(df$facility))
  
  # order the data frame
  df <- df %>% arrange(date, facility)
  
  W_star = make_district_W2_matrix(df)
  
  # make the complete model matrix
  df2 <- df; df2$y[is.na(df2$y)] <- 0
  X = model.matrix(formula, data = df2)
  
  # get the outcome
  y = df$y
  
  # # comparing missingness.
  # if(!identical(as.integer(rownames(X_obs)), which(!is.na(df$y)))){
  #   stop('mismatch of model matrix and df missing rows')
  # }
  
  # missingness data
  N_miss = sum(is.na(y))
  N_obs = sum(!is.na(y))
  ind_miss = which(is.na(y))
  ind_obs = which(!is.na(y))
  y_obs = y[ind_obs]
  
  # compute the priors
  lm_fit <- glm(formula, family = 'poisson', data = df)
  coef_mat <- summary(lm_fit)$coefficients
  prior_mean_beta <- coef_mat[,1]
  sigma_beta = 10*vcov(lm_fit)
  
  # make the stan data frame
  stan_data <- list(
    N = N, # number of observations
    p = ncol(X_obs), # number of variables
    N_F = N_F, # number of facilities
    N_T = N_T, # number of time points
    N_miss = N_miss,
    N_obs = N_obs,
    ind_miss = ind_miss,
    ind_obs = ind_obs,
    X = X, # design matrix
    y_obs = y_obs, # outcome variable 
    mu = prior_mean_beta, # prior mean
    Sigma = sigma_beta, # prior variance
    W_star = W_star, 
    I = diag(1.0, N_F))
  
  return(stan_data)
}

stan_data <- prep_stan_data_rushworth(df, formula)

m8 <- stan(file = "regression_rushworth.stan", data = stan_data, 
           # Below are optional arguments
           iter = 2000, 
           warmup = 1000,
           chains = 1, 
           init = '0',
           cores = min(parallel::detectCores(), 4))

stan_est = extract(m8, pars = 'beta', permuted = F)
ESS <- effective_sample(m8)
ESS[grep('beta', ESS[,1]),]
plot(density(ESS[grep('beta', ESS[,1]),2]))


# BTW I'm fitting with the 2020 data, but that can be dealt with
res <- get_posterior_mean(m8)
betas <- res[grep('beta',rownames(res)),1]
names(betas) <- colnames(stan_data[['X']])
names(betas) <- gsub('\\(|\\)','',names(betas))

fac_beta_list <- list()
beta_ref <- betas[c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")]
for(f in  levels(df$facility)){
  # if this is the reference facility (likely A1 in my simulations)
  if(sum(grepl(f, names(betas))) == 0){
    beta_f <- beta_ref
  }else{
    cols <- paste0('facility', f, c("",":year", ":cos1", ":sin1", ":cos2", ":sin2", ":cos3", ":sin3"))
    beta_f <- betas[cols] + beta_ref 
    names(beta_f) <- c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")
  }
  fac_beta_list[[f]] <- beta_f
}

beta_df <- t(data.frame(fac_beta_list))

plot(beta_df[,1], lst$betas[,1], main = 'intercepts of facilities', xlab = 'mean posterior intercept', ylab = 'true intercept')
# ok good
plot(beta_df[,2], lst$betas[,2], main = 'year terms of facilities', xlab = 'mean posterior trend', ylab = 'true trend')
# ok also good

par(mfrow = c(3,2))
for(i in 3:8){
  plot(beta_df[,i], lst$betas[,i])
  print(cor(beta_df[,i], lst$betas[,i]))
}
# ok not bad. Unsurprising that these are more off.

#save(m8, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/m8_rushworth_missingdata_09012023.RData')
