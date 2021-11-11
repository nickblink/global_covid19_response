### Now doing this as an R script rather than Rmd because it's easier to work with.
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)


# Next to try to just fit the data without any missingness

# res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = F, seed = 10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

df = res$df_list[[1]]

par(mfrow = c(2,2))
for(f in unique(df$facility)){
  tmp = df %>% filter(facility == f)
  plot(tmp$date, tmp$y - tmp$y_seasonal, type = 'l', ylim = c(0,max(tmp$y)))
  lines(tmp$date, tmp$y_seasonal, col = 'red')
}

#### No missingness ####
res_lst = list()

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, scale_by_num_neighbors = T) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}


# save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_testFit_results_v2_10042021.RData')

#
#### No missingness - cheating with initial values ####
res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$phi, res$betas[i,])
}


for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits, scale_by_num_neighbors = T) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

#save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_testFit_results_10052021.RData')

#
#### MCAR 0.2 and parameter estimates ####

## parameter estimates
res_lst = list()

for(i in 1:100){
  df_miss = MCAR_sim(res$df_list[[i]], p = 0.2, by_facility = T)
  freqGLMepi_list = freqGLMepi_imputation(df_miss, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}
# lots of iterations limits reached. Uh-oh
# save(res_lst, res, file = 'results/freqGLM_epi_testFit_results_09232021.RData')


# ah no wonder the fitting is difficult. The y.neighbors term is highly replaceable with the beta terms. As in the y.neighbors term picks up the betas from the neighbors, which are highly correlated. So the seasonal values for these models are going to be highly correlated already and that will make this harder to converge.

# need a better way to explain this.

#### plotting the parameter estimates ####
# load('results/freqGLM_epi_testFit_results_09232021.RData')
# load('results/freqGLM_epi_noMISS_testFit_results_09232021.RData')
load('results/freqGLM_epi_noMISS_true_inits_testFit_results_09232021.RData')

warning('only looking at facility 1. Could look at other facilities too')

# creating the necessary data frame, an ugly way
params_true = as.data.frame(t(res$betas[1,,drop = F]))
colnames(params_true) = 'value'
params_true$parameter = paste0('B', rownames(params_true))
params_true = rbind(data.frame(value = c(res$lambda, res$phi), parameter = c('By.AR1', 'By.neighbors')),params_true)

new_df = NULL
for(i in 1:length(res_lst)){
  new_df = rbind(new_df, t(res_lst[[i]]$params[,1,drop = F]))
}

# should be true as a check
identical(colnames(new_df), params_true$parameter)

# organize the parameters data frmae
params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
params = merge(params, params_true)
params$residual = params$estimate - params$value
params$residual_prop = params$residual/abs(params$value)
params$parameter = gsub('By.AR1','lambda', params$parameter)
params$parameter = gsub('By.neighbors', 'phi', params$parameter)

# params_res = tidyr::gather(as.data.frame(resid_df),  parameter, residual, By.AR1:Bsin3)
# test = merge(params_res, params_true)
# 
# params_res$parameter = gsub('By.AR1','lambda', params_res$parameter)
# params_res$parameter = gsub('By.neighbors', 'phi', params_res$parameter)

# ggplot(params, aes(x = parameter, y = estimate)) + 
#   geom_boxplot() # + 
#   geom_point(data = params_true, aes(x = parameter, y = estimate)) +
#   ylim(-15, 15)

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')

# Good plot to see. So something is off with all of them,
# but especially with lambda, which is surprising because I would expect phi to be the most off. Why is the autoregressive term so bad? I will need to explore more.


### Plotting the fits
tmp = res$df_list[[1]]
plot(res$df_list[[1]]$y, type = 'l')
lines(res_lst[[1]]$df$y_imp, col = 'blue')

par(mfrow = c(2,2))
for(f in unique(df$facility)){
  tmp = res_lst[[1]]$df %>% arrange(date) %>% filter(facility == f)
  plot(tmp$date, tmp$y_true, type = 'l')
  lines(tmp$date, tmp$y_pred_freqGLMepi, col = 'red')
}




#### No missingness and parameter estimates ####

load('results/freqGLM_epi_noMISS_true_inits_testFit_results_09232021.RData')

#par_true = res$
  
# creating the necessary data frame, an ugly way
params_true = as.data.frame(t(res$betas))
rownames(params_true) = paste0('B', rownames(params_true))
params_true = rbind(t(data.frame(By.AR1 = rep(res$lambda, ncol(params_true)), By.neighbors = rep(res$phi, ncol(params_true)), row.names = colnames(params_true))), params_true)

t(data.frame(By.AR1 = rep(res$lambda, 4), By.neighbors = rep(res$phi, 4)))

facs = c('A1','A2','A3','A4')  

D = NULL #as.data.frame(param_norm = NULL, prediction_norm = NULL, iter = NULL)
for(i in 1:length(res_lst)){
  # compare the norms of the estimates
  par_est = res_lst[[i]]$params
  df = res_lst[[i]]$df
  
  # norm_off = norm(as.matrix(par_est - params_true))/norm(as.matrix(params_true))
  # diff = df$y_true - df$y_pred_freqGLMepi
  # pred_off = norm(na.omit(diff))/norm(na.omit(df$y_true))
  # 
  # new_row = data.frame(param_norm = norm_off, prediction_norm = pred_off, iter = i)
  # D = rbind(D, new_row)

  for(j in 1:4){
    norm_off = norm(as.matrix(par_est[,j] - params_true[,j]))/norm(as.matrix(params_true[,j]))
    
    tmp = df %>% filter(facility == facs[j])
    diff = tmp$y_true - tmp$y_pred_freqGLMepi
    pred_off = norm(na.omit(diff))/norm(na.omit(tmp$y_true))
    convergence = res_lst[[i]]$convergence[j]
    
    new_row = data.frame(param_norm = norm_off, prediction_norm = pred_off, iter = i, facility = facs[j], convergence = convergence)
    D = rbind(D, new_row)
  }
  
  # compare the normalized norms (?)
  # After, investigate the difference in the parameters
}

lm.fit = lm(D, formula = log(prediction_norm) ~ log(param_norm))

ggplot(data = D, aes(x = log(param_norm), y = log(prediction_norm))) +
  geom_point() + 
  xlab('log(proportion norm off of parameters)') + 
  ylab('log(proportion norm off of prediction)') + 
  geom_abline(intercept = lm.fit$coefficients[1], slope = lm.fit$coefficients[2], col = 'red')

# hm ok so there actually is some relationship. Why were these picked? Are the mean values of the true parameters ever better? That should be an issue. I can't imagine because they can't be better in terms of the log likelihood, otherwise they would have been picked


res_lst[[40]]$params
test = res_lst[[40]]$df

# why is it always A1 that's the most off? The higher intercept? Or is there some weird bug in how I simulate the data that makes it worse. 



#### Testing sample sizes - 20 years ####

res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$phi, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

# save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_20yrs_testFit_results_10072021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = c(res$lambda, res$phi), parameter = c('By.AR1', 'By.neighbors')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  params$parameter = gsub('By.AR1','lambda', params$parameter)
  params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#### Testing sample sizes - 100 years ####

res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '1970-01-01', b1_mean = -.01, b1_sd = .001)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$phi, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

HERE AT 312pm on 10/01

TRY THIS WITH LOWER LAMBDA AND PHI VALUES? THAT MIGHT MAKE IT CONVERGE BETTER

# save(res_lst, res, file = 'results/freqGLM_epi_noMISS_true_inits_100yrs_testFit_results_10012021.RData')


# creating the necessary data frame, an ugly way
params_true = as.data.frame(t(res$betas[1,,drop = F]))
colnames(params_true) = 'value'
params_true$parameter = paste0('B', rownames(params_true))
params_true = rbind(data.frame(value = c(res$lambda, res$phi), parameter = c('By.AR1', 'By.neighbors')),params_true)

new_df = NULL
for(i in 1:length(res_lst)){
  new_df = rbind(new_df, t(res_lst[[i]]$params[,1,drop = F]))
}

# should be true as a check
identical(colnames(new_df), params_true$parameter)

# organize the parameters data frmae
params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
params = merge(params, params_true)
params$residual = params$estimate - params$value
params$residual_prop = params$residual/abs(params$value)
params$parameter = gsub('By.AR1','lambda', params$parameter)
params$parameter = gsub('By.neighbors', 'phi', params$parameter)


ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#### Testing no lambda (no autocorrelation) ####
source('R/freqGLM_epi_fxns_noLambda.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = 10, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2016-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$phi, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits, scale_by_num_neighbors = T) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}



# save(res_lst, res, file = 'results/freqGLM_epi_noMISS_true_inits_noLambda_testFit_results_10072021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$phi, parameter = c('By.neighbors')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.neighbors:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  #params$parameter = gsub('By.AR1','lambda', params$parameter)
  params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')

#
#### Testing no lambda (no autocorrelation) - 20 years ####
source('R/freqGLM_epi_fxns_noLambda.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = 10, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$phi, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

#save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_noLambda_20yrs_testFit_results_10072021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$phi, parameter = c('By.neighbors')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.neighbors:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  #params$parameter = gsub('By.AR1','lambda', params$parameter)
  params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')



#### Testing no lambda (no autocorrelation) - 20 years - higher spatial ####
source('R/freqGLM_epi_fxns_noLambda.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = 10, phi = -1.3, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$phi, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits, scale_by_num_neighbors = T) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

#save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_noLambda_phi1.3_20yrs_testFit_results_10072021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$phi, parameter = c('By.neighbors')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.neighbors:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  #params$parameter = gsub('By.AR1','lambda', params$parameter)
  params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#### Testing no phi (no spatial correlation) - 4 years ####
source('R/freqGLM_epi_fxns_noPhi.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = 10, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2016-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}


#save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_nophi_testFit_results_10012021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$lambda, parameter = c('By.AR1')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  params$parameter = gsub('By.AR1','lambda', params$parameter)
  #params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#
#### Testing no phi (no spatial correlation) - 20 years ####
source('R/freqGLM_epi_fxns_noPhi.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = 10, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}


# save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_nophi_20yrs_testFit_results_10012021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$lambda, parameter = c('By.AR1')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  params$parameter = gsub('By.AR1','lambda', params$parameter)
  #params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#
#### Testing no phi (no spatial correlation) - 20 years - higher lambda ####
source('R/freqGLM_epi_fxns_noPhi.R')
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -1.3, phi = 10, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01', b1_mean = -0.1, b1_sd = 0.1)

res_lst = list()

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$betas[i,])
}

for(i in 1:100){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}


# save(res_lst, res, file = 'results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_nophi_lambda13_20yrs_testFit_results_10012021.RData')

params_full = NULL
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = res$lambda, parameter = c('By.AR1')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  params$parameter = gsub('By.AR1','lambda', params$parameter)
  #params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

ggplot(params, aes(x = parameter, y = residual)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  ggtitle('parameter estimate residuals')

ggplot(params, aes(x = parameter, y = residual_prop)) + 
  geom_boxplot() +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_hline(yintercept = -1, color = 'red') +
  ggtitle('parameter estimate residuals/abs(true value)')


#
#### Plotting ####

load('results/freqGLM_epi_debugging/freqGLM_epi_noMISS_true_inits_testFit_results_10052021.RData')


# for the model with lambda and phi
params_full = NULL
param_norm_diff = data.frame(A1 = rep(as.numeric(NA), 100),
                             A2 = rep(as.numeric(NA), 100),
                             A3 = rep(as.numeric(NA), 100),
                             A4 = rep(as.numeric(NA), 100))
for(j in 1:4){
  # creating the necessary data frame, an ugly way
  params_true = as.data.frame(t(res$betas[j,,drop = F]))
  colnames(params_true) = 'value'
  params_true$parameter = paste0('B', rownames(params_true))
  params_true = rbind(data.frame(value = c(res$lambda, res$phi), parameter = c('By.AR1', 'By.neighbors')),params_true)
  
  new_df = NULL
  for(i in 1:length(res_lst)){
    new_df = rbind(new_df, t(res_lst[[i]]$params[,j,drop = F]))
  }
  
  # should be true as a check
  print(identical(colnames(new_df), params_true$parameter))
  
  param_norm_diff[,j] = apply(new_df, 1, function(xx){
    return(norm(xx - params_true$value))
  })
  
  # organize the parameters data frmae
  params = tidyr::gather(as.data.frame(new_df),  parameter, estimate, By.AR1:Bsin3)
  params = merge(params, params_true)
  params$residual = params$estimate - params$value
  params$residual_prop = params$residual/abs(params$value)
  params$parameter = gsub('By.AR1','lambda', params$parameter)
  params$parameter = gsub('By.neighbors', 'phi', params$parameter)
  
  params$facility = i
  
  params_full = rbind(params_full, params)
}

param_norm_diff$sum = rowSums(param_norm_diff)

indices = c()
tt = param_norm_diff$sum
for(i in 1:10){
  ind = which(tt == max(tt, na.rm = T))
  tt[ind] = NA
  indices = c(indices, ind)
}

par(mfrow = c(2,2))
for(ii in indices){
  df = res_lst[[ii]]$df
  print(res_lst[[ii]]$params)
  for(f in unique(df$facility)){
    tmp = df %>% filter(facility == f)
    plot(tmp$date, tmp$y, type = 'l', ylim = c(0,max(tmp$y)), main = paste0(ii, ': ',f))
    lines(tmp$date, tmp$y_pred_freqGLMepi, col = 'red')
  }
}


#### Testing parametric bootstrap ####

res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 2, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

tmp = res$df_list[[1]]
tmp$y_true = tmp$y
freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 2, scale_by_num_neighbors = T) 

A = tmp %>% filter(facility == 'A1')

## Try this with 20 years
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 2, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01')

tmp = res$df_list[[1]]
tmp$y_true = tmp$y
freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 2, scale_by_num_neighbors = T) 


### Counting the number of errors
set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

num_errors = sum(sapply(res_lst, function(xx) xx$num_errors))
# 116
prop_erros = num_errors/400
# 29%!

# Now with 20 years
set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01')

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

num_errors = sum(sapply(res_lst, function(xx) xx$num_errors)) # 71
prop_erros = num_errors/400 # 18%!

#
model.cols = c(16, 3, 5, 7, 8, 9, 10, 11, 12)



corrplot::corrplot(cor(tmp[tmp$facility == 'A1', model.cols], use = 'complete.obs'))

#### Testing parametric bootstrap - cheating with initial values ####

### Counting the number of errors
set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

optim_inits = list()
for(i in 1:nrow(res$betas)){
  optim_inits[[rownames(res$betas)[i]]] = c(res$lambda, res$phi, res$betas[i,])
}

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

sum(sapply(res_lst, function(xx) xx$num_errors))
# 61

sum(sapply(res_lst, function(xx) xx$num_errors))/400
# .15/ 15% - so significantly better

# Now with 20 years
set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, start_date = '2000-01-01')

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

num_errors = sum(sapply(res_lst, function(xx) xx$num_errors)) 
prop_erros = num_errors/400 

#
model.cols = c(16, 3, 5, 7, 8, 9, 10, 11, 12)



corrplot::corrplot(cor(tmp[tmp$facility == 'A1', model.cols], use = 'complete.obs'))

#### Testing parametric bootstrap on diff. models ####

source('R/freqGLM_epi_fxns_noLambda.R')

set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 5, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10, )

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

# Still not good. Are the y's between facilities just too highly correlated?
cors = c()
for(i in 1:length(res$df_list)){
  tmp = res$df_list[[i]]
  tt = tmp %>% 
    select(date, facility, y) %>%
    reshape(., idvar = 'date', timevar = 'facility', direction = 'wide')
  
  cors = c(cors, c(cor(tt[,-1])))
}
cors = cors[cors != 1]

# Darn should I mess around with the data?

## No phi
source('R/freqGLM_epi_fxns_noPhi.R')

set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 5, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

# Still not good. Are the y's between facilities just too highly correlated?
cors = c()
for(i in 1:length(res$df_list)){
  tmp = res$df_list[[i]]
  tt = tmp %>% 
    select(date, facility, y) %>%
    reshape(., idvar = 'date', timevar = 'facility', direction = 'wide')
  
  cors = c(cors, c(cor(tt[,-1])))
}
cors = cors[cors != 1]

#### Testing parametric bootstrap no lambda or phi ####
source('R/freqGLM_epi_fxns_noPhi_noLambda.R')

set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

sum(sapply(res_lst, function(xx) xx$num_errors))
# phew. But also quite low covariances

tt = unlist(lapply(res_lst, function(xx){
  # print(unique(xx$df$district));
  lapply(xx$vcov_list, function(yy){
    #browser()
    res<- ifelse(is.null(yy), -1, det(yy));
    return(res)
  })}))

mean(tt == -1)
mean(tt)
median(tt)
# so very very low. Is that a problem? It's the exponentiating that is potentially an issue here.




#### Testing parametric bootstrap on real data ####

D = readRDS('data/liberia_cleaned_NL.rds')

facility_miss = D %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            ari_miss_count = sum(is.na(indicator_count_ari_total)))

D = D %>% 
  filter(facility %in% facility_miss$facility[facility_miss$ari_miss <= 0.5])

# fix name for smoother name matching
D$facility = gsub("'","", D$facility)
D$facility = gsub("-"," ", D$facility)

tt = table(D$district)
multi_dists = names(tt[tt > 56])

res_lst = list()
for(dist in multi_dists){
  # prep the data
  tmp = D %>% filter(district == dist) %>%
    select(date, district, facility, y = indicator_count_ari_total)
  tmp = tmp %>% 
    add_periodic_cov() %>%
    add_autoregressive() %>%
    add_neighbors()
  tmp$y_true = tmp$y
  
  # run the model
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[dist]] = freqGLMepi_list
}

lst = unlist(lapply(res_lst, function(xx) xx$convergence))
table(lst)

tt = unlist(lapply(res_lst, function(xx){
 # print(unique(xx$df$district));
  lapply(xx$vcov_list, function(yy){
    #browser()
    res<- ifelse(is.null(yy), -1, det(yy));
    return(res)
  })}))

mean(tt == -1)



lapply(res_lst, function(xx){print(length(xx$vcov_list))})
#### Testing parametric bootstrap with nlminb ####
source('R/freqGLM_epi_fxns_expPhiLamVar3.R')

set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

sum(sapply(res_lst, function(xx) xx$num_errors)) # Only 1 Out of 400!

# for testing
tt <- sapply(res_lst, function(xx) xx$params[,1])

# a little bias. That could have to do with the initial values

save(res_lst, file = 'results/freGLM_100sim_nlminb.RData')

### Want to compare some of the variances when it does converge
source('R/imputation_functions.R')
res_lst2 = list()
for(i in 1:10){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst2[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

diag(res_lst2[[1]]$vcov_list[[2]])
diag(res_lst[[1]]$vcov_list[[2]])
# oh ya we wouldn't expect these to be the exact same because of the exponentiating duh. But shouldn't the other diagonals be the same?

# Should I try rerunning this all with nlminb?

#### Testing parametric bootstrap with exponentiated first two params ####
source('R/freqGLM_epi_fxns_expPhiLamVar.R')

set.seed(10)
res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = T, seed = 10)

res_lst = list()
for(i in 1:length(res$df_list)){
  #df_miss = MCAR_sim(, p = 0.2, by_facility = T)
  tmp = res$df_list[[i]]
  tmp$y_true = tmp$y
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'parametric_bootstrap', R_PI = 10, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

sum(sapply(res_lst, function(xx) xx$num_errors)) # only 16/400? What's up with that? So maybe it does have to do with finding the inverse of a high determinant matrix

#### Testing the fits of other packages ####
# addreg approach - nvm each x needs to be non-negative this doesn't work. Also this doesn't really fit what my model is doing anyway. I'm doing some weird function of the predictors that's not even linear
y = A$y
x = A[,c('y.neighbors', 'y.AR1', 'year', "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")]

fit.1 <- addreg::nnpois(y, x)






















