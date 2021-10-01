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

res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = F, seed = 10)

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
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

# save(res_lst, res, file = 'results/freqGLM_epi_noMISS_true_inits_testFit_results_09232021.RData')

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
  freqGLMepi_list = freqGLMepi_imputation(tmp, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F, optim_init = optim_inits) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}

save(res_lst, res, file = 'results/freqGLM_epi_noMISS_true_inits_testFit_results_09232021.RData')

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

