### The aim here is to assess the results of CARBayes and the other models in ways
### that could be useful to identifying issues with the models
library(ggplot2)
library(dplyr)
library(cowplot)
library(CARBayesST)
library(spdep)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('..')

source('R/imputation_functions.R')

#
#### freqEpi DGP analysis ####
load('results/simulation_epi_MCARp2_R500_res_07142022.RData')

# prep data
{
imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CB_intercept", "y_CB_facility")
rename_vec = c('glmFreq','glmFreq_epi','CARBayes_int', 'CARBayes_facility')
color_vec = c('red','blue','lightgreen', 'forestgreen')

df <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, rm_ARna = T, use_point_est = F, imputed_only = T)

# replace the names
if(!is.null(rename_vec)){
  for(i in 1:length(imp_vec)){
    df$method = gsub(imp_vec[i], rename_vec[i],df$method)
  }
}

df$method = factor(df$method, levels = rename_vec)
}
  
# plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])
  
# (0) Plot some fits from different runs
plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[3]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[4]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

# (1)
# plot all models' mean/median estimates vs truth
ggplot(data = df) +
  geom_line(aes(x = date, y = y_missing)) + 
  geom_line(aes(x = date, y = median, color = method)) + 
  facet_wrap(~facility) +
  ggtitle('mean model fits across simulations for missing points') + 
  ylab('y')

# so it looks like on average the models do pretty well. This matches up with the bias being relatively low to the true numbers in what I am running. Now, what's going on with the prediction intervals?

# (2) 
# plot the prediction intervals of methods vs truth

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = y_true)) + 
  geom_line(aes(x = date, y = point_est, color = method)) + 
  geom_ribbon(aes(x = date, ymin = lower_025, ymax = upper_975, fill = method, colour = method), alpha = 0.1) + 
  facet_wrap(~facility) + 
  ggtitle('Median upper 97.5% and 2.5% prediction quantiles for R = 500 simulations')
## OK this doesn't really make sense to look at. Of course it's different than the interval width. It's just showing that there's not really systematic bias in the upper and lower prediction intervals. It's not showing that they can get super wide.

# (3) 
# for a given metric 
# and a given model
# plot the metric over time for each facility
# use facets!
ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = bias, color = method)) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~facility) +
  ggtitle('bias (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = RMSE, color = method)) + 
  facet_wrap(~facility) +
  ggtitle('RMSE (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = interval_width, color = method)) + 
  facet_wrap(~facility) +
  ggtitle('interval width (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = prop_interval_width, color = method)) + 
  facet_wrap(~facility) +
  ggtitle('proportional interval width (imputed points only)')


ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = coverage95, color = method)) + 
  geom_hline(yintercept = 0.95) + 
  facet_wrap(~facility) +
  ggtitle('coverage95 (imputed points only)')


#### CARBayes DGP Analysis ####
load('results/simulation_ST_MCARp2_R500_res_08092022.RData')

# prep data
{
  imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CB_intercept", "y_CB_facility")
  rename_vec = c('glmFreq','glmFreq_epi','CARBayes_int', 'CARBayes_facility')
  color_vec = c('red','blue','lightgreen', 'forestgreen')
  
  df <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, rm_ARna = T, use_point_est = F, imputed_only = T)
  
  # replace the names
  if(!is.null(rename_vec)){
    for(i in 1:length(imp_vec)){
      df$method = gsub(imp_vec[i], rename_vec[i],df$method)
    }
  }
  
  df$method = factor(df$method, levels = rename_vec)
}

# plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

# (0) Plot some fits from different runs
plot_facility_fits(imputed_list[[1]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[3]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

plot_facility_fits(imputed_list[[4]], imp_vec = imp_vec[c(1,2,4)], imp_names = rename_vec[c(1,2,4)], color_vec = color_vec[c(1,2,4)])

# (1)
# plot all models' mean/median estimates vs truth
ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = median, color = method), size = 1) + 
  geom_line(aes(x = date, y = y_missing), size = 1) + 
  facet_wrap(~facility) +
  ggtitle('mean model fits across simulations for missing points') + 
  ylab('y')

# so it looks like on average the models do pretty well. This matches up with the bias being relatively low to the true numbers in what I am running. Now, what's going on with the prediction intervals?

# (2) 
# plot the prediction intervals of methods vs truth

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = y_true)) + 
  geom_line(aes(x = date, y = point_est, color = method)) + 
  geom_ribbon(aes(x = date, ymin = lower_025, ymax = upper_975, fill = method, colour = method), alpha = 0.1) + 
  facet_wrap(~facility) + 
  ggtitle('Median upper 97.5% and 2.5% prediction quantiles for R = 500 simulations')
## OK this doesn't really make sense to look at. Of course it's different than the interval width. It's just showing that there's not really systematic bias in the upper and lower prediction intervals. It's not showing that they can get super wide.

# (3) 
# for a given metric 
# and a given model
# plot the metric over time for each facility
# use facets!
ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = bias, color = method), size = 1) + 
  geom_hline(yintercept = 0) + 
  facet_wrap(~facility) +
  ggtitle('bias (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = RMSE, color = method), size = 1) + 
  facet_wrap(~facility) +
  ggtitle('RMSE (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = interval_width, color = method), size = 1) + 
  ylim(0,2500) + 
  facet_wrap(~facility) +
  ggtitle('interval width (imputed points only)')

ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = prop_interval_width, color = method), size = 1) + 
  facet_wrap(~facility) +
  ggtitle('proportional interval width (imputed points only)')


ggplot(data = df %>% filter(method != 'CARBayes_int')) +
  geom_line(aes(x = date, y = coverage95, color = method),  size = 1) + 
  geom_hline(yintercept = 0.95) + 
  facet_wrap(~facility) +
  ggtitle('coverage95 (imputed points only)')

#### Investigating Convergence ####

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 5, rho = 0.5, alpha = 0.3, tau = 0.5)

par(mfrow = c(3,3))

for(i in 1:3){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
  
  CAR_list3 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')
  
  mc <- CAR_list3$model_chain
  
  rho <- mc$samples$rho
  rho_means <- data.frame(rho.S = rep(NA, 500), rho.T = rep(NA, 500))
  for(j in 1:2){
    rho_means[,j] <- sapply(1:nrow(rho), function(ii) mean(rho[1:ii,j]))
  }
  
  tau2 <- mc$samples$tau2
  tau2_means <- sapply(1:length(tau2), function(ii) mean(tau2[1:ii]))
  
  plot((1:500)*10,rho_means[,1], main = 'rolling average of spatial parameter (true is 0.5)', xlab = 'iteration')
  
  plot((1:500)*10,rho_means[,2], main = 'rolling average of temporal parameter (true is 0.3)', xlab = 'iteration')
  
  plot((1:500)*10,tau2_means, main = 'rolling average of R.E. variance parameter (true is 0.25)', xlab = 'iteration')

}

# ggplot() + 
#   geom_abline(intercept = 0.25, slope = 0, col = 'red') + 
#   geom_line(aes(y = tau2_means, x = 1:500)) + ylim(0,0.5)


beta <- mc$samples$beta
beta_means <- as.data.frame(matrix(NA, nrow = nrow(beta), ncol = ncol(beta)))
for(j in 1:ncol(beta)){
  beta_means[,j] <- sapply(1:nrow(beta), function(ii) mean(beta[1:ii,j]))
}

#### Investigating Convergence - higher sample size ####

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 5, rho = 0.5, alpha = 0.3, tau = 0.5)

par(mfrow = c(3,3))

for(i in 1:3){
  df = lst$df_list[[i]]
  
  # simulation function!
  df_miss = MCAR_sim(df, p = 0.2, by_facility = T)
  
  CAR_list3 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 5000, n.sample = 100000, prediction_sample = T, model = 'facility_fixed')
  
  mc <- CAR_list3$model_chain
  
  tau2 <- mc$samples$tau2
  tau2_means <- sapply(1:length(tau2), function(ii) mean(tau2[1:ii]))
  
  rho <- mc$samples$rho
  rho_means <- data.frame(rho.S = rep(NA, length(tau2)), rho.T = rep(NA, length(tau2)))
  for(j in 1:2){
    rho_means[,j] <- sapply(1:nrow(rho), function(ii) mean(rho[1:ii,j]))
  }
  

  plot((1:length(tau2))*10,rho_means[,1], main = 'rolling average of spatial parameter (true is 0.5)', xlab = 'iteration')
  
  plot((1:length(tau2))*10,rho_means[,2], main = 'rolling average of temporal parameter (true is 0.3)', xlab = 'iteration')
  
  plot((1:length(tau2))*10,tau2_means, main = 'rolling average of R.E. variance parameter (true is 0.25)', xlab = 'iteration')
  
}

HERE YO
# ggplot() + 
#   geom_abline(intercept = 0.25, slope = 0, col = 'red') + 
#   geom_line(aes(y = tau2_means, x = 1:500)) + ylim(0,0.5)


beta <- mc$samples$beta
beta_means <- as.data.frame(matrix(NA, nrow = nrow(beta), ncol = ncol(beta)))
for(j in 1:ncol(beta)){
  beta_means[,j] <- sapply(1:nrow(beta), function(ii) mean(beta[1:ii,j]))
}


#### Investigating MCMC properties ####
source('R/personal_CARBayes_fxn.R')

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 5, rho = 0.5, alpha = 0.3, tau = 0.5)

df = lst$df_list[[1]]

# simulation function!
df_miss = MCAR_sim(df, p = 0.2, by_facility = T)

CAR_list1 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 500, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')

folder <- '../CARBayesST_personal/R'

for(f in list.files(folder, full.names = T)){
  print(f)
  source(f)
}
setwd('../CARBayesST_personal/src')
Rcpp::sourceCpp('CARBayesST.cpp')
Rcpp::sourceCpp('RcppExports.cpp')

CAR_list2 = CARBayes_imputation(df_miss, col = "y", return_type = 'all', burnin = 500, n.sample = 10000, prediction_sample = T, model = 'facility_fixed')

## DOESNT WORK. WHY?