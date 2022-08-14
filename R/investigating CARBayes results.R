### The aim here is to assess the results of CARBayes and the other models in ways
### that could be useful to identifying issues with the models
library(ggplot2)
library(dplyr)
library(cowplot)

setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')
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
  geom_line(aes(x = date, y = point_est, color = method)) + 
  facet_wrap(~facility) +
  ggtitle('median model fits across simulations for missing points') + 
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

