### The aim here is to assess the results of CARBayes and the other models in ways
### that could be useful to identifying issues with the models
library(ggplot2)
library(dplyr)


setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')
load('results/simulation_epi_MCARp2_R500_res_07142022.RData')

df <- calculate_metrics_by_point(imputed_list, imp_vec = c("y_pred_harmonic", "y_pred_freqGLMepi", "y_CB_intercept", "y_CB_facility"), rm_ARna = T, use_point_est = F)

# (1)
# plot all models' mean/median estimates vs truth
ggplot(data = df %>% filter(method == "y_CB_facility")) +
  geom_line(aes(x = date, y = y_missing)) + 
  geom_line(aes(x = date, y = point_est, color = 'red')) + 
  facet_wrap(~facility)

# so it looks like on average the models do pretty well. This matches up with the bias being relatively low to the true numbers in what I am running. Now, what's going on with the prediction intervals?

# (2) 
# plot the prediction intervals of methods vs truth

ggplot(data = df) +
  geom_line(aes(x = date, y = y_missing)) + 
  geom_line(aes(x = date, y = point_est, color = method)) + 
  facet_wrap(~facility)

p1 <- ggplot() +
  geom_line(data = tmp, aes(x = date, y = y_true), size = 1) +
  geom_line(data = df_f, aes(x = date, y = y, group = method, color = method)) +
  geom_ribbon(data = df_f, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.1) +
  scale_color_manual(values = c(color_vec)) + 
  scale_fill_manual(values = c(color_vec)) + 
  ylim(c(0,1.5*max(tmp$y_true))) + 
  ggtitle(sprintf('facility %s', f)) + 
  ylab('y') +
  theme_bw() +
  theme(text = element_text(size = 10))

# (2.5) plot all models' mean intervals vs truth

# (3) 
# for a given metric 
# and a given model
# plot the metric over time for each facility
# use facets!