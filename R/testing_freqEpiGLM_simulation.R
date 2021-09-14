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

res <- simulate_data_freqGLM_epi(district_sizes = 4, R = 100, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = F)

df = res$df_list[[1]]

par(mfrow = c(2,2))
for(f in unique(df$facility)){
  tmp = df %>% filter(facility == f)
  plot(tmp$date, tmp$y, type = 'l', ylim = c(0,max(tmp$y)))
  lines(tmp$date, tmp$y_seasonal, col = 'red')
}



## parameter estimates
res_lst = list()

for(i in 1:100){
  df_miss = MCAR_sim(res$df_list[[i]], p = 0.2, by_facility = T)
  freqGLMepi_list = freqGLMepi_imputation(df_miss, prediction_intervals = 'bootstrap', R_PI = 2, verbose = F) 
  res_lst[[i]] = freqGLMepi_list
  #print(freqGLMepi_list$params)
}
# lots of iterations limits reached. Uh-oh
save(res_lst, res, file = 'results/freqGLM_epi_testFit_results_09072021.RData')


# ah no wonder the fitting is difficult. The y.neighbors term is highly replaceable with the beta terms. As in the y.neighbors term picks up the betas from the neighbors, which are highly correlated. So the seasonal values for these models are going to be highly correlated already and that will make this harder to converge.


### plotting the parameter estimates

# creating the necessary data frame, an ugly way
new_df = NULL
for(i in 1:length(res_lst)){
  new_df = rbind(new_df, t(res_lst[[i]]$params[,1,drop = F]))
}
new_df = as.data.frame(new_df)
params = tidyr::gather(new_df,  parameter, estimate, By.AR1:Bsin3)
params$parameter = gsub('By.AR1','lambda', params$parameter)
params$parameter = gsub('By.neighbors', 'alpha', params$parameter)

#params_true = as.data.frame(t(res$betas[1,,drop = F]))
#colnames(params_true) = 'estimate'
#params_true$parameter = paste0('B', rownames(params_true))

ggplot(params, aes(x = parameter, y = estimate)) + 
  geom_boxplot() + 
  geom_point(data = params_true, aes(x = parameter, y = estimate)) +
  ylim(-15, 15)



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





