current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')
library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)


#### 12/16/2022: Plotting across methods (lost work) ####
file_MCAR <- grep('mcar', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

plot_all_methods(file_MCAR[1:7])


#
#### 12/14/2022: Analyzing Outbreak Detection results ####
# file_MCAR <- grep('mcar', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MCAR <- grep('mcar06_nost_beta6', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MCAR <- NULL
for(file in file_MCAR){
  p <- stringr::str_match(file, 'mcar(.*?)_')[[2]]
  print(p)
  print(file)
  
  # if a directory, combine results. If already combined, load them
  if(dir.exists(file)){
    lst_full <- combine_results(input_folder = file, return_lst = T, results_file = NULL)
  }else{
    load(file)
  }
  
  HERE AT 321pm - NEED TO COMBINE THE RESULTS PROPERLY GIVEN THE NEW FORMAT
  
  #res <- calculate_metrics_by_point(lst_full, imp_vec =  c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 
  res$method = paste0(res$method, sprintf('_MCAR_p%s_', p))
  res_MCAR <- rbind(res_MCAR, res)
}

### WF MCAR
res_WF <- res_MCAR[grep('CCA_WF', res_MCAR$method),]
imp_vec <- unique(res_WF$method)
rename_vec <- c('prop=0','prop=0.1','prop=0.2','prop=0.3','prop=0.4','prop=0.5','prop=0.6')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue2', 'blue4')

test = res_WF %>% 
  filter(date == '2020-01-01',
         method == 'y_pred_CCA_WF_MCAR_p06_')

test2 = res_WF %>% 
  filter(method == 'y_pred_CCA_WF_MCAR_p06_')

cor(test2$y_exp,
    test2$outbreak_detection3)
# -0.009

plot(log(test2$y_exp),
    test2$outbreak_detection3)
# -0.12

plot(test2$y_exp,
     test2$outbreak_detection3)

plot(test$y_exp,
     test$outbreak_detection3)

plot(test$y_exp,
    log(test$outbreak_detection3))

cor(test$y_exp,
     test$outbreak_detection3)
# -0.12

cor(log(test$y_exp),
    test$outbreak_detection3)

cor(log(test$y_exp),
    log(test$outbreak_detection3))

cor(test$y_exp,
    log(test$outbreak_detection3))


#### 11/21/2022: Plotting across missingness values ####
file_MCAR <- grep('mcar', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MCAR <- NULL
for(file in file_MCAR){
  p <- stringr::str_match(file, 'mcar(.*?)_')[[2]]
  print(p)
  print(file)
  
  # if a directory, combine results. If already combined, load them
  if(dir.exists(file)){
    lst_full <- combine_results(input_folder = file, return_lst = T, results_file = NULL)
  }else{
    load(file)
  }
  res <- calculate_metrics_by_point(lst_full, imp_vec =  c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 
  res$method = paste0(res$method, sprintf('_MCAR_p%s_', p))
  res_MCAR <- rbind(res_MCAR, res)
}

### WF MCAR
res_WF <- res_MCAR[grep('CCA_WF', res_MCAR$method),]
imp_vec <- unique(res_WF$method)
rename_vec <- c('prop=0','prop=0.1','prop=0.2','prop=0.3','prop=0.4','prop=0.5','prop=0.6')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue2', 'blue4')

plot_metrics_by_point(res = res_WF, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec, min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))


#
### CAR MCAR
res_CAR <- res_MCAR[grep('CCA_CAR', res_MCAR$method),]
imp_vec <- unique(res_CAR$method)
rename_vec <- c('prop=0','prop=0.1','prop=0.2','prop=0.3','prop=0.4','prop=0.5','prop=0.6')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue2', 'blue4')

plot_metrics_by_point(res = res_CAR, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec, min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))

### freqEpi MCAR
res_freq <- res_MCAR[grep('CCA_freqGLMepi', res_MCAR$method),]
imp_vec <- unique(res_freq$method)
rename_vec <- c('prop=0','prop=0.1','prop=0.2','prop=0.3','prop=0.4','prop=0.5','prop=0.6')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue2', 'blue4')

plot_metrics_by_point(res = res_freq, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec, min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))

#### MNAR ####
file_MNAR <- grep('simulation_noST_mnar.*11022022', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MNAR <- NULL
for(file in file_MNAR){
  p <- stringr::str_match(file, 'mnarp(.*?)_')[[2]]
  print(file)
  load(file)
  res <- calculate_metrics_by_point(lst_full, imp_vec =  c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 
  res$method = paste0(res$method, sprintf('_MNAR_p%s', p))
  res_MNAR <- rbind(res_MNAR, res)
}

### WF MNAR
res_WF <- res_MNAR[grep('CCA_WF|baseline_WF_MNAR_p1', res_MNAR$method),]
imp_vec <- unique(res_WF$method)
rename_vec <- c('p=0','p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_WF, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width'))

### CAR MNAR
res_CAR <- res_MNAR[grep('CCA_CAR', res_MNAR$method),]
imp_vec <- unique(res_CAR$method)
rename_vec <- c('p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_CAR, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width'))

### freqEpi MNAR
res_freq <- res_MNAR[grep('CCA_freqGLMepi', res_MNAR$method),]
imp_vec <- unique(res_freq$method)
rename_vec <- c('p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_freq, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0, max_intW_lim = 80,  metric_list = c('bias','RMSE','coverage95','interval_width'))


#### MAR ####
file_MAR <- grep('simulation_noST_mar.*11022022', dir('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MAR <- NULL
for(file in file_MAR){
  p <- stringr::str_match(file, 'marp(.*?)_')[[2]]
  print(file)
  load(file)
  res <- calculate_metrics_by_point(lst_full, imp_vec =  c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 
  res$method = paste0(res$method, sprintf('_MAR_p%s', p))
  res_MAR <- rbind(res_MAR, res)
}

### WF MAR
res_WF <- res_MAR[grep('CCA_WF|baseline_WF_MAR_p1', res_MAR$method),]
imp_vec <- unique(res_WF$method)
rename_vec <- c('p=0','p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('black','deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_WF, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width'))

### CAR MNAR
res_CAR <- res_MAR[grep('CCA_CAR', res_MAR$method),]
imp_vec <- unique(res_CAR$method)
rename_vec <- c('p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_CAR, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0,  metric_list = c('bias','RMSE','coverage95','interval_width'))

### freqEpi MNAR
res_freq <- res_MAR[grep('CCA_freqGLMepi', res_MAR$method),]
imp_vec <- unique(res_freq$method)
rename_vec <- c('p=0.1','p=0.2','p=0.3','p=0.4','p=0.5')
color_vec <-c('deepskyblue','deepskyblue3','deepskyblue4','blue1','blue4')

plot_metrics_by_point(res = res_freq, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', rename_vec = rename_vec,
                      min_missing = 0, max_intW_lim = 80,  metric_list = c('bias','RMSE','coverage95','interval_width'))

#### Plotting CI's ####
head(res_CAR)

plot_list = list()
i = 0
for(metric in c('bias', 'RMSE', 'coverage95', 'interval_width','outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10')){
  i = i + 1
  tmp <- res_MCAR %>%
  group_by(method) %>% 
  summarize(median = median(get(metric)),
            lower = stats::quantile(get(metric), probs = 0.025),
            upper = stats::quantile(get(metric), probs = 0.975))
  str_lst <- strsplit(tmp$method, '_p0')
  
  tmp$method = unlist(lapply(str_lst, '[[', 1))
  tmp$p = as.numeric(gsub('_','',unlist(lapply(str_lst, '[[', 2))))/10
  tmp$p[is.na(tmp$p)] <- 0
  tmp$method <- gsub('y_pred_|_MCAR', '', tmp$method)
  
  p1 <- ggplot(tmp, aes(x = p, y = median, ymin = lower, ymax = upper, color = method)) + 
    geom_point(position = position_dodge(width = 0.1)) + 
    geom_errorbar(position = position_dodge(width = 0.1)) +
    ylab(metric) +
    theme_bw()
  
  legend = get_legend(p1 + theme(legend.position = 'bottom'))
  p1 <- p1 + theme(legend.position = 'none')
  
  plot_list[[i]] <- p1
}

final_plot <- plot_grid(plot_grid(plotlist = plot_list, nrow = 2), legend, ncol = 1, rel_heights = c(10,1))


#### 10/21/2022: Comparing WF MCAR, MAR, MNAR ####
imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

# MCAR
load('results/simulation_noST_MCARp2_R100_10162022.RData')
res_MCAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') %>%
  filter(method %in% c('y_pred_baseline_WF', 'y_pred_CCA_WF'))

res_MCAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MCAR', res_MCAR$method)

# MAR
rm(imputed_list)
load('results/simulation_noST_MARp2_R80_10162022.RData')
res_MAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') %>%
  filter(method == 'y_pred_CCA_WF')

res_MAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MAR', res_MAR$method)

# MNAR
rm(imputed_list)
load('results/simulation_noST_MNARp2_R100_10202022.RData')
res_MNAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') %>%
  filter(method == 'y_pred_CCA_WF')

res_MNAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MNAR', res_MNAR$method)

# combine and plot
res <- rbind(res_MCAR, res_MNAR, res_MAR)

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF_MCAR",  "y_pred_CCA_WF_MAR",  "y_pred_CCA_WF_MNAR")
rename_vec = c('WF full data', 'WF CCA MCAR', 'WF CCA MAR', 'WF CCA MNAR')
color_vec = c('red', 'blue', 'green', 'orange')

plot_metrics_by_point(res = res, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))
