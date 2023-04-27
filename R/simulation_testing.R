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


#### 4/21/2023: Comparing size of training data on results ####
load('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/WF_n_comparisons_04162023.RData')

imp_vec = c("y_pred_CCA_WF")
lst4 <- lapply(imputed_list4, '[[', 1)
res4 <- calculate_metrics_by_point(lst4, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 

lst8 <- lapply(imputed_list8, '[[', 1)
res8 <- calculate_metrics_by_point(lst8, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 

lst12 <- lapply(imputed_list12, '[[', 1)
res12 <- calculate_metrics_by_point(lst12, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') 

res4$method = 'WF_4yr'
res8$method = 'WF_8yr'
res12$method = 'WF_12yr'

res <- rbind(res4, res8, res12)
res$method = factor(res$method, levels = c('WF_4yr', 'WF_8yr', 'WF_12yr'))
res$prop_missing = .01

plot_all_methods(res = res, fix_axis = F)

###

beta4 <- lapply(imputed_list4, '[[', 3)
beta8 <- lapply(imputed_list8, '[[', 3)
beta12 <- lapply(imputed_list12, '[[', 3)

plot_beta_bias <- function(beta_estimates, true_betas, title = 'coefficient bias across facilities'){
  est_betas <- Reduce("+", beta_estimates)/length(beta_estimates)
  
  bias <- est_betas - true_betas
  colnames(bias)[1] <- 'intercept'
  
  res2 <- NULL
  for(col in colnames(bias)){
    
    x <- bias[,col]
    
    res2 <- rbind(res2,
                  data.frame(beta = col,
                             bias = mean(x),
                             lower = stats::quantile(x, probs = 0.25),
                             upper = stats::quantile(x, probs = 0.75)))
  }
  
  res2$beta <- factor(res2$beta, levels = res2$beta)
  p1 <- ggplot() + 
    geom_point(data = res2, aes(x = beta, y = bias), position = position_dodge(width = 0.1)) +
    geom_errorbar(data = res2, aes(x = beta, y = bias, ymin = lower, ymax = upper), position = position_dodge(width = 0.1)) +
    geom_hline(yintercept = 0) + 
    ylab('coefficient bias') + 
    ylim(-.01, .002) + 
    guides(alpha = 'none') +
    ggtitle(title) +
    theme_bw()
  
  return(p1)
}

plot_beta_bias(beta4, true_betas4, title = 'bias from 4 years training data')

plot_beta_bias(beta8, true_betas8, title = 'bias from 8 years training data')

plot_beta_bias(beta12, true_betas12, title = 'bias from 12 years training data')


  
#
#### 3/25/2023: Analyzing coefficients' results ####
file <- "C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/mcar01_betatest_nost_beta43_n025_2023_03_16"

lst <- combine_results(input_folder = file, return_lst = T, results_file = NULL)

est_betas <- Reduce("+", lst$WF_lst)/length(lst$WF_lst)

bias <- est_betas - lst$true_betas
colnames(bias)[1] <- 'intercept'

res2 <- NULL
for(col in colnames(bias)){
  
  x <- bias[,col]
  
  res2 <- rbind(res2,
                data.frame(beta = col,
                           bias = mean(x),
                           lower = stats::quantile(x, probs = 0.25),
                           upper = stats::quantile(x, probs = 0.75)))
}

res2$beta <- factor(res2$beta, levels = res2$beta)
ggplot() + 
  geom_point(data = res2, aes(x = beta, y = bias), position = position_dodge(width = 0.1)) +
  geom_errorbar(data = res2, aes(x = beta, y = bias, ymin = lower, ymax = upper), position = position_dodge(width = 0.1)) +
  geom_hline(yintercept = 0) + 
  ylab('coefficient bias') + 
  guides(alpha = 'none') +
  ggtitle('coefficient bias across facilities') +
  theme_bw()


plot(lst$true_betas[,1], bias[,1], xlab = 'true beta', ylab = 'bias')

#
#### 3/25/2023: Plotting bias by date ####
file <- "C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/mcar02_nost_beta6_n025_2023_03_16" 

res <- grab_results(file, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF_MCAR','CAR_MCAR','freqGLM_MCAR')) %>%
  filter(method == 'WF_MCAR')

tmp <- res %>%
  group_by(method, prop_missing, date) %>% 
  summarize(median = median(bias),
            lower = stats::quantile(bias, probs = 0.25),
            upper = stats::quantile(bias, probs = 0.75)) 

ggplot() + 
  geom_point(data = tmp, aes(x = date, y = median, color = method), position = position_dodge(width = 0.1)) +
  geom_errorbar(data = tmp, aes(x = date, y = median, ymin = lower, ymax = upper, color = method), position = position_dodge(width = 0.1)) +
  geom_hline(yintercept = 0) + 
  ylab('bias') + 
  guides(alpha = 'none') +
  ggtitle('bias from all facilities') 
  theme_bw()

  tmp <- res %>%
    group_by(method, prop_missing, date, facility) %>% 
    summarize(median = median(bias),
              lower = stats::quantile(bias, probs = 0.25),
              upper = stats::quantile(bias, probs = 0.75)) 
  
  ggplot() + 
    facet_wrap(vars(facility)) +
    geom_point(data = tmp, aes(x = date, y = median, color = method), position = position_dodge(width = 0.1)) +
    #geom_errorbar(data = tmp, aes(x = date, y = median, ymin = lower, ymax = upper, color = method), position = position_dodge(width = 0.1)) +
    geom_hline(yintercept = 0) + 
    ylab('bias') + 
    guides(alpha = 'none') +
    ggtitle('bias from each facilities') +
  theme_bw()
  
#
#### 3/24/2023: Getting, grouping, and plotting different beta runs together ####
file_MCAR <- grep('2023_03_16', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MCAR <- grab_results(file_MCAR, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF_MCAR','CAR_MCAR','freqGLM_MCAR')) %>%
  filter(method == 'WF_MCAR')

res_MCAR$method = paste0('WF_b0_', res_MCAR$b0_mean, '_b1_', res_MCAR$b1_mean)

res_1 = res_MCAR %>% filter(method %in% c('WF_b0_4.3_b1_-0.25',
                        'WF_b0_6/4.3_b1_-0.25',
                        'WF_b0_6_b1_-0.25'))

res_2 = res_MCAR %>% filter(method %in% c('WF_b0_4.3_b1_-0.25',
                                          'WF_b0_6_b1_-0.25',
                                          "WF_b0_6_b1_0",
                                          "WF_b0_4.3_b1_0"))
tt <- plot_all_methods(res = res_1, 
                       fix_axis = F) 
                       #metrics = c('bias', 'relative_bias', 'RMSE', 'coverage95', 'interval_width','outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))

tt2 <- plot_all_methods(res = res_2, 
                       fix_axis = F) 

#ggsave(plot = tt, filename = 'figures/WF_beta_comparison_03242023.png', width = 7, height = 4)

#ggsave(plot = tt2, filename = 'figures/WF_beta_comparison2_03242023.png', width = 7, height = 4)

#### 3/3/2023: Getting, grouping, and plotting MCAR, MAR, MNAR together ####

file_MCAR <- grep('mcar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

res_MCAR <- grab_results(file_MCAR[c(1:6,11)], imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF_MCAR','CAR_MCAR','freqGLM_MCAR'))

file_MAR <- grep('mar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MAR <- grep('02_26', file_MAR, value = T)

res_MAR <- grab_results(file_MAR, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF_MAR','CAR_MAR','freqGLM_MAR'))

file_MNAR <- grep('mnar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MNAR <- file_MNAR[c(2,4,6,8,10,12,15)]

res_MNAR <- grab_results(file_MNAR, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF_MNAR','CAR_MNAR','freqGLM_MNAR'))

res_all <- rbind(res_MCAR, res_MAR, res_MNAR)

# save(res_all, file = 'C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/res_combined_03032023.RData')

res_WF <- res_all %>%
  filter(method %in% c('WF_MCAR', 'WF_MAR', 'WF_MNAR'))

tt <- plot_all_methods(res = res_WF, fix_axis = F, metrics = c('bias', 'relative_bias', 'RMSE', 'coverage95', 'interval_width','outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))

ggsave(plot = tt, filename = 'figures/WF_comparison_03062023.png', width = 7, height = 4)
#
#### 2/27/2023: Getting (past) MCAR results ####
file_MCAR <- grep('mcar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MCAR <- file_MCAR[c(1:6,11)]
# file_MCAR <- file_MCAR[1]

tt <- plot_all_methods(file_MCAR, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF','CAR','freqGLM'))
# ggsave(plot = tt, filename = 'figures/MCAR_02272023.png', width = 7, height = 4)


#
#### 2/27/2023: Getting MAR results ####
file_MAR <- grep('mar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MAR <- grep('02_26', file_MAR, value = T)

tt <- plot_all_methods(file_MAR, imp_vec = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF','CAR','freqGLM'))
#ggsave(plot = tt, filename = 'figures/MAR_02272023.png', width = 7, height = 4)

#
#### 2/24/2023: Getting MNAR results ####
file_MNAR <- grep('mnar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file_MNAR <- file_MNAR[c(2,4,6,8,10,12,15)]

tt <- plot_all_methods(file_MNAR)
#ggsave(plot = tt, filename = 'figures/MNAR_02272023.png', width = 7, height = 4)

#
#### 2/21/2023: Checking error handling ####
file_MCAR <- grep('mcar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

file <- 'C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/mcar06_nost_beta43_n025_2023_02_20'

lst_full <- combine_results(input_folder = file, return_lst = T, results_file = NULL)

file <- "C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/mnar06_nost_beta43_n025_2023_02_20" 

lst_full <- combine_results(input_folder = file, return_lst = T, results_file = NULL)

nas <- c()
for(i in 1:320){
  tt <- lst_full$df_lst[[i]]
  print(ncol(tt))
  nas <- c(nas, sum(is.na(tt)))
}

#
#### 01/19/2023: Plotting across methods ####
file_MCAR <- grep('mcar', dir('C:/Users/Admin-Dell/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results', full.names = T), value = T)

p1 <- plot_all_methods(file_MCAR[1:7])

ggsave(plot = p1, filename = 'C:/Users/Admin-Dell/Documents/github_projects/global_covid19_response/figures/MCAR_metrics_01182023.png', height = 4, width = 6)

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

#### 12/14/2022: Testing different beta generation ####
df1 <- simulate_data(district_sizes = c(4, 6, 10), R = 2, end_date = '2020-12-01', b0_mean = 4.3, b1_mean = -0.25)$df_list[[1]]

df2 <- simulate_data(district_sizes = c(4, 6, 10), R = 2, end_date = '2020-12-01', b0_mean = 6, b1_mean = -0.25)$df_list[[1]]

df3 <- simulate_data(district_sizes = c(4, 6, 10), R = 2, end_date = '2020-12-01', b0_mean = 8.05, b1_mean = -1)$df_list[[1]]

df4 <- simulate_data(district_sizes = c(4, 6, 10), R = 2, end_date = '2020-12-01', b0_mean = 3.05, b1_mean = 0)$df_list[[1]]

mean(df2$y_exp)/mean(df1$y_exp)

mean(df2$y_exp)

mean(df1 %>% filter(date == '2020-01-01') %>% pull(y_exp))

mean(df2 %>% filter(date == '2020-01-01') %>% pull(y_exp))

mean(df3 %>% filter(date == '2020-01-01') %>% pull(y_exp))

mean(df4 %>% filter(date == '2020-01-01') %>% pull(y_exp))

# 1, 3, and 4 have the same mean at 1/1/2020. 2 has a higher mean. So how will the outbreak detections compare? That is a good question, my friend. To test! I will compere, my friends

#### combining results ####
dropbox_results <- 'C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results'

res <- combine_results(input_folder = sprintf('%s/mcar07_nost_2022_11_17', dropbox_results), return_lst = T, results_file = NULL)

for(f in grep('mcar',dir(dropbox_results), value = T)){
  res <- combine_results(input_folder = sprintf('%s/%s', dropbox_results,f), return_lst = T, results_file = NULL)
}

for(f in grep('mnar',dir(dropbox_results), value = T)){
  print(f)
  res <- combine_results(input_folder = sprintf('%s/%s', dropbox_results,f), return_lst = T, results_file = NULL)
}
# that is unfortunate

for(f in grep('mar',dir(dropbox_results), value = T)){
  print(f)
  res <- combine_results(input_folder = sprintf('%s/%s', dropbox_results,f), return_lst = T, results_file = NULL)
}

object_sizes <- c()
for(i in 1:1000){
  object_sizes <- c(object_sizes, object.size(res[[i]]))
}

table(object_sizes)
# For MCAR_05, this worked fine for 16 of them, but crashed for the rest, why? It failed for 3, 6, 7, 11
# For MCAR_05, it was a timing error for 11. Not for the others though. Weird 
# For MCAR_07 and 08, it was an error in R. Something "infinite or missing values in 'x'". Idk this is for me to figure out another time.

# For MCAR_07 and MCAR_08, no results

#### 10/21/2022: Comparing WF MCAR, MAR, MNAR ####
res_folder <- 'C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results'
imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

# MCAR
load(sprintf('%s/simulation_noST_MCARp2_R100_10162022.RData', res_folder))
res_MCAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F, min_date = '2020-01-01') %>%
  filter(method %in% c('y_pred_baseline_WF', 'y_pred_CCA_WF'))

res_MCAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MCAR', res_MCAR$method)

# MAR
rm(imputed_list)
load(sprintf('%s/simulation_noST_MARp2_R80_10162022.RData', res_folder))
res_MAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_CCA_WF')

res_MAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MAR', res_MAR$method)

# MNAR
rm(imputed_list)
load(sprintf('%s/simulation_noST_MNARp2_R100_10202022.RData', res_folder))
res_MNAR <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_CCA_WF')

res_MNAR$method <- gsub('y_pred_CCA_WF', 'y_pred_CCA_WF_MNAR', res_MNAR$method)

# combine and plot
res <- rbind(res_MCAR, res_MNAR, res_MAR)

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF_MCAR",  "y_pred_CCA_WF_MAR",  "y_pred_CCA_WF_MNAR")
rename_vec = c('WF full data', 'WF CCA MCAR', 'WF CCA MAR', 'WF CCA MNAR')
color_vec = c('red', 'blue', 'green', 'orange')

plot_metrics_by_point(res = res, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

#### 10/20/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MNAR ####

R = 1

#results_file <- 'results/simulation_noST_MNARp2_R100_10202022.RData'
results_file <- 'C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/simulation_noST_MNARp2_R100_10202022.RData'

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MNAR_sim(df, p = 0, direction = 'upper', gamma = 1, by_facility = T)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
    
    if(i %% 5 == 0){
      print(sprintf('saving results for i = %s',  i))
      save(imputed_list, i, file = results_file)
    }
  }
}) # ~15hrs

# save(imputed_list, file = results_file)

# over-writing expected with the simulated value because I wasn't storing these before. And now I am!
imputed_list <- lapply(1:100, function(ii){
  tmp <- imputed_list[[ii]]
  tmp$y_exp <- tmp$y_true
  tmp
})


imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width','outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  filter(method == 'y_pred_baseline_WF') %>%
  group_by(date) %>%
  summarize(num_missing = sum(num_missing))

# number missing points across simulations
plot(res$date, res$num_missing, ylab = 'number of times missing', xlab = 'date', main = 'missingness from R = 100 with 20 facilities')

for(i in 1:10){
  test <- imputed_list[[i]]  %>% filter(is.na(y))
  print(max(test$date))
}

#### 10/16/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MCAR ####

R = 100

results_file <- 'results/simulation_noST_MCARp2_R200_10162022.RData'

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MCAR_sim(df, p = 0.2, by_facility = T)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF baseline (complete data) imputation
    df_miss <- WF_baseline(df_miss, R_PI = 100)
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
    
    if(i %% 5 == 0){
      print(sprintf('saving results for i = %s',  i))
      save(imputed_list, i, file = results_file)
    }
  }
}) # 15 hrs R = 100

# save(imputed_list, file = results_file)

load('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/results/simulation_noST_MCARp2_R100_10162022.RData')

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, min_date = '2020-01-01', imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  arrange(coverage95) %>%
  filter(method == 'y_pred_baseline_WF')


#### 10/16/2022: WF Baseline, WF CCA, Epi CCAm and CAR CCA with MAR ####

R = 80

system.time({
  lst <- simulate_data(district_sizes = c(4, 6, 10), R = R, end_date = '2020-12-01')
  
  imputed_list = list()
  res_full = res_imputed = NULL
  for(i in 1:R){
    df <- lst$df_list[[i]]
    
    # simulation function!
    df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.7, alpha = 0.7, tau = 3)
    
    # run the freqGLM_epi complete case analysis
    freqGLMepi_list = freqGLMepi_CCA(df_miss, R_PI = 100, verbose = F)
    df_miss <- freqGLMepi_list$df
    
    # run the WF baseline (complete data) imputation
    df_miss <- WF_baseline(df_miss, R_PI = 100)
    
    # run the WF complete case analysis model
    df_miss <- WF_CCA(df_miss, col = "y", family = 'poisson', R_PI = 100)
    
    # run the CAR complete case analysis model
    df_miss <- CARBayes_CCA(df_miss, burnin = 5000, n.sample = 10000, prediction_sample = T, model = 'facility_fixed', predict_start_date = '2016-01-01')
    
    imputed_list[[i]] = df_miss
  }
}) # Started at 8pm on 10/16. Should be 12 hours, so ready by the morn

# save(imputed_list, file = 'results/simulation_noST_MARp2_R200_10162022.RData')

imp_vec = c("y_pred_baseline_WF", "y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi")
rename_vec = c('WF full data', 'WF CCA', 'CAR CCA', 'freqEpi CCA')
color_vec = c('red', 'blue', 'green', 'orange')

plot_facility_fits(imputed_list[[2]], imp_vec = imp_vec, imp_names = rename_vec, color_vec = color_vec)

plot_metrics_by_point(imputed_list, imp_vec = imp_vec, color_vec = color_vec, min_date = '2020-01-01', 
                      min_missing = 0, rename_vec = rename_vec, metric_list = c('bias','RMSE','coverage95','interval_width'))

res <- calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, min_date = '2020-01-01', imputed_only = F, rm_ARna = F, use_point_est = F) %>%
  arrange(coverage95) %>%
  filter(method == 'y_pred_baseline_WF')


#### Analyzing the missingness pattern ####
plot_missingness <- function(df_miss){
  df_spread = df_miss %>%
    select(date, facility,y) %>%
    tidyr::spread(., facility, y)
  
  tmp2 = df_spread[,-c(1)]
  for(col in colnames(tmp2)){
    tmp2[,col] = as.integer(is.na(tmp2[,col]))
  }
  
  tmp2 = as.matrix(tmp2)
  rownames(tmp2) = as.character(df_spread$date)
  heatmap(tmp2, keep.dendro = F, Rowv = NA, )
  
  gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = F, Colv = F, xlab = 'facilities', trace = 'none', key = F)
}

lst <- simulate_data(district_sizes = c(4, 6, 10), R = 1, end_date = '2019-12-01')

df <- lst$df_list[[1]] 

# simulation function!
df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0, alpha = 0, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.9, alpha = 0, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0, alpha = 0.9, tau = 3)

plot_missingness(df_miss)

df_miss <- MAR_spatiotemporal_sim(df, p = 0.2, rho = 0.9, alpha = 0.9, tau = 3)

plot_missingness(df_miss)


