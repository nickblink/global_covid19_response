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

#### Results from quasipoisson DGP ####
load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/quasipoisson_results_11202023.RData")

## theta = 2
df = res_theta2$df_miss
# run baseline model (no missing data)
df = WF_baseline(df)
# plot em!
plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), plot_missing_points = F, include_legend = T, PIs = F)
# not much

## theta = 9
df = res_theta9$df_miss
# run baseline model (no missing data)
df = WF_baseline(df)
# plot em!
plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), plot_missing_points = F, include_legend = T, PIs = F)

ggsave('figures/baseline_vs_MNAR_WF_QPtheta9_11192023.png', height = 10, width = 14)

## theta = 100
df = res_theta100$df_miss
# run baseline model (no missing data)
df = WF_baseline(df)
# plot em!
plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), plot_missing_points = F, include_legend = T, PIs = F)

ggsave('figures/baseline_vs_MNAR_WF_QPtheta100_11192023.png', height = 10, width = 14)

## theta = 100, b0 = 6, b1 = 0
df = res_theta100_b6_0$df_miss
df = WF_baseline(df)
# plot em!
plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), plot_missing_points = F, include_legend = T, PIs = F)

ggsave('figures/baseline_vs_MNAR_WF_QPtheta100_b6_0_11192023.png', height = 10, width = 14)

#
#### Results from plots with full data and missing data ####
# So I want the plot to have the full data, the fit of the WF model on the full data, and the fit of the WF model on the missing data.

# pull in data
file_1 = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar03_wf_beta6_n025_2023_10_21"
lst <- combine_results(file_1, return_raw_list = T)
df = lst[[1]]$df_miss

# run baseline model (no missing data)
df = WF_baseline(df)

plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), 
                   plot_missing_points = F, include_legend = T, PIs = F) #fac_list = c('A001','A002','A003','A004'),

ggsave('figures/baseline_vs_MNAR_WF_11192023.png', height = 10, width = 14)

# Now let's look at this with different betas
file_1 = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar03_wf_beta4_0_2023_11_08/"
lst <- combine_results(file_1, return_raw_list = T)
df = lst[[1]]$df_miss

# run baseline model (no missing data)
df = WF_baseline(df)

plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), 
                   plot_missing_points = F, include_legend = T, PIs = F) #fac_list = c('A001','A002','A003','A004'),

ggsave('figures/baseline_vs_MNAR_WF_B4_0_11192023.png', height = 10, width = 14)

# Now with gamma = 2
file_1 = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar03_wf_gamma2_beta6_n025_2023_11_08/"
lst <- combine_results(file_1, return_raw_list = T)
df = lst[[1]]$df_miss

# run baseline model (no missing data)
df = WF_baseline(df)

plot_facility_fits(df, methods = c('y_pred_WF', 'y_pred_baseline_WF'), imp_names = c('MNAR','Full Data'), 
                   plot_missing_points = F, include_legend = T, PIs = F) #fac_list = c('A001','A002','A003','A004'),

ggsave('figures/baseline_vs_MNAR_WF_gamma2_11192023.png', height = 10, width = 14)

#
#### Results and plots of different DGPs and MGPs ####
files <- grep('2023_11_08',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

file_CAR <- grep('car33025', files, value = T)
res_CAR <- combine_results_wrapper(files = file_CAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

file_freqGLM <- grep('freqglm', files, value = T)
res_freq <- combine_results_wrapper(files = file_freqGLM, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_WF4_0 = grep('beta4_0', files, value = T)
res_WF4_0 <- combine_results_wrapper(files = files_WF4_0, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_WF4_n025 = grep('beta4_n025', files, value = T)
res_WF4_n025 <- combine_results_wrapper(files = files_WF4_n025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_WF6_n025 = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar03_wf_beta6_n025_2023_10_21"
res_WF6_n025 <- combine_results_wrapper(files = files_WF6_n025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_WF8_0 = grep('beta8_0', files, value = T)
res_WF8_0 <- combine_results_wrapper(files = files_WF8_0, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_WF8_n025 = grep('beta8_n025', files, value = T)
res_WF8_n025 <- combine_results_wrapper(files = files_WF8_n025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

files_gamma2 = grep('gamma2', files, value = T)
res_gamma2 <- combine_results_wrapper(files = files_gamma2, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))


## plots of the different DGPs stacked.
{
  WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
  p0 <- plot_all_methods(res = res_gamma2$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF DGP and MNAR; gamma = 2', include_legend = F)
  
  p1 <- plot_all_methods(res = res_WF6_n025$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF DGP and MNAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_freq$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'freqGLM DGP and MNAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_CAR$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP and MNAR', include_legend = F)
  
  cowplot::plot_grid(p0$plot, p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,3,1))
  ggsave('figures/MNAR_WFfreqCAR_11172023.png', height = 10, width = 10)
}
## ^ K these are nearly identical to the MCAR results. Weird

## plots of the different beta params stacked
{
  WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
  p1 <- plot_all_methods(res = res_WF4_n025$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'E[B0] = 4; E[B1] = -0.25', include_legend = F)
  
  p2 <- plot_all_methods(res = res_WF4_0$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'E[B0] = 4; E[B1] = 0', include_legend = F)
  
  p3 <- plot_all_methods(res = res_WF6_n025$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'E[B0] = 6; E[B1] = -0.25', include_legend = F)
  
  p4 <- plot_all_methods(res = res_WF8_n025$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'E[B0] = 8; E[B1] = -0.25', include_legend = F)
  
  p5 <- plot_all_methods(res = res_WF8_0$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'E[B0] = 8; E[B1] = 0', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p4$plot, p5$plot, p1$legend, ncol = 1, rel_heights = c(rep(2,5),1))
  ggsave('figures/MNAR_WFbetas_11172023.png', height = 12.5, width = 10)
}

#
#### Plots of missing assumption types ####
# MAR 03
load(dir("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar03_wf_beta6_n025_2023_10_21", full.names = T)[1])

df_MAR = imputed_list[[1]]$df_miss

plot_facility_fits(df_MAR, methods = NULL, plot_missing_points = F,  include_legend = F)

ggsave('figures/MAR03_y_plots_11082023.png', width = 14, height = 10)

# MCAR 03
load(dir("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar03_wf_beta6_n025_2023_10_19", full.names = T)[1])

df_MCAR = imputed_list[[1]]$df_miss

plot_facility_fits(df_MCAR, methods = NULL, plot_missing_points = F,  include_legend = F)

ggsave('figures/MCAR03_y_plots_11082023.png', width = 14, height = 10)

# MNAR 03
load(dir("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar03_wf_beta6_n025_2023_10_21", full.names = T)[1])

df_MNAR = imputed_list[[1]]$df_miss

plot_facility_fits(df_MNAR, methods = NULL, plot_missing_points = T,  include_legend = F)

ggsave('figures/MNAR03_y_plots_11082023.png', width = 14, height = 10)

## All together now!
p1 <- plot_facility_fits(df_MCAR, methods = NULL, plot_missing_points = T,  include_legend = F, fac_list = c('A001','A002','A003','A004'), ncol = 4)

p2 <- plot_facility_fits(df_MAR, methods = NULL, plot_missing_points = T,  include_legend = F, fac_list = c('A001','A002','A003','A004'), ncol = 4)

p3 <- plot_facility_fits(df_MNAR, methods = NULL, plot_missing_points = T,  include_legend = F, fac_list = c('A001','A002','A003','A004'), ncol = 4)

cowplot::plot_grid(p1, p2, p3, ncol = 1, labels = c('MCAR','MAR','MNAR'))

ggsave('figures/districtA_missing_y_plots_11082023.png', width = 14, height = 8)
#
#### Investigating freqGLM - why so bad? ####
# I'm going to just do the WF MCAR data first, 
# Then look at bias
# Then look at the parameter estimates from freqGLM

## Bias analysis
load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_CAR_freq_results_10212023.RData')

p1 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MCAR'), add_lines = list(F, F, F, F),  metrics = c('bias','RMSE'), results_by_point = F, rows = 1, title = 'Data Generated with MCAR Missing Data', include_legend = F)

## Parameter estimates
files_MCAR <- grep('2023_10_19',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

res_MCAR <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), return_unprocessed = T)

# pull the alpha and rho params by simulation
params = NULL
for(p in c('0','01','02','03','04','05')){
  tt = res_MCAR$results[[p]]$freqGLM_lst
  
  param_res = do.call('rbind', lapply(tt, function(xx){
    yy = xx %>% filter(param %in% c('rho', 'alpha'))
    yy
  }))
  
  param_res$p = as.numeric(p)/10
  
  params = rbind(params, param_res)
  
}
params$full_estimate = as.numeric(params$full_estimate)

tmp <- params %>%
  group_by(param, p) %>% 
  summarize(median = median(full_estimate),
            lower = stats::quantile(full_estimate, probs = 0.25),
            upper = stats::quantile(full_estimate, probs = 0.75)) 

ggplot(data = tmp) +
  geom_point(aes(x = p, y = median, color = param), position = position_dodge(width = 0.1)) + 
  geom_errorbar(aes(x = p, ymin = lower, ymax = upper, color = param), position = position_dodge(width = 0.1)) + 
  ggtitle('freqGLM parameter estimates for WF DGP; MCAR')


#
#### Plots of facility fits ####
load(dir("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar02_wf_beta6_n025_2023_10_21", full.names = T)[1])

df = imputed_list[[1]]$df_miss %>%
  filter(facility %in% c('A001','A002'), date <= '2020-01-01') %>%
  mutate(facility = gsub('facility|00','',facility))

# add in the "zoomed in" plots
tmp = df %>% filter(date >= '2019-06-01')
tmp$facility = paste0(tmp$facility, ' zoomed in')
df <- rbind(df, tmp)

plot_facility_fits(df, methods = c('y_pred_WF'), fac_list = c('A1', 'A1 zoomed in','A2', 'A2 zoomed in'), plot_missing_points = F, outbreak_points = c(3,5,10), include_legend = T)

ggsave('figures/sample_WF_fits_10312023.png', height = 5, width = 7)

#
#### Metrics plots! ####
load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_CAR_freq_results_10212023.RData')
#load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_CAR_freq_results_10122023.RData')

res_facility$method = gsub('CARstan', 'CARBayes',res_facility$method)
res_district$method = gsub('CARstan', 'CARBayes',res_district$method)

table(res_district$sim, res_district$prop_missing)
table(res_facility$sim, res_facility$prop_missing)
# CAR_33025   CAR_331   freqGLM    WF_MAR   WF_MCAR   WF_MNAR

res_facility$method = factor(res_facility$method, levels = c('WF','freqGLM','CARBayes'))
res_district$method = factor(res_district$method, levels = c('WF','freqGLM','CARBayes'))

## facility plots of the 3 WF methods stacked
{
  WF_ylims = list(ylim(0.5,1), ylim(0.5,1), ylim(0.5, 1), ylim(0.5,1))
  p1 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MCAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MCAR Missing Data', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MAR Missing Data', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MNAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MNAR Missing Data', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_facility_plots_10302023.png', height = 7.5, width = 10)
}

## Single WF MCAR plot
{
  p1 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MCAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF Model and MCAR Missing Data', include_legend = F)
  cowplot::plot_grid(p1$plot, p1$legend, ncol = 1, rel_heights = c(3,1))
  
  ggsave('figures/WF_MCAR_plot_11012023.png',height = 3, width = 10)
}

## district plots of the 3 WF methods stacked
{
  WF_ylims = list(ylim(0.25,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1))
  p1 <- plot_all_methods(res = res_district %>% filter(sim == 'WF_MCAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MCAR Missing Data', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district %>% filter(sim == 'WF_MAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MAR Missing Data', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district %>% filter(sim == 'WF_MNAR'), fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with MNAR Missing Data', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_district_plots_10302023.png', height = 7.5, width = 10)
}

## facility plots of WF, CAR, freqGLM stacked
{
  DGP_ylims = list(ylim(0.25,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1))
  p1 <- plot_all_methods(res = res_facility %>% filter(sim == 'WF_MCAR'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility %>% filter(sim == 'freqGLM'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility %>% filter(sim == 'CAR_33025'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CARBayes model', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_10302023.png', height = 7.5, width = 10)
}

## district plots of WF, CAR, freqGLM stacked
{
  DGP_ylims = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1))
  p1 <- plot_all_methods(res = res_district %>% filter(sim == 'WF_MCAR'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district %>% filter(sim == 'freqGLM'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district %>% filter(sim == 'CAR_33025'), fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CARBayes model', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_district_comparison_plots_10302023.png', height = 7.5, width = 10)
}


#
#### Combining all results (WF w/ MCAR, MAR, MNAR; freqGLM MCAR; CAR_331 MCAR, CAR_33025 MCAR) - changing 10/21/2023 ####
files <- grep('2023_10_21',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

## Pull in WF MCAR results
{
  files_MCAR <- grep('2023_10_19',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)
  
  res_MCAR <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))
  
  res_MCAR$results$sim = 'WF_MCAR' 
  res_MCAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 
  
  res_MCAR_district <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)
  
  res_MCAR_district$results$sim = 'WF_MCAR' 
  res_MCAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 
}

## Pull in MNAR results
{
  files_MNAR = grep('mnar', files, value = T)
  
  res_MNAR <- combine_results_wrapper(files = files_MNAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))
  
  res_MNAR$results$sim = 'WF_MNAR'
  res_MNAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MNAR' 
  
  res_MNAR_district <- combine_results_wrapper(files = files_MNAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)
  
  res_MNAR_district$results$sim = 'WF_MNAR'
  res_MNAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MNAR' 
}

## Pull in MAR results
{
  files_MAR = grep('mar', files, value = T)
  
  res_MAR <- combine_results_wrapper(files = files_MAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))
  
  res_MAR$results$sim = 'WF_MAR' 
  res_MAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MAR' 
  
  res_MAR_district <- combine_results_wrapper(files = files_MAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)
  
  res_MAR_district$results$sim = 'WF_MAR' 
  res_MAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MAR' 
}
# had to remove a few results, but not enough to warrant concern (between 1-2/1000)

## Pull in CAR results
{
  files_33025 = grep('car33025', files, value = T)
  
  # res_CAR331 <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))
  # 
  # res_CAR331$results$sim = 'CAR_331'
  # res_CAR331$results$sim_long = 'CAR_DGP(0.3,0.3,1,b0=6,b1=-0.25);MCAR' 
  # res_CAR331_district <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'), district_results = T)
  # 
  # res_CAR331_district$results$sim = 'CAR_331'
  # res_CAR331_district$results$sim_long = 'CAR_DGP(0.3,0.3,1,b0=6,b1=-0.25);MCAR' 
  # 
  res_CAR33025 <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM','CARstan'))
  
  res_CAR33025$results$sim = 'CAR_33025'
  res_CAR33025$results$sim_long = 'CAR_DGP(0.3,0.3,0.25,b0=6,b1=-0.25);MCAR' 
  res_CAR33025_district <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM','CARstan'), district_results = T)
  
  res_CAR33025_district$results$sim = 'CAR_33025'
  res_CAR33025_district$results$sim_long = 'CAR_DGP(0.3,0.3,0.25,b0=6,b1=-0.25);MCAR' 
}

## Pull in freqGLM results 
{
  files_freqGLM <- grep('freqglm', files, value = T)
  
  res_freqGLM <- combine_results_wrapper(files = files_freqGLM, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))
  
  res_freqGLM$results$sim = 'freqGLM' 
  res_freqGLM$results$sim_long = 'freqGLM_DGP(b0=5,b1=-0.25);MCAR' 
  
  res_freqGLM_district <- combine_results_wrapper(files = files_freqGLM, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)
  
  res_freqGLM_district$results$sim = 'freqGLM' 
  res_freqGLM_district$results$sim_long = 'WF_DGP(b0=5,b1=-0.25);MCAR' 
}

## Join these together
res_facility = rbind(res_MCAR$results, res_MNAR$results, res_CAR33025$results, res_freqGLM$results, res_MAR$results)
res_district = rbind(res_MCAR_district$results, res_MNAR_district$results, res_CAR33025_district$results, res_freqGLM_district$results, res_MAR_district$results)

## Copy WF MCAR p = 0 to MAR p = 0 and MNAR p = 0 (since these are all the same)
tmp = res_facility %>% filter(sim == 'WF_MCAR', prop_missing == 0)
tmp2 = tmp; tmp2$sim = 'WF_MAR'
tmp3 = tmp; tmp3$sim = 'WF_MNAR'
res_facility = rbind(res_facility, tmp2, tmp3)

tmp = res_district %>% filter(sim == 'WF_MCAR', prop_missing == 0)
tmp2 = tmp; tmp2$sim = 'WF_MAR'
tmp3 = tmp; tmp3$sim = 'WF_MNAR'
res_district = rbind(res_district, tmp2, tmp3)

# save(res_facility, res_district, file = 'C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_CAR_freq_results_10212023.RData')

#
#### Combining all results part 2 (WF w/ MCAR, ...) ####

files <- grep('2023_10_19',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

## Pull in MCAR results
{
  files_MCAR = grep('mcar', files, value = T)
  
  res_MCAR <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan', 'y_CARBayesST'), rename_vec = c('WF','freqGLM', 'CARstan', 'CARBayes'))
  
  res_MCAR$results$sim = 'WF_MCAR' 
  res_MCAR$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 
  
  res_MCAR_district <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), district_results = T)
  
  res_MCAR_district$results$sim = 'WF_MCAR' 
  res_MCAR_district$results$sim_long = 'WF_DGP(b0=6,b1=-0.25);MCAR' 
}


plot_all_methods(res = res_MCAR$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF DGP: MCAR; facility-level')

#
#### Analyzing facility naming debug results ####
files <- grep('2023_10_18',dir('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

## WF
files_WF = grep('wf', files, value = T)

res_WF <- combine_results_wrapper(files = files_WF, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

plot_all_methods(res = res_WF$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F), results_by_point = F, rows = 1, title = 'WF DGP: MCAR; facility-level')

## CAR
files_CAR = grep('car_', files, value = T)

res_CAR <- combine_results_wrapper(files = files_CAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

plot_all_methods(res = res_CAR$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F), results_by_point = F, rows = 1, title = 'CAR DGP: MCAR; facility-level')

## freqGLMepi
files_freq = grep('glm_', files, value = T)

res_freq <- combine_results_wrapper(files = files_freq, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

plot_all_methods(res = res_freq$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F), results_by_point = F, rows = 1, title = 'freqGLM DGP: MCAR; facility-level')

#
#### Analyzing results by point ####
files <- grep('2023_10_11',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

## Pull in MCAR results
files_MCAR = grep('mcar', files, value = T)

res_MCAR <- combine_results_wrapper(files = files_MCAR, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), results_by_point = T)
  
tt = res_MCAR$results

tmp = tt %>% filter(prop_missing == 0, method == 'CARstan')
View(tmp)

tmp2 = tt %>% filter(prop_missing == 0.5, method == 'CARstan')
View(tmp2)

p1 <- plot_all_methods(res = res_MCAR$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF DGP: MCAR; facility-level')


files <- grep('2023_10_12',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)
files_freqGLM <- grep('freqglm', files, value = T)

res_freqGLM <- combine_results_wrapper(files = files_freqGLM, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'), results_by_point = T)

tt = res_freqGLM$results

tmp3 = tt %>% filter(prop_missing == 0.5, method == 'CARstan')
View(tmp3)

## Time to look at some betas
dir(files_MCAR[1], full.names = T)[1] -> tmp
load(tmp)
imputed_list[[1]]$WF_betas
imputed_list[[2]]$WF_betas


#
#### Processing CAR results ####
files <- grep('2023_10_06',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

files_331 = grep('car331', files, value = T)
files_33025 = grep('car33025', files, value = T)

res1 <- combine_results_wrapper(files = files_331, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res2 <- combine_results_wrapper(files = files_33025, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

table(res1$params$tau2_DGP)
table(res2$params$tau2_DGP)
# Good

plot_all_methods(res = res1$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=1): MCAR')

plot_all_methods(res = res2$results, fix_axis = list(ylim(0,1), ylim(0.25,1), ylim(0.5, 1), ylim(0.5,1)),  add_lines = list(0.95, F, F, F), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=0.25): MCAR')


res1_district <- combine_results_wrapper(files = files_331, district_results = T, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

res2_district <- combine_results_wrapper(files = files_33025, district_results = T, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

plot_all_methods(res = res1_district$results, fix_axis = list(ylim(0,1), ylim(0,1), ylim(0.3, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=1): MCAR')

plot_all_methods(res = res2_district$results, fix_axis = list(ylim(0,1), ylim(0,1), ylim(0.3, 1), ylim(0.5,1)), add_lines = list(0.95, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'CAR DGP (rho=0.3; alpha=0.3; tau2=.025): MCAR')

#save(res1, res2, res1_district, res2_district, file = 'C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/results_CARDGP_10102023.RData')

     
#
#### (one time run) Fix file naming ####
files <- grep('2023_10_06',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

res <- combine_results_wrapper(files = files, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARBayesST','y_CARstan'), rename_vec = c('WF','freqGLM','CARBayes','CARstan'))

params = res$params

ind <- grep('\\(1\\)', params$file)
table(params$tau2_DGP[ind])
table(params$tau2_DGP[-ind])

# files_to_move = params$file[ind]
# for(f in files_to_move){
#   new_file = gsub('car331','car33025', f)
#   new_file = gsub('\\(1\\)','', new_file)
#   file.rename(f, new_file)
# }

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


