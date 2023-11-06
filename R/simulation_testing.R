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
load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_CAR_freq_results_10212023.RData')
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
files <- grep('2023_10_21',dir('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results', full.names = T), value = T)

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


