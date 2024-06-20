current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')
source('R/imputation_functions.R')
library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)

if(file.exists('C:/Users/Admin-Dell')){
  res_dir = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results"
}else{
  res_dir = "C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results"
}

# get the results and name the simulations
get_results <- function(file_names, sim_name, district_results = F, ...){
  res <- combine_results_wrapper(files = file_names, methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARBayes'), district_results = district_results, ...)
  
  res$results$sim = sim_name 
  
  return(res)
}

# aggregate the results across simulations
aggregate_results <- function(res, bar_quants = c(0.25, 0.75), metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), full_results = T){
  # rename the metrics
  if(!is.null(metric_rename)){
    for(i in 1:length(metrics)){
      colnames(res)[which(colnames(res) == metrics[i])] <- metric_rename[i]
    }
    metrics <- metric_rename
  }
  
  res_lst <- list()
  
  for(metric in metrics){
    tmp <- res %>%
      group_by(method, prop_missing) %>% 
      summarize(median = median(get(metric)),
                lower = stats::quantile(get(metric), probs = bar_quants[1]),
                upper = stats::quantile(get(metric), probs = bar_quants[2])) 
    if(full_results){
      res_lst[[metric]] <- tmp
    }else{
      res_lst[[metric]] <- tmp %>% group_by(method) %>%
        summarize(min = min(median), 
                  max = max(median))
    }
  }
  
  return(res_lst)
}

# get CAR convergence values from the results in a simulation.
get_CARconvergence <- function(res){
  # get CAR summary vals.
  CARsums <- res$results[[1]]$CARstan_summary
  
  # get ESS values.
  ESS <- sapply(CARsums, function(xx) xx[,'n_eff']) %>%
    apply(1, median)
  
  # make ESS more concise.
  ESS <- c(ESS[1:3], betas = median(ESS[4:length(ESS)]))
  
  # get rhat values.
  rhat <- sapply(CARsums, function(xx) xx[,'Rhat']) %>%
    apply(1, median)
  
  # make rhat concise.
  rhat <- c(rhat[1:3], betas = median(rhat[4:length(rhat)]))
  
  # print the results
  print('Median effective sample size'); print(ESS)
  print('Median rhat value'); print(rhat)
  #print(sprintf('Median effective sample sizes are %s. \n Median rhat values are %s', ESS, rhat))
  
  return(list(ESS = ESS, rhat = rhat))
}

# load the full main paper results
load(paste0(res_dir,'/full_paper_results_12262023.RData'))

#### Comparing new and old freqGLM - why so different? ####
files = grep('2023_12_14',dir(res_dir, full.names = T), value = T)

# freqGLM0202 EB0 = 5.5, EB1 = -0.25: MCAR (not actually the appendix though)
files_freqglm0202_MCAR_beta55_n025 <- grep('beta055_beta1n025',grep('freqglm', files, value = T), value = T)

load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_freqglm0202_beta055_beta1n025_id618184_2023_12_14/sim_results_p0.0_mcar_1(50).RData')
params_OG <- params; rm(params)
arguments_OG <- arguments; rm(arguments)
imputed_list_OG <- imputed_list; rm(imputed_list)
true_betas_OG <- true_betas;rm(true_betas)

load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_freqglm0202_beta055_beta1n025_id447944_2024_06_17/sim_results_p0.0_mcar_1(50).RData')
params_new <- params; rm(params)
arguments_new <- arguments; rm(arguments)
imputed_list_new <- imputed_list; rm(imputed_list)
true_betas_new <- true_betas;rm(true_betas)

names(params_OG)
names(params_new)
params_OG
params_new

#
#### Making new DGP plots ####
load(paste0(res_dir, '/DGP_comparison_results_higherCARnsample_06192024.RData'))

DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_beta6_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_beta6_beta10"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model (BETA1 = 0), assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_nsample10000_burnin5000_06192024.pdf', height = 7.5, width = 10)
}

#
#### Getting new DGP results with higher n.sample and burnin ####
files = grep('2024_06_17',dir(res_dir, full.names = T), value = T)
files = files[-grep('91167', files)]

for(d in files){
  ff = dir(d)
  if(length(ff) < 51){
    print(d)
    print(length(ff))
    print('-------')
    #unlink(d, recursive = T) # deletes the directory
  }
}

files_WF_MCAR_beta6_n025 <- grep('wf_beta06', files, value = T)

files_freqglm0202_MCAR_beta55_n025 <- grep('beta055_beta1n025',grep('freqglm', files, value = T), value = T)

files_CAR33025_MCAR_beta6_beta10 <- grep('car0303025_beta06_beta10', files, value = T)


file_name_str <- c('files_WF_MCAR_beta6_n025',
                   'files_freqglm0202_MCAR_beta55_n025',
                   'files_CAR33025_MCAR_beta6_beta10')

all_file_names <- unlist(sapply(file_name_str, get))
setdiff(files, all_file_names)
setdiff(all_file_names, files)
tt = table(all_file_names); tt[tt>1]
# ok all gravy

# initialize
res_facility = list()
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  res_facility[[name]] <- get_results(get(name), name, expected_sims = NULL)
  res_district[[name]] <-  get_results(get(name), name,  district_results = T, expected_sims = NULL)
}

# check the burnin and nsample
load(dir(files_CAR33025_MCAR_beta6_beta10[1], full.names = T)[1])
params$CARburnin
params$CARnsample

# check the CAR convergence!
tmp1 <- get_results(files_WF_MCAR_beta6_n025, 'files_WF_MCAR_beta6_n025', return_unprocessed = T) %>% get_CARconvergence()
tmp2 <- get_results(files_freqglm0202_MCAR_beta55_n025, 'files_freqglm0202_MCAR_beta55_n025', return_unprocessed = T, expected_sims = NULL) %>% get_CARconvergence()

# save the results
save(res_facility, res_district, file = paste0(res_dir, '/DGP_comparison_results_higherCARnsample_06192024.RData'))


#
#### Get the results by point for each district ####
# getting the file names
{
  files = grep('2023_12_09|2023_12_10',dir(res_dir, full.names = T), value = T)
  files2 = grep('2023_12_14',dir(res_dir, full.names = T), value = T)
  
  files_MCAR <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)
  
  files_CAR <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)
  
  files <- unique(c(files, files_MCAR, files_CAR))
  
  for(d in files){
    ff = dir(d)
    if(length(ff) < 51){
      print(d)
      print(length(ff))
      print('-------')
      #unlink(d, recursive = T) # deletes the directory
    }
  }
  # all good
  
  # (1) WF MCAR; poisson; beta = 6, -0.25
  files_WF_MCAR_6n025 <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)
  
  # (2) freqGLM0202 EB0 = 5.5, EB1 = -0.25: MCAR
  files_freqglm0202_MCAR_55n025 <- grep('beta055_beta1n025',grep('freqglm', files2, value = T), value = T)
  
  # (3) CAR0303025 MCAR; poisson; beta = 6, -0.25
  files_CAR33025_MCAR_6n025 <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)
  
  # (4) WF MCAR; QP theta = 4; beta = 6,-0.25
  files_WF_MCAR_QP4_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # (5) WF MAR; QP theta = 4; beta = 6,-0.25
  files_WF_MAR_QP4_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # (6) WF MNAR; QP theta = 4; beta = 6,-0.25
  files_WF_MNAR_QP4_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # Putting all the files together
  file_name_str <- c('files_WF_MCAR_6n025', 'files_freqglm0202_MCAR_55n025', 'files_CAR33025_MCAR_6n025', 'files_WF_MCAR_QP4_6n025', 'files_WF_MAR_QP4_6n025','files_WF_MNAR_QP4_6n025')
}

# initialize
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  if(grepl('QP', name)){
    res_district[[name]] <-  get_results(get(name), name, QP_variance_adjustment = 4, district_results = T, results_by_point = T)
  }else{
    res_district[[name]] <-  get_results(get(name), name,  district_results = T, results_by_point = T)
  }
}

save(res_district, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/main_paper_district_by_point_results_06182024.RData')

text_size = 12
## DGP district plots
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_6n025"]]$results, fix_axis = ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F, plot_indiv_points = T, legend_text = text_size)
  
  p2 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_55n025"]]$results, fix_axis = ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F, plot_indiv_points = T, legend_text = text_size)
  
  p3 <- plot_all_methods(res = res_district[["files_CAR33025_MCAR_6n025"]]$results, fix_axis = ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F, plot_indiv_points = T, legend_text = text_size)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_district_bypoint_comparison_plots_06182024.pdf', height = 7.5, width = 10)
}


#
#### Printing all the CAR results! ####
# getting the file names
{
  files = grep('2023_12_09|2023_12_10',dir(res_dir, full.names = T), value = T)
  files2 = grep('2023_12_14',dir(res_dir, full.names = T), value = T)
  
  files_MCAR <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)
  
  files_CAR <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)
  
  files <- unique(c(files, files_MCAR, files_CAR))
  
  for(d in files){
    ff = dir(d)
    if(length(ff) < 51){
      print(d)
      print(length(ff))
      print('-------')
      #unlink(d, recursive = T) # deletes the directory
    }
  }
  # all good
  
  # (1) WF MCAR; poisson; beta = 6, -0.25
  files_WF_MCAR_6n025 <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)
  
  # (2) freqGLM0202 EB0 = 5.5, EB1 = -0.25: MCAR
  files_freqglm0202_MCAR_55n025 <- grep('beta055_beta1n025',grep('freqglm', files2, value = T), value = T)
  
  # (3) CAR0303025 MCAR; poisson; beta = 6, -0.25
  files_CAR33025_MCAR_6n025 <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)
  
  # (4) WF MCAR; QP theta = 4; beta = 6,-0.25
  files_WF_MCAR_QP4_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # (5) WF MAR; QP theta = 4; beta = 6,-0.25
  files_WF_MAR_QP4_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # (6) WF MNAR; QP theta = 4; beta = 6,-0.25
  files_WF_MNAR_QP4_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)
  
  # Putting all the files together
  file_name_str <- c('files_WF_MCAR_6n025', 'files_freqglm0202_MCAR_55n025', 'files_CAR33025_MCAR_6n025', 'files_WF_MCAR_QP4_6n025', 'files_WF_MAR_QP4_6n025','files_WF_MNAR_QP4_6n025')
}

all_file_names <- unlist(sapply(file_name_str, get))
setdiff(files, all_file_names)
setdiff(all_file_names, files)
tt = table(all_file_names); tt[tt>1]
# great

### Processing results and combining
# initialize
res_facility = list()
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  if(grepl('QP', name)){
    res_facility[[name]] <- get_results(get(name), name, QP_variance_adjustment = 4, return_unprocessed = T)  %>% get_CARconvergence()
    res_district[[name]] <-  get_results(get(name), name,  district_results = T, QP_variance_adjustment = 4, return_unprocessed = T)  %>% get_CARconvergence()
  }else{
    res_facility[[name]] <- get_results(get(name), name, return_unprocessed = T) %>% get_CARconvergence()
    res_district[[name]] <-  get_results(get(name), name,  district_results = T, return_unprocessed = T)  %>% get_CARconvergence()
  }
}

# saving CAR convergence from all the main paper runs.
TO DO

#
#### Testing new simulation_main function ####
file <- sprintf('%s/mcar02_freqglm0202_beta055_beta1n025_id91167_2024_06_17',res_dir)

res <- get_results(file, sim_name = 'test', expected_sims = 10)

#
#### Testing CAR convergence ####
res <- combine_results_wrapper(files = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar05_wf_qptheta9_beta06_beta1n025_2023_12_01/', methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARBayes'), return_unprocessed = T)

get_CARconvergence(res)



#
#### WF NB and Poisson Comparison ####

## WF DGP
# load(paste0(res_dir, '/WF_Poisson_NB_comparison_WFDGP_R100_04072024.RData'))
# # 
# # Called all_Model_results and not Imputation_lst. Ugh
# imputed_list <- all_model_results
# save(imputed_list, params, file = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_Poisson_NB_comparison_WFDGP_R100_04072024/sim_results_1(1).RData")

# load(paste0(res_dir, '/WF_Poisson_NB_comparison_NBDGP_R100_04082024.RData'))
# imputed_list <- all_model_results
# save(imputed_list, params, file = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_Poisson_NB_comparison_NBDGP_R100_04082024/sim_results_1(1).RData")

#load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_Poisson_NB_comparison_WFDGP_R100_04082024/sim_results_1(1).RData")

res_WF <- combine_results_wrapper(files = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_Poisson_NB_comparison_WFDGP_R100_04072024', methods = c("y_pred_WF", "y_pred_WF_negbin"), rename_vec = c('WF_Pois','WF_NB'), district_results = F, QP_variance_adjustment = NULL, expected_sims = 100)

res_WF_NB <- combine_results_wrapper(files = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/WF_Poisson_NB_comparison_NBDGP_R100_04082024', methods = c("y_pred_WF", "y_pred_WF_negbin"), rename_vec = c('WF_Pois','WF_NB'), district_results = F, QP_variance_adjustment = NULL, expected_sims = 100, family = 'negbin')

{
WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))

  p1 <- plot_all_methods(res = res_WF$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF Poisson DGP: Empirical Betas', include_legend = F)
  # ok obvi wrong
  
  p2 <- plot_all_methods(res = res_WF_NB$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'WF NB DGP: Empirical Betas', include_legend = F)

  cowplot::plot_grid(p1$plot, p2$plot, p1$legend, ncol = 1, rel_heights = c(3,3,1))
}
ggsave('figures/WF_DGP_Poisson_NB_04102024.png', height = 5, width = 8)

#
#### Making main MDM figure with four panels ####
load(paste0(res_dir,'/full_paper_results_12262023.RData'))
names(res_facility)

WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
## Missing data facility plots
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_WF_MAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_WF_MNAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_facility_plots_fourpanels_03292024_QP_adjust.png', height = 7.5, width = 10)
  # ggsave('figures/WF_missingness_facility_plots_03292024_QP_adjust.pdf', height = 7.5, width = 7)
}

#### Making main DGP facility figure with four panels ####
DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_fourpanels_03222024.png', height = 7.5, width = 10)
}


#
#### Plot MAR Missingness ####
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
  
  
  png(file = 'figures/MAR_heatmap_12282023.png')
  gplots::heatmap.2(1 - tmp2, dendrogram = 'none', Rowv = F, Colv = F, xlab = 'facilities', trace = 'none', key = F)
  dev.off()

  return(p1)
}

lst <- simulate_data(district_sizes = c(10), R = 1, end_date = '2019-12-01')

df <- lst$df_list[[1]] 

# simulation function!
df_miss <- MAR_spatiotemporal_sim(df, p = 0.3, rho = 0.7, alpha = 0.7, tau2 = 16)

p1 <- plot_missingness(df_miss)

# ggsave(p1, path = 'figures/MAR_missing_plot_12282023.png')

#
#### Testing Figures with two parts ####
## DGP facility plots
DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3'), metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3'), metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with CARBayes model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  #ggsave('figures/DGP_facility_comparison_plots_12192023.png', height = 7.5, width = 10)
}

#
#### Updating results names #### 
load(paste0(res_dir,'/full_paper_results_12192023.RData'))
names(res_facility)

for(i in 1:length(res_facility)){
  res_facility[[i]]$results$method <- factor(gsub('CARBayes','CAR', as.character(res_facility[[i]]$results$method)), levels = c('WF','freqGLM','CAR'))
  res_district[[i]]$results$method <- factor(gsub('CARBayes','CAR', as.character(res_district[[i]]$results$method)), levels = c('WF','freqGLM','CAR'))
}

# save(res_facility, res_district, file = paste0(res_dir,'/full_paper_results_12262023.RData'))

#
#### Making appendix plots ####
#load(paste0(res_dir,'/full_paper_results_12192023.RData'))
load(paste0(res_dir,'/full_paper_results_12262023.RData'))
names(res_facility)

# facility DGP B0 = 6/5.5: B1 = 0
DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_beta6_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_beta6_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_beta6_0_12192023.pdf', height = 7.5, width = 10)
}

# district DGP B0 = 6/5.5: B1 = 0
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_beta6_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_beta55_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_CAR33025_MCAR_beta6_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_district_comparison_plots_beta6_0_12192023.pdf', height = 7.5, width = 10)
}

# facility DGP B0 = 2/1.5: B1 = 0
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_beta2_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta15_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_beta2_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_beta2_0_12192023.pdf', height = 7.5, width = 10)
}

# district DGP B0 = 2/1.5: B1 = 0
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_beta2_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_beta15_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_CAR33025_MCAR_beta2_0"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_district_comparison_plots_beta2_0_12192023.pdf', height = 7.5, width = 10)
}

WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
# facility MGP: B0 = 6: B1 = 0: Theta = 4
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_WF_MAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_WF_MNAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_facility_beta6_0_QPtheta4_12192023.pdf', height = 7.5, width = 10)
}

# district MGP: B0 = 6: B1 = 0: Theta = 4
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_WF_MAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_WF_MNAR_QPtheta4_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_district_beta6_0_QPtheta4_12192023.pdf', height = 7.5, width = 10)
}

# facility MGP: B0 = 6: B1 = 0: Theta = 16
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_WF_MAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_WF_MNAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_facility_beta6_0_QPtheta16_12192023.pdf', height = 7.5, width = 10)
}

# district MGP: B0 = 6: B1 = 0: Theta = 16
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_WF_MAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_WF_MNAR_QPtheta16_beta6_0"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_district_beta6_0_QPtheta16_12192023.pdf', height = 7.5, width = 10)
}

#
#### Combining main results and appendix together ####
load(paste0(res_dir,'/appendix_results_12142023.RData'))
res_facility_new = res_facility; res_district_new = res_district
rm(res_facility, res_district)
load(paste0(res_dir,'/main_paper_results_12102023_QP_variance_adjusted.RData'))
res_facility_old = res_facility; res_district_old = res_district
rm(res_facility, res_district)

res_facility = c(res_facility_old, res_facility_new)
res_district = c(res_district_old, res_district_new)

table(table(names(res_facility)))
table(table(names(res_district)))
# good

# save(res_facility, res_district, file = paste0(res_dir,'/full_paper_results_12192023.RData'))

#
#### Comparing freq0550202 with freq050202 ####
load(paste0(res_dir,'/full_paper_results_12192023.RData'))

## Comparison  plots
DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'beta = 5: Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'beta = 5.5: Data Generated with CARBayes model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p1$legend, ncol = 1, rel_heights = c(3,3,1))
  #ggsave('figures/DGP_facility_comparison_plots_12102023.png', height = 7.5, width = 10)
}

{
  p1 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'beta = 5: District level', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), results_by_point = F, rows = 1, title = 'beta = 5.5: District level', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p1$legend, ncol = 1, rel_heights = c(3,3,1))
  #ggsave('figures/DGP_facility_comparison_plots_12102023.png', height = 7.5, width = 10)
}

#
#### Getting newer results: Appendix and freqGLM0550202 ####
files = grep('2023_12_14',dir(res_dir, full.names = T), value = T)

# freqGLM0202 EB0 = 5.5, EB1 = -0.25: MCAR (not actually the appendix though)
files_freqglm0202_MCAR_beta55_n025 <- grep('beta055_beta1n025',grep('freqglm', files, value = T), value = T)

#	freqGLM0202 EB0 = 5.5, EB1 = 0: MCAR
files_freqglm0202_MCAR_beta55_0 <- grep('beta055_beta10',grep('freqglm', files, value = T), value = T)

#	WF EB0 = 6, EB1 = 0: MCAR
files_WF_MCAR_beta6_0 <- grep('wf_beta06', files, value = T)

#	CAR0303025 EB0 = 6, EB1 = 0: MCAR
files_CAR33025_MCAR_beta6_0 <- grep('car0303025_beta06', files, value = T)

#	freqGLM0202 EB0 = 1.5, EB1 = 0: MCAR
files_freqglm0202_MCAR_beta15_0 <- grep('freqglm0202_beta015', files, value = T)

#	WF EB0 = 2, EB1 = 0: MCAR
files_WF_MCAR_beta2_0 <- grep('wf_beta02', files, value = T)

#	CAR0303025 EB0 = 2, EB1 = 0: MCAR
files_CAR33025_MCAR_beta2_0 <- grep('car0303025_beta02', files, value = T)

#	WF EB0 = 6, EB1 = 0: QP 4: MCAR
files_WF_MCAR_QPtheta4_beta6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta4', files, value = T)

#	WF EB0 = 6, EB1 = 0: QP 4: MAR
files_WF_MAR_QPtheta4_beta6_0 <- grep('mar[0-9]{1,2}_wf_qptheta4', files, value = T)

# WF EB0 = 6, EB1 = 0: QP 4: MNAR
files_WF_MNAR_QPtheta4_beta6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta4', files, value = T)

#	WF EB0 = 6, EB1 = 0: QP 16: MCAR
files_WF_MCAR_QPtheta16_beta6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta16', files, value = T)

#	WF EB0 = 6, EB1 = 0: QP 16: MAR
files_WF_MAR_QPtheta16_beta6_0 <- grep('mar[0-9]{1,2}_wf_qptheta16', files, value = T)

# WF EB0 = 6, EB1 = 0: QP 16: MNAR
files_WF_MNAR_QPtheta16_beta6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta16', files, value = T)

for(d in files){
  ff = dir(d)
  if(length(ff) < 51){
    print(d)
    print(length(ff))
    print('-------')
    #unlink(d, recursive = T) # deletes the directory
  }
}

# Putting all the files together
file_name_str <- c('files_freqglm0202_MCAR_beta55_n025', 'files_freqglm0202_MCAR_beta55_0', 'files_WF_MCAR_beta6_0', 'files_CAR33025_MCAR_beta6_0', 'files_freqglm0202_MCAR_beta15_0', 'files_WF_MCAR_beta2_0', 'files_CAR33025_MCAR_beta2_0', 'files_WF_MCAR_QPtheta4_beta6_0', 'files_WF_MAR_QPtheta4_beta6_0', 'files_WF_MNAR_QPtheta4_beta6_0', 'files_WF_MCAR_QPtheta16_beta6_0', 'files_WF_MAR_QPtheta16_beta6_0', 'files_WF_MNAR_QPtheta16_beta6_0')

all_file_names <- unlist(sapply(file_name_str, get))
setdiff(files, all_file_names)
setdiff(all_file_names, files)
tt = table(all_file_names); tt[tt>1]
# ok all gravy

# initialize
res_facility = list()
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  res_facility[[name]] <- get_results(get(name), name)
  res_district[[name]] <-  get_results(get(name), name,  district_results = T)
}

# save(res_facility, res_district, file = paste0(res_dir,'/appendix_results_12142023.RData'))


#
#### Get the aggregated results numbers (updated 12262023) ####
#load(paste0(res_dir,'/main_paper_results_12102023_QP_variance_adjusted.RData'))
# load(paste0(res_dir,'/full_paper_results_12192023.RData'))
load(paste0(res_dir,'/full_paper_results_12262023.RData'))

res_f_agg <- res_f_agg2 <- list()
names(res_facility)
for(sim in names(res_facility)){
  res_f_agg[[sim]] <- aggregate_results(res_facility[[sim]]$results)
  res_f_agg2[[sim]] <- aggregate_results(res_facility[[sim]]$results, full_results = F, metrics = c('specificity', 'outbreak_detection3'), metric_rename = c('specificity', 'sensitivity-3'))
}

res_d_agg <- res_d_agg2 <- list()
for(sim in names(res_district)){
  res_d_agg[[sim]] <- aggregate_results(res_district[[sim]]$results)
  res_d_agg2[[sim]] <- aggregate_results(res_district[[sim]]$results, full_results = F)
}

# sensitivity-10
for(xx in res_f_agg){
  print(xx[[4]])
}

for(xx in res_d_agg){
  print(xx[[4]])
}

# sensitivity-5
for(xx in res_f_agg){
  print(names(xx))
  print(xx[[3]])
}

for(xx in res_d_agg){
  print(names(xx))
  print(xx[[3]])
}

### Appendix 1
res_f_agg2[["files_WF_MCAR_beta2_0"]]
res_f_agg2[['files_freqglm0202_MCAR_beta15_0']]
res_f_agg2[['files_CAR33025_MCAR_beta2_0']]

res_f_agg2[["files_WF_MCAR_beta6_0"]]
res_f_agg2[['files_freqglm0202_MCAR_beta55_0']]
res_f_agg2[['files_CAR33025_MCAR_beta6_0']]

### Appendix 2
res_f_agg2[["files_WF_MCAR_QP4_6n025"]]
res_f_agg2[['files_WF_MCAR_QPtheta16_beta6_0']]

res_f_agg2[["files_WF_MAR_QP4_6n025"]]
res_f_agg2[['files_WF_MAR_QPtheta16_beta6_0']]

res_f_agg2[["files_WF_MNAR_QP4_6n025"]]
res_f_agg2[['files_WF_MNAR_QPtheta16_beta6_0']]

### Main
res_f_agg[["files_WF_MCAR_6n025"]]
res_f_agg[["files_freqglm0202_MCAR_6n025"]]
res_f_agg[["files_freqglm0202_MCAR_beta55_n025"]]
res_f_agg[["files_CAR33025_MCAR_6n025"]]

res_d_agg[["files_freqglm0202_MCAR_6n025"]]
res_d_agg[["files_freqglm0202_MCAR_beta55_n025"]]

res_f_agg[["files_WF_MCAR_QP4_6n025"]][['specificity']]
res_f_agg[["files_WF_MAR_QP4_6n025"]][['specificity']]
res_f_agg[["files_WF_MNAR_QP4_6n025"]][['specificity']]

res_f_agg[["files_WF_MCAR_QP4_6n025"]][['sensitivity-3']]
res_f_agg[["files_WF_MAR_QP4_6n025"]][['sensitivity-3']]
res_f_agg[["files_WF_MNAR_QP4_6n025"]][['sensitivity-3']]
#
#### Adjust the variance on QP runs where it wasn't accounted for ####

res = get_results(files[1], 'test', QP_variance_adjustment = 4)
res2 = get_results(files[1], 'test')

a = res$results; b = res2$results
sum(a$specificity == b$specificity)
sum(a$outbreak_detection3 == b$outbreak_detection3)
sum(a$outbreak_detection3 > b$outbreak_detection3)
sum(a$outbreak_detection3 < b$outbreak_detection3)

sum(a$outbreak_detection5 == b$outbreak_detection10)
# perfecto


#
#### Plot main paper results 12/10/2023 ####
# load(paste0(res_dir,'/main_paper_results_12102023.RData'))
# load(paste0(res_dir,'/main_paper_results_12102023_QP_variance_adjusted.RData'))
#load(paste0(res_dir,'/full_paper_results_12192023.RData'))
load(paste0(res_dir,'/full_paper_results_12262023.RData'))
names(res_facility)

WF_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
## Missing data facility plots
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_WF_MAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_WF_MNAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  #ggsave('figures/WF_missingness_facility_plots_12102023_QP_adjust.png', height = 7.5, width = 10)
 # ggsave('figures/WF_missingness_facility_plots_12102023_QP_adjust.pdf', height = 7.5, width = 7)
}

## Missing data district plots
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MCAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_WF_MAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_WF_MNAR_QP4_6n025"]]$results, fix_axis = WF_ylims, add_lines = list(F, F, F, F), metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'MNAR, with Data Generated by WF Quasi-Poisson', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/WF_missingness_district_plots_12102023_QP_adjust.pdf', height = 7.5, width = 7)
}

## DGP facility plots
DGP_ylims = list(ylim(0,1), ylim(0,1), ylim(0, 1), ylim(0,1))
{
  p1 <- plot_all_methods(res = res_facility[["files_WF_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_facility[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_facility[["files_CAR33025_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_facility_comparison_plots_12192023.pdf', height = 7.5, width = 7)
}

## DGP district plots
{
  p1 <- plot_all_methods(res = res_district[["files_WF_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with WF model, assuming MCAR', include_legend = F)
  
  p2 <- plot_all_methods(res = res_district[["files_freqglm0202_MCAR_beta55_n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with freqGLM model, assuming MCAR', include_legend = F)
  
  p3 <- plot_all_methods(res = res_district[["files_CAR33025_MCAR_6n025"]]$results, fix_axis = DGP_ylims, add_lines = list(F, F, F, F),  metrics = c('specificity', 'outbreak_detection3'),  metric_rename = c('specificity', 'sensitivity-3'), results_by_point = F, rows = 1, title = 'Data Generated with CAR model, assuming MCAR', include_legend = F)
  
  cowplot::plot_grid(p1$plot, p2$plot, p3$plot, p1$legend, ncol = 1, rel_heights = c(3,3,3,1))
  ggsave('figures/DGP_district_comparison_plots_12192023.pdf', height = 7.5, width = 7)
}


#
#### Processing results 12/10/2023 ####
files = grep('2023_12_09|2023_12_10',dir(res_dir, full.names = T), value = T)

files_MCAR <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)

files_CAR <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)

files <- unique(c(files, files_MCAR, files_CAR))

# ideally, should print nothing.
for(d in files){
  ff = dir(d)
  if(length(ff) < 51){
    print(d)
    print(length(ff))
    print('-------')
    #unlink(d, recursive = T) # deletes the directory
  }
}
# all good

# (1) WF MCAR; poisson; beta = 6, -0.25
files_WF_MCAR_6n025 <- grep('2023_10_19',dir(res_dir, full.names = T), value = T)

# (2) freqglm0202 MCAR; poisson; beta = 6,-0.25
files_freqglm0202_MCAR_6n025 <- grep('mcar[0-9]{1,2}_freqglm0202_beta05_beta1n025', files, value = T)

# (3) CAR0303025 MCAR; poisson; beta = 6, -0.25
files_CAR33025_MCAR_6n025 <- grep('car33025', grep('2023_10_21',dir(res_dir, full.names = T), value = T), value = T)

# (4) WF MCAR; QP theta = 4; beta = 6,-0.25
files_WF_MCAR_QP4_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)

# (5) WF MAR; QP theta = 4; beta = 6,-0.25
files_WF_MAR_QP4_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)

# (6) WF MNAR; QP theta = 4; beta = 6,-0.25
files_WF_MNAR_QP4_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta4_beta06_beta1n025', files, value = T)

# Putting all the files together
file_name_str <- c('files_WF_MCAR_6n025', 'files_freqglm0202_MCAR_6n025', 'files_CAR33025_MCAR_6n025', 'files_WF_MCAR_QP4_6n025', 'files_WF_MAR_QP4_6n025','files_WF_MNAR_QP4_6n025')

all_file_names <- unlist(sapply(file_name_str, get))
setdiff(files, all_file_names)
setdiff(all_file_names, files)
tt = table(all_file_names); tt[tt>1]
# great

### Processing results and combining
# initialize
res_facility = list()
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  if(grepl('QP', name)){
    res_facility[[name]] <- get_results(get(name), name, QP_variance_adjustment = 4)
    res_district[[name]] <-  get_results(get(name), name,  district_results = T, QP_variance_adjustment = 4)
  }else{
    res_facility[[name]] <- get_results(get(name), name)
    res_district[[name]] <-  get_results(get(name), name,  district_results = T)
  }
}

# save(res_facility, res_district, file = paste0(res_dir,'/main_paper_results_12102023_QP_variance_adjusted.RData'))

# Checking for NA vals
lapply(res_facility, function(xx){sum(is.na(xx$params))})
lapply(res_district, function(xx){sum(is.na(xx$params))})

# combine results into data frames
res_df <- do.call('rbind',lapply(res_facility, function(xx) xx$results))
res_district_df <- do.call('rbind',lapply(res_district, function(xx) xx$results))

# test the table of simulation types
table(res_df$sim)
table(res_district_df$sim)
table(res_df$sim, res_df$prop_missing)
table(res_district_df$sim, res_district_df$prop_missing)
# Ok all good

# test the uniqueness of rows (independent of names) - they should be unique
length(unique(res_df$relative_bias)) # Bueno

# are they unique?
tt = table(res_df$relative_bias)
table(tt) 
# Yay bueno.

# Should get 3 sets of these to be the same, for WF MCAR, MNAR, MAR Poisson, right?

#
#### Testing results from 12/10/2023 ####
load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar01_wf_qptheta4_beta06_beta1n025_id147928_2023_12_09/sim_results_p0.1_mcar_1(50).RData')
imputed_list1 <- imputed_list; params1 <- params

load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mcar0_wf_beta6_n025_2023_10_19/sim_results_p0.0_mcar_1(50).RData")
imputed_list2 <- imputed_list; params2 <- params

sum(imputed_list1[[1]]$df_miss$y, na.rm = T)
sum(imputed_list2[[1]]$df_miss$y, na.rm = T)

#uh-oh there's a params issue. Ending with \r. No bueno

#
#### Processing results from test run 12/08/2023 ####
files = grep('2023_12_08',dir(res_dir, full.names = T), value = T)

load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar0_wf_beta06_beta1n025_id326672_2023_12_08/sim_results_p0.0_mar_1(50).RData")
imputed_list1 <- imputed_list; params1 <- params

load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar0_wf_qptheta100_beta06_beta1n025_id120292_2023_12_08/sim_results_p0.0_mar_1(50).RData")
imputed_list100 <- imputed_list; params100 <- params

load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar0_wf_qptheta9_beta06_beta1n025_id530083_2023_12_08/sim_results_p0.0_mar_1(50).RData")
imputed_list9 <- imputed_list; params9 <- params

sum(imputed_list1[[1]]$df_miss$y)
sum(imputed_list100[[1]]$df_miss$y)
sum(imputed_list9[[1]]$df_miss$y)
# ok all different. Good. Good.

res1 <- combine_results_wrapper(files = files[1], methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

res100 <- combine_results_wrapper(files = files[2], methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

res9 <- combine_results_wrapper(files = files[3], methods = c("y_pred_WF", "y_pred_freqGLMepi", 'y_CARstan'), rename_vec = c('WF','freqGLM', 'CARstan'))

p1 <- plot_all_methods(res = res1$results, rows = 1)
p2 <- plot_all_methods(res = res9$results, rows = 1)
p3 <- plot_all_methods(res = res100$results, rows = 1)

cowplot::plot_grid(p1, p2, p3)
# ok it can't do the cowplot but whatever. This is clearly different. So then we're all good, right? Jesus Christo.

#
#### Comparing the results from two different data sets ####
files = grep('2023_12_05',dir(res_dir, full.names = T), value = T)
files_WF_MNAR_QP9_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)
files_WF_MNAR_QP100_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar0_wf_qptheta100_beta06_beta1n025_id798896_2023_12_05/sim_results_p0.0_mnar_1(50).RData')
imputed_list100 <- imputed_list; rm(imputed_list)

load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mnar0_wf_qptheta9_beta06_beta1n025_id738974_2023_12_05/sim_results_p0.0_mnar_1(50).RData')
imputed_list9 <- imputed_list; rm(imputed_list)

a <- imputed_list100[[1]]$df_miss
b <- imputed_list9[[1]]$df_miss

sum(a$y == b$y) # identical!

# is QP happening at all?
a$residual = abs(a$y - a$y_exp)

plot(a$y, a$residual)
lines(1:1500, sqrt(1:1500), col = 'red')

# no QP. WTF.

#
#### Processing the 12/05 results (I know some will be missing) ####
# get file names
files = grep('2023_12_05',dir(res_dir, full.names = T), value = T)
# 81 (after deleting incomplete ones)

# for each directory, scan for the # of files and delete improper ones
for(d in files){
  ff = dir(d)
  if(length(ff) < 51){
    print(d)
    print(length(ff))
    print('-------')
    #unlink(d, recursive = T) # deletes the directory
  }
}
# all good!

# subsets of file names
# (1) WF MAR QP 9 (beta = 6, -0.25)
files_WF_MAR_QP9_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

# (2) WF MAR QP 100 (beta = 6, -0.25)
files_WF_MAR_QP100_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

# (3) WF MNAR QP 9 (beta = 6, -0.25)
files_WF_MNAR_QP9_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

# (4) WF MNAR QP 100 (beta = 6, -0.25)
files_WF_MNAR_QP100_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

# (5) WF MCAR QP 9 (beta = 6, -0.25)
files_WF_MCAR_QP9_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

# (6) WF MCAR QP 100 (beta = 6, -0.25)
files_WF_MCAR_QP100_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

# (7) WF MAR QP 9 (beta = 6, 0)
files_WF_MAR_QP9_6_0 <- grep('mar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)

# (8) WF MAR QP 100 (beta = 6, 0)
files_WF_MAR_QP100_6_0 <- grep('mar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)

# (9) WF MNAR QP 9 (beta = 6, 0)
files_WF_MNAR_QP9_6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)

# (10) WF MNAR QP 100 (beta = 6, 0)
files_WF_MNAR_QP100_6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)

# (11) WF MCAR QP 9 (beta = 6, 0)
files_WF_MCAR_QP9_6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)

# (12) WF MCAR QP 100 (beta = 6, 0)
files_WF_MCAR_QP100_6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)

# (13) WF MCAR (beta = 2, 0)
files_WF_MCAR_2_0 <- grep('mcar[0-9]{1,2}_wf_beta02_beta10', files, value = T)

# (14) freqGLM MCAR (beta = 2, 0)
files_freqGLM_MCAR_2_0 <- grep('mcar[0-9]{1,2}_freqglm0202_beta02_beta10', files, value = T)

# (15) CAR MCAR (beta = 2, 0)
files_CAR_MCAR_2_0 <- grep('mcar[0-9]{1,2}_car0303025_beta02_beta10', files, value = T)

# QC the file names
file_names <- c(files_WF_MAR_QP9_6n025, files_WF_MAR_QP100_6n025, files_WF_MNAR_QP9_6n025, files_WF_MNAR_QP100_6n025, files_WF_MCAR_QP9_6n025, files_WF_MCAR_QP100_6n025, files_WF_MAR_QP9_6_0, files_WF_MAR_QP100_6_0, files_WF_MNAR_QP9_6_0, files_WF_MNAR_QP100_6_0, files_WF_MCAR_QP9_6_0, files_WF_MCAR_QP100_6_0, files_WF_MCAR_2_0, files_freqGLM_MCAR_2_0, files_CAR_MCAR_2_0)

setdiff(files, file_names)
tt = table(file_names); tt[tt>1]
# great

### Getting all the results together
file_name_str <- c('files_WF_MAR_QP9_6n025', 'files_WF_MAR_QP100_6n025', 'files_WF_MNAR_QP9_6n025', 'files_WF_MNAR_QP100_6n025', 'files_WF_MCAR_QP9_6n025', 'files_WF_MCAR_QP100_6n025', 'files_WF_MAR_QP9_6_0', 'files_WF_MAR_QP100_6_0', 'files_WF_MNAR_QP9_6_0', 'files_WF_MNAR_QP100_6_0', 'files_WF_MCAR_QP9_6_0', 'files_WF_MCAR_QP100_6_0', 'files_WF_MCAR_2_0', 'files_freqGLM_MCAR_2_0', 'files_CAR_MCAR_2_0')

# initialize
res_facility = list()
res_district = list()

# grab the results from all runs
for(name in file_name_str){
  res_facility[[name]] <- get_results(get(name), name)
  res_district[[name]] <-  get_results(get(name), name,  district_results = T)
}

# save(res_facility, res_district_df, file = paste0(res_dir,'/tmp_results_12082023.RData'))

# Checking for NA vals
lapply(res_facility, function(xx){sum(is.na(xx$params))})
lapply(res_district, function(xx){sum(is.na(xx$params))})

# combine results into data frames
res_df <- do.call('rbind',lapply(res_facility, function(xx) xx$results))
res_district_df <- do.call('rbind',lapply(res_district, function(xx) xx$results))

# test the table of simulation types
table(res_df$sim)
table(res_district_df$sim)
table(res_df$sim, res_df$prop_missing)
table(res_district_df$sim, res_district_df$prop_missing)
# Ok so some glitches with CAR 2_0. 9/3000 at p = 0.5. Whatever

# test the uniqueness of rows (independent of names) - they should be unique
length(res_df$RMSE)
length(unique(res_df$RMSE)) # hmm
length(unique(res_df$prop_interval_width))
length(unique(res_df$relative_bias)) # That's a concern

# are they unique?
tt = table(res_df$relative_bias)
table(tt) 

options(digits = 22)
vals = as.numeric(names(tt[tt==6]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

res_df[ind,]

# so what needs to get rerun?

vals = as.numeric(names(tt[tt==5]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

res_df[ind,]

vals = as.numeric(names(tt[tt==2]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

res_df[ind,]
# HOWWWWW! Oh my god. Jesus.

for(v in sample(vals, 20)){
  ind = which(abs(res_df$relative_bias - v) < 1e-10)
  print(res_df[ind,'prop_missing'])
}

vals = as.numeric(names(tt[tt==1]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

for(v in sample(vals, 20)){
  ind = which(abs(res_df$relative_bias - v) < 1e-10)
  print(res_df[ind,'sim'])
}

#
#### Processing all the 12/01 and 12/04 results ####
files = grep('2023_12_01|2023_12_04',dir(res_dir, full.names = T), value = T)

# initialize
res_facility = list()
res_district = list()

# (1) WF MAR QP 9 (beta = 6, -0.25)
files_MAR_QP9_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

res_facility[['WF_MAR_QP9_6n025']] <- get_results(files_MAR_QP9_6n025, 'WF_MAR_QP9_6n025')
res_district[['WF_MAR_QP9_6n025']] <-  get_results(files_MAR_QP9_6n025, 'WF_MAR_QP9_6n025', district_results = T)

# (2) WF MAR QP 100 (beta = 6, -0.25)
files_MAR_QP100_6n025 <- grep('mar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

res_facility[['WF_MAR_QP100_6n025']] <- get_results(files_MAR_QP100_6n025, 'WF_MAR_QP100_6n025')
res_district[['WF_MAR_QP100_6n025']] <-  get_results(files_MAR_QP100_6n025, 'WF_MAR_QP100_6n025', district_results = T)

# (3) WF MNAR QP 9 (beta = 6, -0.25)
files_MNAR_QP9_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

res_facility[['WF_MNAR_QP9_6n025']] <- get_results(files_MNAR_QP9_6n025, 'WF_MNAR_QP9_6n025')
res_district[['WF_MNAR_QP9_6n025']] <-  get_results(files_MNAR_QP9_6n025, 'WF_MNAR_QP9_6n025', district_results = T)

# (4) WF MNAR QP 100 (beta = 6, -0.25)
files_MNAR_QP100_6n025 <- grep('mnar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

res_facility[['WF_MNAR_QP100_6n025']] <- get_results(files_MNAR_QP100_6n025, 'WF_MNAR_QP100_6n025')
res_district[['WF_MNAR_QP100_6n025']] <-  get_results(files_MNAR_QP100_6n025, 'WF_MNAR_QP100_6n025', district_results = T)

# (5) WF MCAR QP 9 (beta = 6, -0.25)
files_MCAR_QP9_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta9_beta06_beta1n025', files, value = T)

res_facility[['WF_MCAR_QP9_6n025']] <- get_results(files_MCAR_QP9_6n025, 'WF_MCAR_QP9_6n025')
res_district[['WF_MCAR_QP9_6n025']] <-  get_results(files_MCAR_QP9_6n025, 'WF_MCAR_QP9_6n025', district_results = T)

# (6) WF MCAR QP 100 (beta = 6, -0.25)
files_MCAR_QP100_6n025 <- grep('mcar[0-9]{1,2}_wf_qptheta100_beta06_beta1n025', files, value = T)

res_facility[['WF_MCAR_QP100_6n025']] <- get_results(files_MCAR_QP100_6n025, 'WF_MCAR_QP100_6n025')
res_district[['WF_MCAR_QP100_6n025']] <-  get_results(files_MCAR_QP100_6n025, 'WF_MCAR_QP100_6n025', district_results = T)

# # (7) WF MAR QP 9 (beta = 6, 0)
# files_MAR_QP9_6_0 <- grep('mar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)
# 
# res_facility[['WF_MAR_QP9_6_0']] <- get_results(files_MAR_QP9_6_0, 'WF_MAR_QP9_6_0')
# res_district[['WF_MAR_QP9_6_0']] <-  get_results(files_MAR_QP9_6_0, 'WF_MAR_QP9_6_0', district_results = T)
# 
# # (8) WF MAR QP 100 (beta = 6, 0)
# files_MAR_QP100_6_0 <- grep('mar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)
# 
# res_facility[['WF_MAR_QP100_6_0']] <- get_results(files_MAR_QP100_6_0, 'WF_MAR_QP100_6_0')
# res_district[['WF_MAR_QP100_6_0']] <-  get_results(files_MAR_QP100_6_0, 'WF_MAR_QP100_6_0', district_results = T)

# (9) WF MNAR QP 9 (beta = 6, 0)
files_MNAR_QP9_6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)

res_facility[['WF_MNAR_QP9_6_0']] <- get_results(files_MNAR_QP9_6_0, 'WF_MNAR_QP9_6_0')
res_district[['WF_MNAR_QP9_6_0']] <-  get_results(files_MNAR_QP9_6_0, 'WF_MNAR_QP9_6_0', district_results = T)

# (10) WF MNAR QP 100 (beta = 6, 0)
files_MNAR_QP100_6_0 <- grep('mnar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)

res_facility[['WF_MNAR_QP100_6_0']] <- get_results(files_MNAR_QP100_6_0, 'WF_MNAR_QP100_6_0')
res_district[['WF_MNAR_QP100_6_0']] <-  get_results(files_MNAR_QP100_6_0, 'WF_MNAR_QP100_6_0', district_results = T)

# (11) WF MCAR QP 9 (beta = 6, 0)
files_MCAR_QP9_6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta9_beta06_beta10', files, value = T)

res_facility[['WF_MCAR_QP9_6_0']] <- get_results(files_MCAR_QP9_6_0, 'WF_MCAR_QP9_6_0')
res_district[['WF_MCAR_QP9_6_0']] <-  get_results(files_MCAR_QP9_6_0, 'WF_MCAR_QP9_6_0', district_results = T)

# (12) WF MCAR QP 100 (beta = 6, 0)
files_MCAR_QP100_6_0 <- grep('mcar[0-9]{1,2}_wf_qptheta100_beta06_beta10', files, value = T)

res_facility[['WF_MCAR_QP100_6_0']] <- get_results(files_MCAR_QP100_6_0, 'WF_MCAR_QP100_6_0')
res_district[['WF_MCAR_QP100_6_0']] <-  get_results(files_MCAR_QP100_6_0, 'WF_MCAR_QP100_6_0', district_results = T)

# (13) WF MCAR (beta = 2, 0)
files_WF_MCAR_2_0 <- grep('mcar[0-9]{1,2}_wf_beta02_beta10', files, value = T)

# (14) freqGLM MCAR (beta = 2, 0)
files_freqGLM_MCAR_2_0 <- grep('mcar[0-9]{1,2}_freqglm0202_beta02_beta10', files, value = T)

# (15) CAR MCAR (beta = 2, 0)
files_CAR_MCAR_2_0 <- grep('mcar[0-9]{1,2}_car0303025_beta02_beta10', files, value = T)

### Combine the results
lapply(res_facility, function(xx){sum(is.na(xx$params))})
lapply(res_district, function(xx){sum(is.na(xx$params))})

res_df <- do.call('rbind',lapply(res_facility, function(xx) xx$results))
res_district_df <- do.call('rbind',lapply(res_district, function(xx) xx$results))

### Unit testing
# combine all the file names. Test the uniqueness
file_names <- c(files_MAR_QP9_6n025, files_MAR_QP100_6n025, files_MNAR_QP9_6n025, files_MNAR_QP100_6n025, files_MCAR_QP9_6n025, files_MCAR_QP100_6n025, files_MAR_QP9_6_0, files_MAR_QP100_6_0, files_MNAR_QP9_6_0, files_MNAR_QP100_6_0, files_MCAR_QP9_6_0, files_MCAR_QP100_6_0, files_MCAR_2_0)

tt = table(file_names); tt[tt>1]

# test the table of simulation types
table(res_df$sim)
table(res_district_df$sim)
table(res_df$sim, res_df$prop_missing)
table(res_district_df$sim, res_district_df$prop_missing)

# test the uniqueness of rows (indepedent of names) - they should be unique
length(res_df$RMSE)
length(unique(res_df$RMSE)) # hmm
length(unique(res_df$prop_interval_width))
length(unique(res_df$relative_bias)) # 87000....

tt = table(res_df$relative_bias)
table(tt) # hm ok there's something not right here. Ah is it the p = 0! Yes it should be.
options(digits = 22)
vals = as.numeric(names(tt[tt==6]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

res_df[ind,]
# interesting. I get the p = 0. But why QP = 9 and 100?

vals = as.numeric(names(tt[tt==2]))
ind = which(abs(res_df$relative_bias - vals[1]) < 1e-10)

res_df[ind,]
# no difference! Argh. What is going on?

tmp <- res_df %>% filter(sim %in% c('WF_MNAR_QP9_6n025', 'WF_MNAR_QP100_6n025'))
table(table(tmp$relative_bias)) # all the exact same!! Wtf?

vals = as.numeric(names(tt[tt==4]))
ind = which(abs(res_df$relative_bias - vals[3000]) < 1e-10)

res_df[ind,]

#
#### Processing all the 12/01 results ####
# So I named the output files wrong. Damn
files = grep('2023_12_01',dir(res_dir, full.names = T), value = T)

# check that there are enough files (as in, none overwrote each other)
for(d in files){
  tt = dir(d)
  print(length(tt))
}

# check the parameters of each file in each folder
param_list = list()

# cycle through each file, pull the params, and store
for(d in files){
  tt = dir(d, full.names = T)
  tt = tt[-grep('simulated_data', tt)]
  
  df = data.frame(root_dir = NULL, 
                  parenthses = NULL,
                  theta = NULL)
  for(file in tt){
    load(file)
    df = rbind(df, data.frame(
      root_dir = d,
      parentheses = str_count(file, '\\('),
      theta = params$theta
    ))
  }
  
  param_list[[d]] = df
}

# make sure the params are distinctly separate
for(d in names(param_list)){
  print(d)
  xx = param_list[[d]]
  print(table(xx$parentheses, xx$theta))
}

# delete the ish
for(d in files){
  if(grepl('theta100', d)){
    # remove files with 1 parentheses
    for(file in dir(d, full.names = T)){
      parentheses = str_count(file, '\\(')
      if(parentheses == 1){
        file.remove(file)
      }
    }
    # rename the files
    for(file in dir(d, full.names = T)){
      new_file = gsub('\\)\\([0-9]\\)','\\)',file)
      file.rename(file, new_file)
    }
  }else{
    for(file in dir(d, full.names = T)){
      parentheses = str_count(file, '\\(')
      if(parentheses == 2){
        file.remove(file)
      }
    }
  }
}

# run the parameter check again to confirm this was ok. 
# check the parameters of each file in each folder
param_list = list()

# cycle through each file, pull the params, and store
for(d in files){
  tt = dir(d, full.names = T)
  tt = tt[-grep('simulated_data', tt)]
  
  df = data.frame(root_dir = NULL, 
                  parenthses = NULL,
                  theta = NULL)
  for(file in tt){
    load(file)
    df = rbind(df, data.frame(
      root_dir = d,
      parentheses = str_count(file, '\\('),
      theta = params$theta
    ))
  }
  
  param_list[[d]] = df
}

# make sure the params are distinctly separate
for(d in names(param_list)){
  print(d)
  xx = param_list[[d]]
  print(table(xx$parentheses, xx$theta))
}

# good! Phew! Jesus.

#
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


