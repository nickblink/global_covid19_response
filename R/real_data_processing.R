library(MASS)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)
library(cowplot)

source('R/imputation_functions.R')

if(file.exists('C:/Users/Admin-Dell')){
  res_dir = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}else{
  res_dir = "C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}


generate_outbreak_columns <- function(df, models = c('y_pred_WF_negbin','y_pred_freqGLMepi_negbin','y_CAR_phifit'), rename_vec = NULL){
  
  # setting names of outbreak cols.
  if(is.null(rename_vec)){
    rename_vec <- models
  }
  
  # making outbreak columns
  for(i in 1:length(models)){
    col <- paste0(models[i], '_0.975')
    out_col <- paste0('outbreak_', rename_vec[i])
    print(col)
    df[,out_col] <- ifelse(df$y > df[,col], 1, 0)
  }
  
  return(df)
}


#### 05/12/2025: Getting CAR stan diagnostics ####
load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_07232024.rdata')

# They were fit on rolling windows for 11 months so there are 11 model fits.
CAR_res <- lapply(results_list, function(xx){
  tt <- xx$res_list$CAR_phifit_negbin$CARstan_summary
  as.data.frame(tt)
})

# Step 1: Define a function to tag each parameter with a group
label_groups <- function(df) {
  df$param <- rownames(df)
  
  df$group <- "other"
  df$group[df$param == "tau2"] <- "tau2"
  df$group[df$param == "rho"] <- "rho"
  df$group[df$param == "alpha"] <- "alpha"
  df$group[grepl("^theta\\[[0-9]+\\]$", df$param)] <- "theta"
  
  return(df[, c("group", "n_eff", "Rhat")])
}

# Step 2: Apply the labeling function to each df and combine them
all_data <- do.call(rbind, lapply(df_list, label_groups))

# Step 3: Compute the median across all data frames for each group
summary_medians <- all_data %>%
  group_by(group) %>%
  summarize(
    median_n_eff = median(n_eff, na.rm = TRUE),
    median_Rhat = median(Rhat, na.rm = TRUE),
    .groups = "drop"
  )



#### 07/23/2024: Pulling in newest results ####
load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_07232024.rdata')

df_predict <- NULL 
district_df <- NULL
for(i in 1:length(results_list)){
  # facility df
  tt <- results_list[[i]]$df
  max_date <- max(tt$date)
  print(max_date)
  df_predict <- rbind(df_predict,
                      tt %>% filter(date == max_date))
  
  # district df
  dist_list <- lapply(results_list[[i]]$res_list, function(xx) xx$district_df %>% filter(date == max_date))
  district_df <- Reduce(function(x, y) merge(x, y, by = c('district', 'date')), dist_list) %>%
    rbind(district_df)
}

district_df$y_pred_WF_negbin <- district_df$y_pred_WF_negbin_0.5
district_df$y_pred_freqGLMepi_negbin <- district_df$y_pred_freqGLMepi_negbin_0.5
district_df$y_CAR_phifit_negbin <- district_df$y_CAR_phifit_negbin_0.5

tmp = df_predict
tmp$y_rollup <- ifelse(is.na(tmp$y), tmp$y_pred_WF_negbin, tmp$y)
district_df <- tmp %>%
  group_by(district, date) %>%
  summarize(y = sum(y_rollup)) %>%
  merge(district_df, by = c('district', 'date'))

df_predict <- generate_outbreak_columns(df_predict, 
                                        models = c('y_pred_WF_negbin','y_pred_freqGLMepi_negbin','y_CAR_phifit_negbin'), 
                                        rename_vec = c('WF', 'freqGLM', 'CAR'))

df_o <- df_predict %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            CAR = mean(outbreak_CAR, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, freqGLM, CAR), names_to='method', values_to = 'proportion_outbreak')
df_o$method = factor(df_o$method, levels = c('WF','freqGLM','CAR'))

district_df <- generate_outbreak_columns(district_df, 
                                        models = c('y_pred_WF_negbin','y_pred_freqGLMepi_negbin','y_CAR_phifit_negbin'), 
                                        rename_vec = c('WF', 'freqGLM', 'CAR'))
district_df_o <- district_df %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            CAR = mean(outbreak_CAR, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, freqGLM, CAR), names_to='method', values_to = 'proportion_outbreak')
district_df_o$method = factor(district_df_o$method, 
                              levels = c('WF','freqGLM','CAR'))

# save(df_predict, district_df, df_o, district_df_o, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_PROCESSED_07232024.RData')

#
#### 07/23/2024: Plotting district and facility fits from new newest results ####
load('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_PROCESSED_07232024.RData')

dist_n <- df_predict %>% group_by(district) %>%
  summarize(n = length(unique(facility)))
dist_n$new_name = paste0(dist_n$district, '(', dist_n$n, ')')
district_df$district_nn <- dist_n$new_name[match(district_df$district, dist_n$district)]
district_rename <- district_df %>%
  mutate(facility = district_nn)

p1 <- plot_facility_fits(district_rename, methods = c('y_pred_WF_negbin', 'y_pred_freqGLMepi_negbin', 'y_CAR_phifit_negbin'), PI = F, imp_names = c('WF','freqGLM', 'CAR'), upper_lim = T, color_vec = c('orange3', 'forestgreen', 'blue'))

ggsave(plot = p1, filename = 'figures/Maryland_analysis_district_fits_NB_07232024.pdf', height = 5, width = 10)

p1 <- ggplot(df_o %>% filter(date >= '2020-01-01', 
                             method %in% c('WF','freqGLM','CAR'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(labels = c('WF','freqGLM','CAR'), values = c('orange3', 'forestgreen', 'blue')) +
  ylab('proportion') + 
  ggtitle('facility outbreaks') +  
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0,1)) + 
  scale_x_date(date_breaks="1 months", date_labels="%b", minor_breaks = NULL) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# legend <- get_plot_component(p1, 'guide-box', return_all = T)[[3]]
# p1 <- p1 + theme(legend.position = 'none')

p2 <- ggplot(district_df_o %>% filter(date >= '2020-01-01', 
                                      method %in% c('WF','freqGLM','CAR'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(labels = c('WF','freqGLM','CAR'), values = c('orange3', 'forestgreen', 'blue')) +
  #ylim(c(0,1)) + 
  ylab('proportion') + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0,1)) + 
  scale_x_date(date_breaks="1 months", date_labels="%b", minor_breaks = NULL) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('district outbreaks')

legend = get_plot_component(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=12)),
                            'guide-box',
                            return_all = T)

final_plot <- plot_grid(plot_grid(p1 + theme(legend.position = 'none'), p2 + theme(legend.position = 'none')),
                        legend, ncol = 1, rel_heights = c(5,1))

# ggsave(plot = final_plot, filename = 'figures/Maryland_analysis_proportion_outbreaks_NB_07232024.pdf', height = 3, width = 7)

#
#### 07/03/2024: Pulling in results with district ####
load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_07022024.rdata')

df_predict <- NULL 
district_df <- NULL
for(i in 1:length(results_list)){
  # facility df
  tt <- results_list[[i]]$df
  max_date <- max(tt$date)
  print(max_date)
  df_predict <- rbind(df_predict,
                      tt %>% filter(date == max_date))
  
  # district df
  res_tmp <- 
  dist_list <- lapply(results_list[[i]]$res_list, function(xx) xx$district_df %>% filter(date == max_date))
  district_df <- Reduce(function(x, y) merge(x, y, by = c('district', 'date')), dist_list) %>%
    rbind(district_df)
}

district_df$y_pred_WF <- district_df$y_pred_WF_0.5
district_df$y_pred_WF_negbin <- district_df$y_pred_WF_negbin_0.5
district_df$y_pred_freqGLMepi <- district_df$y_pred_freqGLMepi_0.5
district_df$y_pred_freqGLMepi_negbin <- district_df$y_pred_freqGLMepi_negbin_0.5
district_df$y_CAR_sample <- district_df$y_CAR_sample_0.5
district_df$y_CAR_phifit <- district_df$y_CAR_phifit_0.5

tmp = df_predict
tmp$y_rollup <- ifelse(is.na(tmp$y), tmp$y_pred_WF, tmp$y)
district_df <- tmp %>%
  group_by(district, date) %>%
  summarize(y = sum(y_rollup)) %>%
  merge(district_df, by = c('district', 'date'))

df_predict <- generate_outbreak_columns(df_predict, models = c('y_pred_WF','y_pred_WF_negbin','y_pred_freqGLMepi','y_pred_freqGLMepi_negbin','y_CAR_sample','y_CAR_phifit'), 
                                        rename_vec = c('WF', 'WF_negbin', 'freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))

df_o <- df_predict %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            WF_negbin = mean(outbreak_WF_negbin, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            freqGLM_negbin = mean(outbreak_freqGLM_negbin, na.rm = T),
            CAR_sample = mean(outbreak_CAR_sample, na.rm = T),
            CAR_phifit = mean(outbreak_CAR_phifit, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, WF_negbin, freqGLM, freqGLM_negbin, CAR_sample, CAR_phifit), names_to='method', values_to = 'proportion_outbreak')
df_o$method = factor(df_o$method, levels = c('WF','WF_negbin','freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))

district_df <- generate_outbreak_columns(district_df, models = c('y_pred_WF','y_pred_WF_negbin','y_pred_freqGLMepi','y_pred_freqGLMepi_negbin','y_CAR_sample','y_CAR_phifit'), 
                                         rename_vec = c('WF', 'WF_negbin', 'freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))

district_df_o <- district_df %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            WF_negbin = mean(outbreak_WF_negbin, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            freqGLM_negbin = mean(outbreak_freqGLM_negbin, na.rm = T),
            CAR_sample = mean(outbreak_CAR_sample, na.rm = T),
            CAR_phifit = mean(outbreak_CAR_phifit, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, WF_negbin, freqGLM, freqGLM_negbin, CAR_sample, CAR_phifit), names_to='method', values_to = 'proportion_outbreak')
district_df_o$method = factor(district_df_o$method, levels = c('WF','WF_negbin','freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))


# save(df_predict, district_df, df_o, district_df_o, file = 'C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_PROCESSED_07082024.RData')

#
#### 07/03/2024: Plotting district and facility fits from newest results ####
load('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/real_data_analysis_rolling_PROCESSED_07082024.RData')

dist_n <- df_predict %>% group_by(district) %>%
  summarize(n = length(unique(facility)))
dist_n$new_name = paste0(dist_n$district, '(', dist_n$n, ')')
district_df$district_nn <- dist_n$new_name[match(district_df$district, dist_n$district)]
district_rename <- district_df %>%
  mutate(facility = district_nn)

plot_facility_fits(district_rename, methods = c('y_pred_WF', 'y_pred_freqGLMepi', 'y_CAR_sample'), PI = F, upper_lim = T)

plot_facility_fits(district_rename, methods = c('y_CAR_sample', 'y_CAR_phifit'), PI = F, upper_lim = T)

p1 <- plot_facility_fits(district_rename, methods = c('y_pred_WF_negbin', 'y_pred_freqGLMepi_negbin', 'y_CAR_phifit'), PI = F, imp_names = c('WF','freqGLM', 'CAR'), upper_lim = T)

# ggsave(plot = p1, filename = 'figures/Maryland_analysis_district_fits_NB.pdf', height = 5, width = 10)

p1 <- ggplot(df_o %>% filter(date >= '2020-01-01', 
                       method %in% c('WF_negbin','freqGLM_negbin','CAR_phifit'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(labels = c('WF','freqGLM','CAR'), values = c('orange3', 'forestgreen', 'blue')) +
  ylab('proportion') + 
  ggtitle('facility outbreaks') +  
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0,1)) + 
  scale_x_date(date_breaks="1 months", date_labels="%b", minor_breaks = NULL) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# legend <- get_plot_component(p1, 'guide-box', return_all = T)[[3]]
# p1 <- p1 + theme(legend.position = 'none')

p2 <- ggplot(district_df_o %>% filter(date >= '2020-01-01', 
                       method %in% c('WF_negbin','freqGLM_negbin','CAR_phifit'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(labels = c('WF','freqGLM','CAR'), values = c('orange3', 'forestgreen', 'blue')) +
  #ylim(c(0,1)) + 
  ylab('proportion') + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0,1)) + 
  scale_x_date(date_breaks="1 months", date_labels="%b", minor_breaks = NULL) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('district outbreaks')

legend = get_plot_component(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=12)),
                            'guide-box',
                            return_all = T)[[3]]

final_plot <- plot_grid(plot_grid(p1 + theme(legend.position = 'none'), p2 + theme(legend.position = 'none')),
                        legend, ncol = 1, rel_heights = c(5,1))

ggsave(plot = final_plot, filename = 'figures/Maryland_analysis_proportion_outbreaks_NB.pdf', height = 3, width = 7)

#
#### (OLD) Pulling in test results ####
setwd(res_dir)
load('results/real_data_analysis_test2.rdata')

df_predict <- NULL 
for(i in 1:length(results_list)){
  tt <- results_list[[i]]
  max_date <- max(tt$date)
  print(max_date)
  df_predict <- rbind(df_predict,
                      tt %>% filter(date == max_date))
}

#

#### (OLD) Rolling up district data (ignore for now) ####
{
# make the base district values, imputing missing ones with WF. 
tmp = df
tmp$y_rollup <- ifelse(is.na(tmp$y), tmp$y_pred_WF, tmp$y)
district_df <- tmp %>%
  group_by(district, date) %>%
  summarize(y = sum(y_rollup))

# merge all the district results together
district_df <- merge(res_list[['WF']]$district_df, res_list[['freqGLM']]$district_df, by= c('district', 'date')) %>%
  merge(res_list[['freqGLM_NB']]$district_df, by= c('district', 'date')) %>%
  merge(res_list[['WF_NB']]$district_df, by = c('district', 'date')) %>%
  merge(res_list[['CAR']]$district_df, by = c('district', 'date')) %>%
  merge(district_df, by = c('district', 'date'))
district_df$y_pred_WF <- district_df$y_pred_WF_0.5
district_df$y_pred_WF_negbin <- district_df$y_pred_WF_negbin_0.5
district_df$y_pred_freqGLMepi <- district_df$y_pred_freqGLMepi_0.5
district_df$y_pred_freqGLMepi_negbin <- district_df$y_pred_freqGLMepi_negbin_0.5
district_df$y_CARstan <- district_df$y_CARstan_0.5
}

#save(df, district_df, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/Maryland_county_model_fits(TEST-200RPI)_06302024.RData')

#### (OLD) Plotting all of 2020 at once ####
### Plot 1: Plotting the facility fits
tmp  = df_predict %>% filter(district == 'E')
plot_facility_fits(tmp, methods = c('y_pred_WF', 'y_pred_WF_negbin'), PI = F, upper_lim = T)
plot_facility_fits(tmp, methods = c('y_pred_freqGLMepi', 'y_pred_freqGLMepi_negbin'), PI = F, upper_lim = T)
plot_facility_fits(tmp, methods = c('y_CAR_sample', 'y_CAR_phifit'), PI = F, upper_lim = T)

dist_n <- df %>% group_by(district) %>%
  summarize(n = length(unique(facility)))
dist_n$new_name = paste0(dist_n$district, '(', dist_n$n, ')')
district_df$district_nn <- dist_n$new_name[match(district_df$district, dist_n$district)]

district_rename <- district_df %>%
  mutate(facility = district_nn)
plot_facility_fits(district_rename, methods = c('y_pred_freqGLMepi', 'y_pred_freqGLMepi_negbin', 'y_pred_WF'), PI = F, upper_lim = T)
p1 <- plot_facility_fits(district_rename, 
                         methods = c('y_pred_WF', 'y_pred_freqGLMepi', 'y_CARstan'), 
                         color_vec = c('orange3', 'forestgreen', 'blue'), PI = F, upper_lim = T)

#ggsave(plot = p1, filename = 'figures/Maryland_analysis_district_fits.pdf', height = 5, width = 10)

test <- district_df %>% 
  filter(district == 'C', date >= '2020-01-01')


### Plot 2: Plotting the proportion outbreaks over time.
# maybe I just need to show the actual numbers? 
df_predict$outbreak_WF <- ifelse(df_predict$y > df_predict$y_pred_WF_0.975, 1, 0)
df_predict$outbreak_WF_negbin <- ifelse(df_predict$y > df_predict$y_pred_WF_negbin_0.975, 1, 0)
df_predict$outbreak_freqGLM <- ifelse(df_predict$y > df_predict$y_pred_freqGLMepi_0.975, 1, 0)
df_predict$outbreak_freqGLM_negbin <- ifelse(df_predict$y > df_predict$y_pred_freqGLMepi_negbin_0.975, 1, 0)
df_predict$outbreak_CAR_sample <- ifelse(df_predict$y > df_predict$y_CAR_sample_0.975, 1, 0)
df_predict$outbreak_CAR_phifit <- ifelse(df_predict$y > df_predict$y_CAR_phifit_0.975, 1, 0)
  
df_o <- df_predict %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            WF_negbin = mean(outbreak_WF_negbin, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            freqGLM_negbin = mean(outbreak_freqGLM_negbin, na.rm = T),
            CAR_sample = mean(outbreak_CAR_sample, na.rm = T),
            CAR_phifit = mean(outbreak_CAR_phifit, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, WF_negbin, freqGLM, freqGLM_negbin, CAR_sample, CAR_phifit), names_to='method', values_to = 'proportion_outbreak')
df_o$method = factor(df_o$method, levels = c('WF','WF_negbin','freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))


ggplot(df_o %>% filter(date >= '2020-01-01', 
                       method %in% c('CAR_sample','CAR_phifit'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(values = c('orange3', 'forestgreen', 'blue')) +
  ylim(c(0,1)) + 
  ylab('proportion') + 
  ggtitle('facility outbreaks') +  
  theme_bw() + 
  theme(legend.position = 'bottom') # legend.text=element_text(size=20)

district_df$outbreak_WF <- ifelse(district_df$y > district_df$y_pred_WF_0.975, 1, 0)
district_df$outbreak_WF_negbin <- ifelse(district_df$y > district_df$y_pred_WF_negbin_0.975, 1, 0)
district_df$outbreak_freqGLM <- ifelse(district_df$y > district_df$y_pred_freqGLMepi_0.975, 1, 0)
district_df$outbreak_freqGLM_negbin <- ifelse(district_df$y > district_df$y_pred_freqGLMepi_negbin_0.975, 1, 0)
district_df$outbreak_CAR_sample <- ifelse(district_df$y > district_df$y_CAR_sample_0.975, 1, 0)
district_df$outbreak_CAR_phifit <- ifelse(district_df$y > district_df$y_CAR_phifit_0.975, 1, 0)

district_df_o <- district_df %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            WF_negbin = mean(outbreak_WF_negbin, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            freqGLM_negbin = mean(outbreak_freqGLM_negbin, na.rm = T),
            CAR_sample = mean(outbreak_CAR_sample, na.rm = T),
            CAR_phifit = mean(outbreak_CAR_phifit, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, WF_negbin, freqGLM, freqGLM_negbin, CAR_sample, CAR_phifit), names_to='method', values_to = 'proportion_outbreak')
district_df_o$method = factor(district_df_o$method, levels = c('WF','WF_negbin','freqGLM','freqGLM_negbin','CAR_sample','CAR_phifit'))


ggplot(district_df_o %>% filter(date >= '2020-01-01',
                                method %in% c('CAR_phifit','CAR_sample'))) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(values = c('orange3', 'forestgreen', 'blue')) +
  ylim(c(0,1)) + 
  ylab('proportion') + 
  ggtitle('district outbreaks') + 
  theme_bw()
  # theme(legend.position = 'bottom') # legend.text=element_text(size=20)

legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=12)))

final_plot <- plot_grid(plot_grid(p1 + theme(legend.position = 'none'), p2 + theme(legend.position = 'none')),
          legend, ncol = 1, rel_heights = c(5,1))

ggsave(plot = final_plot, filename = 'figures/Maryland_analysis_proportion_outbreaks.pdf', height = 3, width = 7)


