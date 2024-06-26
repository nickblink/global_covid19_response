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

### Data prep
{
Dfull <- readRDS(sprintf('%s/data/liberia_cleaned_01-06-2021.rds', res_dir))

Dfull %>% group_by(county) %>%
  summarize(n = length(unique(facility)))

# most recent data on github.
D <- Dfull %>%
  filter(county == 'Maryland') %>%
  select(date, district, facility, ari = indicator_count_ari_total, indicator_denom) %>%
  add_periodic_cov()

D %>% 
  filter(date >= '2020-01-01') %>%
  group_by(date) %>% 
  summarize(n_miss = length(unique(is.na(ari))))

D %>% 
  filter(date < '2020-01-01') %>%
  summarize(prop_miss = mean(is.na(ari)))

# get the dates
dates = unique(D$date)
eval_dates = dates[dates >= '2020-01-01']

# anonymize and make name matching fit the code.
dist_n <- D %>% group_by(district) %>%
  summarize(n = length(unique(facility))) %>% 
  arrange(n)
uni_district = unique(dist_n$district)
matid = match(D$district, uni_district)
D$district = toupper(letters[matid])

# rename the facilities
df = NULL
for(d in unique(D$district)){
  tmp = D %>% filter(district == d)
  uni_fac = unique(tmp$facility)
  matid = match(tmp$facility, uni_fac)
  tmp$facility = paste0(d, matid)
  df = rbind(df, tmp)
}

df$y <- df$ari
}
res_list <- list()
R_PI = 2 # 200
burnin = 100 # 1000
nsample = 200 # 2000

#### Rolling surveillance ####
# get dates
all_dates = sort(unique(df$date))
eval_dates = all_dates[all_dates >= '2020-01-01']

# full res_list and prediction
full_res_list <- list()
df_predict <- NULL

for(d in eval_dates){
  res_list <- list()
  # get the dates
  # dates <- all_dates[all_dates >= (eval %m-% months(48)) & all_dates <= eval]
  ind <- which(all_dates == d)
  dates <- all_dates[(ind - 48):ind]
  train_end <- all_dates[ind-1]
  
  # get the data frame
  df_roll <- df %>% 
    filter(date %in% dates)
  
  # run WF model
  res_list[['WF']] <- WF_CCA(df_roll, col = "y", family = 'poisson', R_PI = R_PI, train_end_date = train_end)
  df_roll <- res_list[['WF']]$df
  
  # run WF NB model
  res_list[['WF_NB']] <- WF_CCA(df_roll, col = "y", family = 'negbin', R_PI = R_PI, train_end_date = train_end)
  df_roll <- res_list[['WF_NB']]$df
  
  # run freqGLM
  system.time({
    res_list[['freqGLM']] <- freqGLMepi_CCA(df_roll, R_PI = R_PI, verbose = F, train_end_date = train_end)
  }) # 20m
  df_roll <- res_list[['freqGLM']]$df
  
  # run freqGLM NB
  system.time({
    res_list[['freqGLM_NB']] <- freqGLMepi_CCA(df_roll, R_PI = R_PI, verbose = F, family = 'negbin', train_end_date = train_end)
  }) # 20m
  df_roll <- res_list[['freqGLM_NB']]$df
  
  # run CAR
  system.time({
    res_list[['CAR']] <- CARBayes_wrapper(df_roll, burnin = burnin, n.sample = n.sample, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', train_end_date = train_end)
  }) # 143s
  df_roll <- res_list[['CAR']]$df
  df_roll$y_CARstan <- df_roll$y_CARstan_0.5
  colnames(df_roll) <- gsub('y_CARstan', 'y_CAR_sample', colnames(df_roll))
  
  system.time({
  res_list[['CAR_phifit']] <- CARBayes_wrapper(df_roll, burnin = burnin, n.sample = n.sample, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan', train_end_date = train_end, use_fitted_phi = T)
  })
  df_roll <- res_list[['CAR_phifit']]$df
  df_roll$y_CARstan <- df_roll$y_CARstan_0.5
  colnames(df_roll) <- gsub('y_CARstan', 'y_CAR_phifit', colnames(df_roll))
  
  # update results
  df_predict <- rbind(df_predict, 
                      df_roll[df_roll$date == d,])
  full_res_list[[d]] <- df_roll
}
  
#save(df_roll, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/Maryland_county_model_fits_rolling_06302024.RData')

# tmp  = df_roll %>% filter(district == 'A')
# plot_facility_fits(tmp, methods = c('y_pred_WF','y_CAR_sample','y_CAR_phifit'),
#                    include_legend = T)

#
#### Predicting all of 2020 at once ####
### Running the models
{
# run WF model
res_list[['WF']] <- WF_CCA(df, col = "y", family = 'poisson', R_PI = R_PI)
df <- res_list[['WF']]$df

# run WF NB model
res_list[['WF_NB']] <- WF_CCA(df, col = "y", family = 'negbin', R_PI = R_PI)
df <- res_list[['WF_NB']]$df

# run freqGLM
system.time({
  res_list[['freqGLM']] <- freqGLMepi_CCA(df, R_PI = R_PI, verbose = F)
}) # 20m
df <- res_list[['freqGLM']]$df

# run freqGLM NB
system.time({
  res_list[['freqGLM_NB']] <- freqGLMepi_CCA(df, R_PI = R_PI, verbose = F, family = 'negbin')
}) # 20m
df <- res_list[['freqGLM_NB']]$df

# run CAR
res_list[['CAR']] <- CARBayes_wrapper(df, burnin = 1000, n.sample = 2000, prediction_sample = T, predict_start_date = '2016-01-01', MCMC_sampler = 'stan')
df <- res_list[['CAR']]$df
df$y_CARstan <- df$y_CARstan_0.5
# 430s
}

### Rolling up district data
{
# make the base district values, imputing missing ones with WF. 
tmp = df
tmp$y_rollup <- ifelse(is.na(tmp$y), tmp$y_pred_WF, tmp$y)
district_df <- tmp %>%
  group_by(district, date) %>%
  summarize(y = sum(y_rollup))

# tt = do.call('merge', c(list(district_df,
#                   res_list[['WF']]$district_df, 
#                   res_list[['freqGLM']]$district_df,
#                   res_list[['CAR']]$district_df,
#                   res_list[['WF_NB']]$district_df), by= c('district', 'date')))

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

#### Plotting all of 2020 at once ####
### Plot 1: Plotting the facility fits
tmp  = df %>% filter(district == 'A')
plot_facility_fits(tmp, methods = c('y_pred_WF','y_pred_freqGLMepi', 'y_CARstan', 'y_pred_WF_negbin', 'y_pred_freqGLMepi_negbin'))
# plot_facility_fits(tmp, methods = c('y_CARstan', 'y_pred_WF_negbin', 'y_pred_freqGLMepi_negbin'))
# plot_facility_fits(tmp, methods = c('y_CARstan', 'y_pred_freqGLMepi', 'y_pred_freqGLMepi_negbin'))

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
df$outbreak_WF <- ifelse(df$y > df$y_pred_WF_0.975, 1, 0)
df$outbreak_freqGLM <- ifelse(df$y > df$y_pred_freqGLMepi_0.975, 1, 0)
df$outbreak_CAR <- ifelse(df$y > df$y_CARstan_0.975, 1, 0)
  
df_o <- df %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            CAR = mean(outbreak_CAR, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, freqGLM, CAR), names_to='method', values_to = 'proportion_outbreak')
df_o$method = factor(df_o$method, levels = c('WF','freqGLM','CAR'))


p1 <- ggplot(df_o %>% filter(date >= '2020-01-01')) +
  geom_line(aes(x = date, y = proportion_outbreak, color = method)) + 
  scale_color_manual(values = c('orange3', 'forestgreen', 'blue')) +
  ylim(c(0,1)) + 
  ylab('proportion') + 
  ggtitle('facility outbreaks') +  
  theme_bw() + 
  theme(legend.position = 'bottom') # legend.text=element_text(size=20)


district_df$outbreak_WF <- ifelse(district_df$y > district_df$y_pred_WF_0.975, 1, 0)
district_df$outbreak_freqGLM <- ifelse(district_df$y > district_df$y_pred_freqGLMepi_0.975, 1, 0)
district_df$outbreak_CAR <- ifelse(district_df$y > district_df$y_CARstan_0.975, 1, 0)

district_df_o <- district_df %>%
  group_by(date) %>%
  summarize(WF = mean(outbreak_WF, na.rm = T),
            freqGLM = mean(outbreak_freqGLM, na.rm = T),
            CAR = mean(outbreak_CAR, na.rm = T)) %>%
  tidyr::pivot_longer(c(WF, freqGLM, CAR), names_to='method', values_to = 'proportion_outbreak')
district_df_o$method = factor(district_df_o$method, levels = c('WF','freqGLM','CAR'))


p2 <- ggplot(district_df_o %>% filter(date >= '2020-01-01')) +
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


