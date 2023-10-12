library(MASS)
library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(doRNG)
library(doParallel)
library(rstan)

source('R/imputation_functions.R')

b0_mean = 6
b1_mean = -0.25

system.time({
  lst_WF <- simulate_data(district_sizes = c(4, 6, 10), 
                       R = 1, 
                       end_date = '2020-12-01',
                       type = 'WF',
                       b0_mean = b0_mean, 
                       b1_mean = b1_mean)
})

system.time({
  lst_CAR_331 <- simulate_data(district_sizes = c(4, 6, 10),
                       R = 1, 
                       end_date = '2020-12-01',
                       b0_mean = b0_mean, 
                       b1_mean = b1_mean,
                       type = 'CAR',
                       rho = 0.3, 
                       alpha = 0.3, 
                       tau2 = 1)
})

system.time({
  lst_CAR_33025 <- simulate_data(district_sizes = c(4, 6, 10),
                               R = 1, 
                               end_date = '2020-12-01',
                               b0_mean = b0_mean, 
                               b1_mean = b1_mean,
                               type = 'CAR',
                               rho = 0.3, 
                               alpha = 0.3, 
                               tau2 = 0.25)
})

system.time({
  # The log(1.5) comes from the adjustment in the code based off the sum of the geometric series of rho and alpha adjusting the estimates upward
  # b0_mean_adj = b0_mean - log(1.5)
  lst_freq <- simulate_data(district_sizes = c(4, 6, 10),
                       R = 1, 
                       end_date = '2020-12-01',
                       b0_mean = 5, 
                       b1_mean = b1_mean,
                       type = 'freqGLM',
                       rho = 0.3, 
                       alpha = 0.3)
})
# slow! 53s. This might actually be an issue that warrants improving, either by not 

### Plot them in comparison
for(i in 1:1){
  tmp = lst_WF$df_list[[i]] %>% select(facility, date, y) %>%
    mutate(method = 'WF')
  tmp2 = lst_CAR_331$df_list[[i]] %>% select(facility, date, y) %>%
    mutate(method = 'CAR_T1')
  tmp3 = lst_CAR_33025$df_list[[i]] %>% select(facility, date, y) %>%
    mutate(method = 'CAR_T025')
  tmp4 = lst_freq$df_list[[i]] %>% select(facility, date, y) %>%
    mutate(method = 'freqGLM')
  df = rbind(tmp, tmp2, tmp3, tmp4)
}

df %>% group_by(method) %>% summarize(sum(y))

ggplot(data = df, aes(x = date, y = y, color = method)) + 
  geom_line() + 
  facet_wrap(~facility, scales = 'free')

ggsave(filename = 'figures/DGP_comparison_10122023.png', width = 10, height = 6)
