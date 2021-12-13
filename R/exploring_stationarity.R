library(dplyr)

setwd("C:/Users/nickl/Documents/global_covid19_response")

D = readRDS('data/liberia_cleaned_NL.rds')
D$year <- lubridate::year(D$date) - min(lubridate::year(D$date)) + 1

tt = D %>% 
  filter(facility == 'Beh Town Clinic')

plot(tt$indicator_count_ari_total, type = 'l')

lm.fit <- lm(indicator_count_ari_total ~ year, data = tt)

plot(lm.fit$residuals, type = 'l')


uni_fac <- unique(D$facility)

par(mfrow = c(2,6))

for(i in 1:6){
  tt = D %>% 
    filter(facility == uni_fac[i])
  plot(tt$indicator_count_ari_total, type = 'l', main = uni_fac[i])
}

for(i in 1:6){
  tt = D %>% 
    filter(facility == uni_fac[i])
  lm.fit <- lm(indicator_count_ari_total ~ year, data = tt)
  
  resids = tt$indicator_count_ari_total
  resids[!is.na(resids)] = resids[!is.na(resids)] - lm.fit$fitted.values
  plot(resids, type = 'l', main = 'conditioned on year')
  
}
