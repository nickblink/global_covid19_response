library(tidyverse)
setwd('github_projects/global_covid19_response/')

D <- readRDS('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/liberia_cleaned_01-06-2021.rds')

D2 <- D %>%
  select(date, facility, ari = indicator_count_ari_total)
D2$na <- as.integer(is.na(D2$ari))

D3 <- D2 %>%
  group_by(date) %>%
  summarize(ari_sum = sum(ari, na.rm = T), 
            ari_avg = mean(ari, na.rm = T),
            na = sum(na),
            n = n()) %>%
  arrange(date)
D3$month = lubridate::month(D3$date)

plot(D3$date, D3$ari_avg, type = 'l')

D4 <- D3 %>%
  filter(date < '2020-01-01') %>%
  group_by(month) %>%
  summarize(ari_avg = mean(ari_avg),
            n = n())

# not really fair because this is not removing the time trend.
