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

if(file.exists('C:/Users/Admin-Dell')){
  res_dir = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}else{
  res_dir = "C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance"
}

# most recent data on github.
D = readRDS(sprintf('%s/data/liberia_cleaned_01-06-2021.rds', res_dir)) %>%
  filter(county == 'Maryland') %>%
  select(date, district, facility, ari = indicator_count_ari_total, indicator_denom)

# get the dates
dates = unique(D$date)
eval_dates = dates[dates >= '2020-01-01']

# anonymize and make name matching fit the code.
uni_district = unique(D$district)
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

# run WF model (this only needs to be fit once).

# run freqGLM (this needs to be fit for each date because of rolling baseline).

# run CAR (this needs to be fit for each date because of rolling baseline).

for(d in eval_dates){
}


