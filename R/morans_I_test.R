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

setwd(res_dir)

get_morans_I <- function(df, out_col = 'y', start_date = NULL, end_date = NULL, avg_across_time = T){
  # get the dates
  dates <- unique(df$date)
  if(!is.null(start_date)){dates <- dates[dates >= start_date]}
  if(!is.null(end_date)){dates <- dates[dates <= end_date]}

  # get data adjacency matrix
  #W <- make_district_adjacency(df)
  
  I_vec<- sapply(dates, function(dd){
    # filter by date and non-missingness
    tmp <- df %>% filter(date == dd)
    tmp <- tmp[!is.na(tmp[,out_col]),]
    
    # get adjacency
    W2 <- make_district_adjacency(tmp)
    print(dd)
    
    browser()
    # compute Moran's I.
    y_exp = mean(tmp[,out_col])
    tmp$residual <- tmp[,out_col] - y_exp
    numerator <- sum(tmp$residual * W2%*%tmp$residual)/sum(W2)
    denominator <- mean((tmp[,out_col] - y_exp)^2)
    I_tmp <- numerator/denominator
    return(I_tmp)
  })
  
  if(avg_across_time){
    return(mean(I_vec))
  }else{
    return(I_vec)
  }
}

load('mcar01_car0303025_beta_06_beta1n025_negbin_2024_07_11/sim_results_p0.1_mcar_1(50).RData')

df <- imputed_list[[1]]$df_miss %>%
  select(date, facility, district, year, month, y, y_true)

get_morans_I(df)
get_morans_I(df, out_col= 'y_true')


load('mcar05_car0303025_beta_06_beta1n025_negbin_2024_07_11/sim_results_p0.5_mcar_1(50).RData')

df <- imputed_list[[1]]$df_miss %>%
  select(date, facility, district, year, month, y, y_true)

get_morans_I(df)
get_morans_I(df, out_col= 'y_true')
