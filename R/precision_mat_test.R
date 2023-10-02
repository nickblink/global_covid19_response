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

file_name = 'C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/mar01_nost_beta43_n025_2023_02_26/sim_results_p0.1_mar_1(50).RData'

load(file_name)

df = imputed_list[[1]]$df_miss

rho = 0.5
Q = make_precision_mat(df, rho = rho)

det(Q)

W2 <- make_district_W2_matrix(df)

lambdas = eigen(W2 - diag(1, nrow(W2)))$values

res = 1
for(i in 1:length(lambdas)){
  res = res*(1 + rho*lambdas[i])
}


