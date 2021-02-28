
setwd('C:/Users/nickl/Documents/global_covid19_response/')

# load in functions and data
source("R/model_functions.R")
source("R/model_figures.R")
data <- readRDS("data/data_example_singlecounty.rds")

# exploring data
str(data)

# function to simulate data
simulate_imputation <- function(df, col, p = 0.1, group = NULL){
  set.seed(1)
  # if not worrying about the grouping of randomization
  
  # set the true value to null for the df
  df[,paste0(col, '_true')] = as.integer(NA)
  
  if(is.null(group)){
    # get the number of data points to impute
    num_impute = round(p*nrow(df)) - sum(is.na(df[,col]))
    
    if(num_impute <= 0){
      warning('skipping imputation. Already too many missing')
      return(df)
    }
    
    # sample the indices to impute
    ind_sample = sample(setdiff(1:nrow(df), which(is.na(df[,col]))), num_impute)
    
    # store the true value of the values to replace
    df[ind_sample, paste0(col, '_true')] <- df[ind_sample, col]
    
    # replace the values with NAs
    df[ind_sample, col] = NA
    
  }else{
    # doing grouping
    uni_group = df %>% pull(UQ(group)) %>% unique()
    
    # apply to each group
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      num_impute = round(p*nrow(tt)) - sum(is.na(tt[,col]))
      if(num_impute <= 0){
        print(sprintf('skipping imputation for %s', xx))
        return(tt)
      }
      ind_sample = sample(setdiff(1:nrow(tt), which(is.na(tt[,col]))), num_impute) 
      tt[ind_sample, paste0(col, '_true')] <- tt[ind_sample, col]
      tt[ind_sample, col] = NA
      return(tt)
    })
    
    # collapse results
    df <- data.table::rbindlist(tmp)
  }
  
  return(df)
}

# OG imputation method
periodic_imputation <- function(df, col, group = 'facility', family = 'quasipoisson', period = 12){
  # prep the data with the harmonic functions
  df = df %>%
    dplyr::mutate(year = year(date) - min(year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  
  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  
  if(family == 'quasipoisson'){
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=quasipoisson)
      tt[,paste0(col, '_pred_harmonic')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)

  }else{
    stop('havent coded negative binomial or other ways yet')
  }
  
}

impute_wrapper <- function(df, p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5)){
  # create the res data frame
  res = data.frame(p = p_vec, RMSE = 0, MAE = 0, MAPE = 0)

  for(i in (1:length(p_vec))){
    p = p_vec[i]
    tmp <- simulate_imputation(df, 'indicator_denom', p = p, group = 'facility')
    tmp <- periodic_imputation(tmp, 'indicator_denom')
    eval <- tmp %>% filter(!is.na(indicator_denom_true))
    
    # calculate metrics
    res[i,'RMSE'] = sqrt(mean((eval$indicator_denom_true - eval$indicator_denom_pred_harmonic)^2))
    res[i,'MAE'] = mean(abs(eval$indicator_denom_true - eval$indicator_denom_pred_harmonic))
    res[i,'MAPE'] = mean(abs(eval$indicator_denom_true - eval$indicator_denom_pred_harmonic)/eval$indicator_denom_true)
  }
  
  # return the results
  return(res)
}

impute_wrapper(data)
# yay first pass is done


