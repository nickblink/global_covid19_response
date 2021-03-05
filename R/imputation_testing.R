library(cutoffR)

# load in functions and data
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source("R/model_functions.R")
source("R/model_figures.R")
data <- readRDS("data/data_example_singlecounty.rds")
# removing M! (for now). It's missing tooo much
data = data %>% 
  filter(facility != 'Facility M',
         date <  as.Date("2020-01-01"))

##### My imputation functions #####

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

cutoff_imputation <- function(df, df_spread = NULL, group = 'facility', method = 'number', N = NULL, cutoff = NULL){
  # if data is not wide yet, spread it!
  if(is.null(df_spread)){
    tmp = df %>% dplyr::select(date, UQ(group), indicator_denom)
    df_spread = tidyr::spread(tmp, UQ(group), indicator_denom)
  }

  # impute!
  df_imp = cutoff(df_spread, N = N, cutoff = cutoff, method = method)
  
  # convert data to long form for merging
  df_long = cbind(df_spread[,'date',drop = F], df_imp) %>% pivot_longer(!date, names_to = group, values_to = 'indicator_denom_pred_cutoff')

  # put back in original form
  df = merge(df, df_long, by = c('date',group))
  
  # return it!
  return(df)
}

impute_wrapper <- function(df, p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5), N = 2, cutoff_cor = NULL, cutoff_method = 'number', group = 'facility'){
  # create the res data frame
  res = NULL
  group_res = NULL

  for(i in (1:length(p_vec))){
    p = p_vec[i]
    df_miss <- simulate_imputation(df, 'indicator_denom', p = p, group = 'facility')
    # run two imputations
    df_imp <- tryCatch({
      periodic_imputation(df_miss, 'indicator_denom')
    }, error = function(e){
      print(sprintf('error for periodic with p = %s', p))
      browser()
    })
    #df_imp <- periodic_imputation(df_miss, 'indicator_denom')
    
    df_imp <- tryCatch({
      cutoff_imputation(df_imp, cutoff = cutoff_cor, N = N, method = cutoff_method)
    }, error = function(e){
      print(sprintf('error for cutoff with p = %s', p))
      tmp <- df_imp
      tmp$indicator_denom_pred_cutoff = NA
      return(tmp)
    })

    # get the eval period
    eval <- df_imp %>% filter(!is.na(indicator_denom_true))
    
    ### get the results over all districts
    for(imputation_method in c('indicator_denom_pred_harmonic', 'indicator_denom_pred_cutoff')){
      
      # create the res data frame and the comparison vectors
      y_true = eval$indicator_denom_true
      y_pred = eval %>% pull(imputation_method)
      
      # calculate metrics
      res_row = data.frame(p = p, imputation_method = imputation_method, 
                           RMSE = sqrt(mean((y_true - y_pred)^2)), 
                           MAE = mean(abs(y_true - y_pred)), 
                           MAPE = mean(abs(y_true - y_pred)/y_true), 
                           N_imp = nrow(eval))
      
      # update full results
      res = rbind(res, res_row)
    }
    
    ### get the group level results
    # get the unique groups
    uni_group = eval %>% pull(get(group)) %>% unique()
    
    for(imputation_method in c('indicator_denom_pred_harmonic', 'indicator_denom_pred_cutoff')){
      
      # cycle through each group/facility
      for(g in uni_group){
        eval2 = eval %>% filter(get(group) == g)
        y_true = eval2 %>% pull(indicator_denom_true)
        y_pred = eval2 %>% pull(imputation_method)
        
        # calculate metrics
        res_row = data.frame(p = p, imputation_method = imputation_method, group = g,
                             RMSE = sqrt(mean((y_true - y_pred)^2)), 
                             MAE = mean(abs(y_true - y_pred)), 
                             MAPE = mean(abs(y_true - y_pred)/y_true), 
                             N_imp = nrow(eval2))
        
        # update full results
        group_res = rbind(group_res, res_row)
      }
      
    }
  }
  
  # return the results
  res_lst = list(full_results = res, group_results = group_res)
  return(res_lst)
}

##### Running imputation #####
  
res = impute_wrapper(data, p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
res$full_results


### Plotting - not prettied up at all, but ok for now.
tt = res$group_results %>% filter(imputation_method == 'indicator_denom_pred_harmonic')
tt = tt %>% filter(p <= 0.5)
tt$p = factor(tt$p, levels = sort(unique(tt$p)))

ggplot(tt, aes(x = p, y = RMSE)) +
  geom_boxplot()

ggplot(tt, aes(x = p, y = MAE)) +
  geom_boxplot()

ggplot(tt, aes(x = p, y = MAPE)) +
  geom_boxplot()



# yay first pass is done

impute_wrapper(data, p_vec = c(0.05, 0.1, 0.15), cutoff_method = 'correlation', cutoff_cor = 0.2)$full_results




##### Testing cutoff and stuff ####

data(hqmr.data)
# so it looks like it is in the wide form with all columns being different groups and then one column for data. Not sure if that column is necessary though, is it? Nvm it should be

# ok just needs "date" and everything else is imputed
test = cutoff(hqmr.data)
d2 = cbind(hqmr.data[,79,drop =F], hqmr.data[,1:78])
test = cutoff(d2)

length(unique(hqmr.data$date))
head(hqmr.data$date, 13)
# so it's by month. And this is how cutoff will search for previous years.

# so I already have the date column. want to create a df
tmp = data %>% dplyr::select(date, facility, indicator_denom)
data_spread = tidyr::spread(tmp, facility, indicator_denom)

test = cutoff(data_spread)
# ok it works. Of course, not optimized right now.

cutoff_imputation(data, N = 2)
cutoff_imputation(data, cutoff = 0.5, method = 'correlation')
