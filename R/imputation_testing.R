library(cutoffR)
library(rstanarm)

# load in functions and data
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source("R/model_functions.R")
source("R/model_figures.R")
data <- readRDS("data/data_example_singlecounty.rds")
# data from county "A"
data = data %>%
  filter(facility != 'Facility M',
         date <  as.Date("2020-01-01"))
# (removed M because it's missing too much)

# county-level data
data2 <- readRDS('data/ari_total_county_revisions.rds') %>%
  filter(date <  as.Date("2020-01-01"))
mean(is.na(data)) # ok no missingness

# data for all facilities
D = readRDS('data/liberia_cleaned_NL.rds')

#
##### My imputation functions #####

# function to add periodic terms
add_periodic_cov <- function(df, period = 12){
  df = df %>%
    dplyr::mutate(year = year(date) - min(year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  return(df)
}

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
  df <- add_periodic_cov(df, period = period)

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

  }else if(family %in% c('NB','negative binomial','neg_bin')){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- MASS::glm.nb(formula_col, data = tt)
      tt[,paste0(col, '_pred_harmonic')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)
  }
  
}

# Bayes imputation method
bayes_periodic_imputation <- function(df, df_OG = NULL, col, group = 'facility', family = 'NB', period = 12, iterations = 500, harmonic_priors = F){
  # prep the data with the harmonic functions
  df <- add_periodic_cov(df, period = period) %>% as.data.frame()
  
  # make the target an integer for the Bayesian approach to work.
  df[,col] = as.integer(df[,col])

  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  
  if(harmonic_priors){
    
    if(!is.null(df_OG)){
      # setting the prior to use the un-imputed matrix models
      df_prior = add_periodic_cov(df_OG, period = period) %>% as.data.frame()
      pred_col_name = paste0(col, '_pred_bayes_harmonic_ideal')
    }else{
      # setting the prior to use the imputed matrix models
      df_prior = df
      pred_col_name = paste0(col, '_pred_bayes_harmonic_miss')
    }
    
    # get the max value
    max_val = max(df_prior[,col])
    
    # run the individual model for each group.
    tmp <- lapply(uni_group, function(xx) {
      tt <- df_prior %>% filter(get(group) == xx)
      
      # run the model
      mod_col <- MASS::glm.nb(formula_col, data = tt)
      
      # return the coefficients
      res_df = data.frame(coef(mod_col))
      colnames(res_df)[1] = xx
      
      return(data.frame(t(res_df)))
      
      # scale to get an idea of the intercept (ignoring)
      # tt[,col] = (max_val/max(tt[,col]))*tt[,col]
      
    })
    
    # get a matrix of the betas, including the intercepts.
    betas = do.call('rbind',tmp)
  }
  
  if(family == 'quasipoisson'){
    stop('havent coded Bayes quasipoisson')

    
  }else if(family %in% c('NB','negative binomial','neg_bin')){
   

      # train the model
      if(harmonic_priors){
        tmp <- lapply(uni_group, function(xx) {
          # get the data frame
          tt <- df %>% filter(get(group) == xx)
          
          # get the distribution of the other betas
          betas2 = betas[-which(rownames(betas) == xx),]
          beta_prior = normal(location = colMeans(betas2[,-1], na.rm = T), scale = apply(betas2[,-1], 2, function(xx) {sd(xx, na.rm = T)}))
          intercept_prior = normal(location = mean(betas2[,1]), scale = sd(betas2[,1]))
          
          # run the model with priors
          tryCatch({
            mod_bayes <- stan_glm(formula_col, data = tt, family = neg_binomial_2, refresh = 0, iter = iterations, prior = beta_prior, prior_intercept = intercept_prior)
          }, error = function(e){
            browser()
          })
          
          
          # need to fill in NAs, even if they aren't used
          tmp = tt; tmp[is.na(tmp)] = 1e6
          
          # predict the median output from a posterior draw
          tt[, pred_col_name] = apply(posterior_predict(mod_bayes, tmp, type = 'response'), 2, median)
          
          return(tt)
        })
        
        # prior_summary(mod_bayes)
        # so the prior is adjusted based on scale of values. Ok ok.
        
      }else{
        tmp <- lapply(uni_group, function(xx) {
          # get the data frame
          tt <- df %>% filter(get(group) == xx)
          
          # run the models
          mod_bayes <- stan_glm(formula_col, data = tt, family = neg_binomial_2, refresh = 0, iter = iterations)
          
          # need to fill in NAs, even if they aren't used
          tmp = tt; tmp[is.na(tmp)] = 1e6
          
          # predict the median output from a posterior draw
          tt[, paste0(col, '_pred_bayes_harmonic')] = apply(posterior_predict(mod_bayes, tmp, type = 'response'), 2, median)
          
          return(tt)
      })

    }
    #browser()
    df <- data.table::rbindlist(tmp)
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

impute_wrapper <- function(df, col = 'indicator_denom', p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5), N = 2, cutoff_cor = NULL, harmonic_family = 'NB', bayes_harmonic_family = 'NB', cutoff_method = 'number', group = 'facility', bayes_iterations = 500, bayes_beta_prior = T){
  
  if(col != 'indicator_denom'){
    warning('cutoff not coded to do imputation on other columns')
  }
  
  # create the res data frame
  res = NULL
  group_res = NULL

  for(i in (1:length(p_vec))){
    # simulate the missing data
    p = p_vec[i]
    print(p)
    df_miss <- simulate_imputation(df, col, p = p, group = group)
    
    # run two imputations
    df_imp <- tryCatch({
      periodic_imputation(df_miss, col = col, family = harmonic_family, group = group)
    }, error = function(e){
      print(sprintf('error for periodic with p = %s', p))
      browser()
    })
    
    if(!is.na(bayes_harmonic_family)){
      df_imp <- tryCatch({
        bayes_periodic_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group)
      }, error = function(e){
        print(sprintf('error for bayes periodic with p = %s', p))
        browser()
      })
    }
    
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_periodic_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
      }, error = function(e){
        print(sprintf('error for bayes periodic beta prior with p = %s', p))
        browser()
      })
    }
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_periodic_imputation(df_imp, df_OG = df, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
      }, error = function(e){
        print(sprintf('error for bayes periodic beta prior with p = %s', p))
        browser()
      })
    }
    
    if(!is.na(cutoff_method)){
      df_imp <- tryCatch({
        cutoff_imputation(df_imp, cutoff = cutoff_cor, N = N, method = cutoff_method, group = group)
      }, error = function(e){
        print(sprintf('error for cutoff with p = %s', p))
        tmp <- df_imp
        tmp$indicator_denom_pred_cutoff = NA
        return(tmp)
      })
    }

    # get the "true" column names
    true_col = paste0(col,'_true')
    
    # get the eval period
    eval <- df_imp %>% filter(!is.na(get(true_col)))
    
    ### get the results over all districts
    for(imputation_method in grep(paste0(col,'_pred'), colnames(eval), value = T)){
      
      # create the res data frame and the comparison vectors
      y_true = eval[,get(true_col)]
      y_pred = eval %>% pull(imputation_method)
      
      #if(imputation_method == 'indicator_denom_pred_bayes_harmonic'){browser()}
      
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
    
    for(imputation_method in grep(paste0(col,'_pred'), colnames(eval), value = T)){
      
      # cycle through each group/facility
      for(g in uni_group){
        eval2 = eval %>% filter(get(group) == g)
        y_true = eval2 %>% pull(UQ(paste0(col,'_true')))
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

plot_results <- function(res, subset = NULL){
  tt = res$group_results %>%
    mutate(p = as.factor(p),
           imputation_method = gsub('observed_count_pred_', '', imputation_method)) %>%
    mutate(imputation_method = factor(imputation_method, levels = c("harmonic", "bayes_harmonic", "bayes_harmonic_miss", "bayes_harmonic_ideal")))
  
  # subset data if that's what's asked for
  if(!is.null(subset)){
    tt = tt %>% filter(imputation_method %in% subset)
  }
  
  # RMSE plot
  p1 <- ggplot(tt, aes(x = p, y = RMSE, fill = imputation_method)) +
    geom_boxplot() + 
    ggtitle('RMSE of imputation predictions') + 
    ylim(c(0,1500))
  
  # MAPE plot
  p2 <- ggplot(tt, aes(x = p, y = MAPE, fill = imputation_method)) +
    geom_boxplot() + 
    ggtitle('Mean Absolute Prediction Error of imputation predictions') + 
    ylim(c(0,1))
  
  
  # final plot
  plot_final = cowplot::plot_grid(p1, p2, ncol = 1)
  return(plot_final)
}

##### Comparing my OG imputation to github code #####
tmp <- data %>% filter(facility == 'Facility D')
tt <- periodic_imputation(tmp, 'indicator_count_ari_total') %>%
  dplyr::select(date, indicator_count_ari_total, indicator_count_ari_total_pred_harmonic)
ss <- fit.site.specific.denom.pi(tmp, 'Facility D',extrapolation_date = as.Date("2020-01-01"), 'indicator_count_ari_total', site_var = 'facility', date_var = 'date', counts_only = T) %>%
  dplyr::select(date, est_raw_counts)

test = merge(tt, ss, by = 'date')

# ok good it works. This doesn't account for denominator, but it works.

##### Running imputation #####

res = impute_wrapper(data2, col = 'observed_count', p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), cutoff_method = NA, bayes_iterations = 2000, group = 'county')
plot_results(res)

HERE ON 3/22 at 12:25pm

#save(res, file = 'results/allcounties_bayes_imputation.RData')

# data_miss = simulate_imputation(data2, 'observed_count', p = 0.9, group = 'county')
# 
# tt = bayes_periodic_imputation(data_miss, df_OG = data2, col = 'observed_count', group = 'county', family = 'NB', period = 12, iterations = 500, harmonic_priors = T)
# 
# res = impute_wrapper(data2, col = 'observed_count', p_vec = c(0.8), cutoff_method = NA, bayes_iterations = 200, group = 'county')
##### Modeling the paper #####
# Declare this for all functions
extrapolation_date <- "2020-01-01"

# filter to only look at one county
# Bomi has 0.255 ARI missingness. Not the highest but reasonably high
D2 = D %>% filter(county == 'Bomi')

# check the distribution of facility missingness
#   any over 0.2? 
facility_miss = D2 %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))
# 13/27 have > 20% missing ARI. And 4 (?) have > 20% missing denom (not sure if the code does anything about that, though)

# only keep facilities with < 80% missing
facility_list = facility_miss %>% filter(ari_miss < 0.8) %>% distinct(facility) %>% pull()

D2 = D2 %>%
  filter(facility %in% facility_list)

# full model fit (remember this is excluding > 0.2 missing)
county.fit = fit.aggregate.pi.boot(D2, 
                      indicator_var = "indicator_count_ari_total",
                      denom_var = "indicator_denom",
                      site_var = "facility",
                      date_var = "date",
                      counts_only=FALSE,
                      R=250)

# district fits
uni_districts = unique(D2$district)
district.fits <- lapply(uni_districts, function(xx){
                        fit.cluster.pi(data = D2 %>% filter(district == xx), 
                                       indicator_var = "indicator_count_ari_total",
                                       denom_var = "indicator_denom",
                                       site_var = "facility",
                                       date_var = "date",
                                       counts_only=FALSE,
                                       R=250) %>%
    mutate(site = xx)
  })

# clinic fits
facility.results <- do.call(rbind, lapply(facility_list,function(x){
  fit.site.specific.denom.pi(data=D2,
                             site_name=x,
                             extrapolation_date=extrapolation_date,
                             indicator_var="indicator_count_ari_total",
                             denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                             site_var="facility",
                             date_var="date",
                             R=500)
})
) 

# individual facility results
plot_facet(facility.results,type="count")

# full county results
plot_site(county.fit, "count", title="Bomi Aggregated Results")

###  
# (1) Impute all the values. Now Bayesian imputation! <- version 1
D2 <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility')

# updating the missing values
D2$ari_bayes_imp = D2$indicator_count_ari_total
D2$ari_bayes_imp[is.na(D2$ari_bayes_imp)] = D2$indicator_count_ari_total_pred_bayes_harmonic[is.na(D2$ari_bayes_imp)]

# (1) Impute all the values. Now Bayesian imputation! <- version 2
D2 <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility', harmonic_priors = T)

# updating the missing values
D2$ari_bayes_imp2 = D2$indicator_count_ari_total
D2$ari_bayes_imp2[is.na(D2$ari_bayes_imp2)] = D2$indicator_count_ari_total_pred_bayes_harmonic_miss[is.na(D2$ari_bayes_imp2)]

# (2) Run the same models as before 
facility.results.bayes <- do.call(rbind, lapply(facility_list,function(x){
  fit.site.specific.denom.pi(data=D2,
                             site_name=x,
                             extrapolation_date=extrapolation_date,
                             indicator_var="ari_bayes_imp",
                             denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                             site_var="facility",
                             date_var="date",
                             R=500,
                             counts_only = T)
})
) 
plot_facet(facility.results.bayes, type="count")


facility.results.bayes2 <- do.call(rbind, lapply(facility_list,function(x){
  fit.site.specific.denom.pi(data=D2,
                             site_name=x,
                             extrapolation_date=extrapolation_date,
                             indicator_var="ari_bayes_imp2",
                             denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                             site_var="facility",
                             date_var="date",
                             R=500)
})
) 

plot_facet(facility.results.bayes2, type="count")

# ok some issues here - to return later
tmp = D2 %>% filter(facility == 'Gonjeh Clinic')

# my own plot
plot_list = list()
iter = 0
for(f in facility_list){
  iter = iter + 1
  tmp = D2 %>% filter(facility == f)
  tmp$indicator_count_ari_total_pred_bayes_harmonic_miss[!is.na(tmp$indicator_count_ari_total)] = NA
  p1 <- ggplot(tmp, aes(x = date, y = indicator_count_ari_total)) + 
    geom_line() + 
    geom_line(aes(x = date, y = indicator_count_ari_total_pred_bayes_harmonic_miss), col = 'blue')
  plot_list[[iter]] = p1
}

cowplot::plot_grid(plotlist = plot_list)

# work on filling in the gaps of this visualization. For later.

# note the original bayesian method will probably give really similar betas but with smaller confidence intervals, which is sort of a lie.


