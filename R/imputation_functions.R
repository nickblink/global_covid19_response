### The functions used for imputation in Liberia
library(MASS)
#library(CARBayesST)
library(Matrix)
library(dplyr)
library(lubridate)
library(ggplot2)
library(cowplot)

##### Helper Functions #####
# function to add periodic terms
add_periodic_cov <- function(df, period = 12){
  df = df %>%
    dplyr::mutate(year = lubridate::year(date) - min(lubridate::year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  return(df)
}

# add autoregressive terms to the data
add_autoregressive <- function(df, target_col = 'y', num.terms = 1){
  
  if(!(target_col %in% colnames(df))){
    stop('need a target column to autoregress')
  }
  
  tmp <- lapply(unique(df$facility), function(xx) {
    tt <- df %>% filter(facility == xx) %>% arrange(date)
    #tt[,sprintf('%s.AR1', target_col)] = c(NA, tt[1:(nrow(tt) - 1), target_col, drop = T])
    tt[,'y.AR1'] = c(NA, tt[1:(nrow(tt) - 1), target_col, drop = T])
    return(tt)
  })
  
  df <- do.call('rbind',tmp)
  return(df)
}

# compute the sum of the values at all the neighbors to a given facility
add_neighbors <- function(df, target_col = 'y', lag = 0, scale_by_num_neighbors = F){
  if(lag > 0){
    stop('havent implemented lags for neighbor calculation')
  }
  
  # remove the neighbor column
  if('y.neighbors' %in% colnames(df)){
    df$y.neighbors = NULL
  }
  
  # get the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  W = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # scale the neighbor sum by the number of neighbors
  if(scale_by_num_neighbors){
    W = apply(W, 2, function(xx) xx/sum(xx))
  }
  
  # get the counts for each facility 
  y.counts <- df %>% 
    dplyr::select(date, facility, UQ(target_col)) %>%
    tidyr::spread(facility,get(target_col)) %>% 
    arrange(date)
  

  if(!identical(colnames(W), colnames(y.counts)[-1])){
    # if the colnames do match but are in the wrong order, reorder them
    if(length(setdiff(colnames(W), colnames(y.counts)[-1])) == 0 &
       length(setdiff(colnames(y.counts)[-1], colnames(W))) == 0){
      
      matid = match(colnames(y.counts)[-1], colnames(W))
      y.counts = y.counts[,c(1, matid + 1)]
      
      # for error checking
      if(!identical(colnames(W), colnames(y.counts)[-1])){browser()}
    }else{
      stop('Adjacency and y matrix column names dont match')
    }
    ind = which(colnames(W) != colnames(y.counts)[-1])
  }
  
  df = cbind(y.counts[,'date',drop = F], as.data.frame(as.matrix(y.counts[,-1])%*%W)) %>% 
    tidyr::gather(facility, y.neighbors, -date) %>%
    merge(df, by = c('date','facility'))
  
  return(df)
}

date_facility_check <- function(imputed_list){
  date_fac_mat = imputed_list[[1]][, c("date", "facility")]
  date_fac = paste0(date_fac_mat[,1], '--', date_fac_mat[,2])
  
  for(i in 2:length(imputed_list)){
    date_fac_mat_new = imputed_list[[i]][, c("date", "facility")]
    date_fac_new = paste0(date_fac_mat_new[,1], '--', date_fac_mat[,2])
    if(!identical(date_fac_new, date_fac)){
      stop('mismatch in dates and facilities in the imputed list run')
    }
  }
}

fit_WF_model <- function(data, outcome = 'indicator_count_ari_total', facilities = NULL, family = 'nb'){
  # use all facilities if not provided
  if(is.null(facilities)){
    facilities <- unique(data$facility)
  }
  
  data <- add_periodic_cov(data)
  data$y <- as.data.frame(data)[,outcome]

  formula_col = as.formula(sprintf("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"))
  
  print(sprintf('fitting models for %s family', family))
  if(family %in% c('nb','NB','negative binomial')){
    tmp <- lapply(facilities, function(xx) {
      tt <- data %>% filter(facility == xx)
      
      # run the model
      mod_col <- tryCatch({MASS::glm.nb(formula_col, data = tt)},
                          error = function(e){
                            NA
                          })
      
      return(mod_col)
    })
  }else{
    stop('havent coded for other families')
  }
  
  # combine the results
  df = do.call('rbind', lapply(tmp, function(xx) {
    if(is.na(xx)){
      return(rep(NA, 8))
    }else{
      if(!xx$converged){
        return(rep(NA, 8))
      }else{
        return(xx$coefficients)
      }
    }
  }))
    
  df <- as.data.frame(df) %>%
    mutate(facility = facilities)
  
  facility_miss = as.data.frame(data) %>% 
    group_by(facility) %>%
    summarize(num_missing = sum(is.na(y)))
  
  df <- merge(facility_miss, df) %>%
    arrange(num_missing)
  
  return(df)
}
#
##### Imputation Functions #####

# OG imputation method
WF_imputation <- function(df, col, group = 'facility', family = 'NB', period = 12, R_PI = 500, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)){
  
  # check if this method has already been run
  if(any(grepl('y_pred_WF', colnames(df)))){
    print('previous WF predictions found. Removing them')
    df[,grep('y_pred_WF', colnames(df))] <- NULL
  }
  
  # prep the data with the harmonic functions
  df <- add_periodic_cov(df, period = period)
  
  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  

  if(family == 'quasipoisson'){
    warning('havent coded PIs in quasipoisson or NB yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=quasipoisson)
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    #df <- data.table::rbindlist(tmp)
    df <- do.call('rbind',tmp)
    
  }else if(family %in% c('NB','negative binomial','neg_bin')){
    warning('havent coded PIs in quasipoisson or NB yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      #mod_col <- MASS::glm.nb(formula_col, data = tt, control = glm.control(maxit=200,trace = 3))
      
      tt[,paste0(col, '_pred_WF')] <- tryCatch({
        mod_col <- MASS::glm.nb(formula_col, data = tt)
        predict(mod_col, tt, type = 'response')
      }, error = function(e){
        print(sprintf('glm.nb failed for %s', xx))
        rep(NA, nrow(tt))
      })
      #tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)
  }else if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=poisson)
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      
      if(R_PI > 0){
        # extract information from model fit
        beta_hat <- mod_col$coefficients
        beta_vcov <- vcov(mod_col)
        
        # bootstrap 
        sapply(1:R_PI, function(r){
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (tt %>% 
                          mutate(intercept=1) %>%
                          #dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          dplyr::select(intercept,year,cos1,sin1,cos2,sin2,cos3,sin3) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot) 
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
          x = rpois(n = nrow(tt), pred_boot_exp)
          
          return(x)
          
        }) -> sim.boot
        
        # get the quantiles and store them
        fitted_quants = t(apply(sim.boot, 1, function(xx){
          quantile(xx, probs = quant_probs)
        }))
        fitted_quants = as.data.frame(fitted_quants)
        colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
        
        tt = cbind(tt, fitted_quants)
        
      }
      tmp_lst = list(tt, sim.boot)
      return(tmp_lst)
    })
    
    # combine the individual facility results into a larger data frame
    #df <- data.table::rbindlist(lapply(tmp,'[[',1))
    df <- do.call('rbind',lapply(tmp,'[[',1))
    
    # take the sum of the bootstrap samples of each facility (this returns an n X R matrix itself)
    sim.full = Reduce('+', lapply(tmp, '[[', 2))
    
    # get the quantiles at the county level
    county_quants = t(apply(sim.full, 1, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    county_quants = as.data.frame(county_quants)
    colnames(county_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
    county_quants$date = tmp[[1]][[1]]$date

  }
  
  res_lst = list(df = df, county_fit = county_quants)
  return(res_lst)
}

# Run Weinberger-Fulcher (WF) model using all observed data with no missingness as a baseline comparison.
WF_baseline <- function(df, train_end_date = '2019-12-01', col = "y", family = 'poisson', R_PI = 100){
  # replace missing points with their true values
  tmp <- df
  tmp$y <- tmp$y_true
  tmp$y[tmp$date > train_end_date] <- NA
  
  res <- WF_imputation(tmp, col = col, family = family, R_PI = R_PI)
  
  # store the results and return the original y values
  tmp <- res$df
  tmp$y <- df$y
  
  # rename columns of the results
  colnames(tmp) <- gsub('y_pred_WF', 'y_pred_baseline_WF', colnames(tmp)) 
  
  return(tmp)
}

# Bayes imputation method
bayes_WF_imputation <- function(df, df_OG = NULL, col, group = 'facility', family = 'NB', period = 12, iterations = 500, harmonic_priors = F){
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
      pred_col_name = paste0(col, '_pred_bayes_WF_ideal')
    }else{
      # setting the prior to use the imputed matrix models
      df_prior = df
      pred_col_name = paste0(col, '_pred_bayes_WF_miss')
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
        tt[, pred_col_name] <- tryCatch({
          # fit the model
          mod_bayes <- stan_glm(formula_col, data = tt, family = neg_binomial_2, refresh = 0, iter = iterations, prior = beta_prior, prior_intercept = intercept_prior)
          
          # need to fill in NAs, even if they aren't used
          tmp2 = tt; tmp2[is.na(tmp2)] = 1e6
          
          # predict the median output from a posterior draw
          tt[, pred_col_name] = apply(posterior_predict(mod_bayes, tmp2, type = 'response'), 2, median)
        }, error = function(e){
          print(sprinft('error for facility %s. Returning NAs', xx))
          rep(NA, nrow(tt))
        })
        
        
        # if(xx == "Family Diagnostic Clinic"){
        #   browser()
        # }
        
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
        tt[, paste0(col, '_pred_bayes_WF')] = apply(posterior_predict(mod_bayes, tmp, type = 'response'), 2, median)
        
        return(tt)
      })
      
    }
    #browser()
    df <- data.table::rbindlist(tmp)
  }
  return(df)
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

CARBayes_imputation <- function(df, col, AR = 1, return_type = 'data.frame', model = c('fixed','facility_intercept','facility_fixed'), burnin = 20000, n.sample = 40000, prediction_sample = T){
  
  # check if this method has already been run
  if(any(grepl('y_CARBayes_ST', colnames(df)))){
    print('previous CAR Bayes predictions found. Removing them')
    df[,grep('y_CARBayes_ST', colnames(df))] <- NULL
  }
  
  tt = df %>% filter(date == unique(date)[1]) %>% pull(district) %>% table
  
  if(any(tt == 1)){
    warning('Randomly choosing a district for a facility without that district (this should be fixed later)')
    # replace the single districts with the biggest ones
    for(nn in names(which(tt == 1))){
      df$district = gsub(nn, names(which.max(tt)), df$district)
    }
  }
  
  districts = df %>% group_by(district) %>% summarize(n = length(unique(facility)))
  
  # create the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  W = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # model formula
  if(model == 'facility_intercept'){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + facility", col))
  }else if(model == 'facility_fixed'){
    formula_col = as.formula(sprintf("%s ~ facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3", col))
  }else{
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  }
  
  # run CAR Bayes
  chain1 <- ST.CARar(formula = formula_col, family = "poisson",
                     data = df, W = W, burnin = burnin, n.sample = n.sample,
                     thin = 10, AR = AR)
  
  df[,paste0(col, '_CARBayes_ST')] = chain1$fitted.values
  
  # Poisson sample the fitted values for the posterior predictive distribution
  if(prediction_sample){
    # pull the fitted values and randomly select prediction values based on the Poisson distribution
    tt = chain1$samples$fitted
    tt = apply(tt, 2, function(xx) rpois(n = length(xx), lambda = xx))
    
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(tt, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
    
  # ignore the sampling of the Poisson and take mean fitted value quantiles
  }else{
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(chain1$samples$fitted, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
  }
  
  df = as.data.frame(cbind(df, fitted_quants))

  # mapping matrix from each data point to a date
  dates = sort(unique(df$date))
  mat = as.matrix(sparseMatrix(i = 1:nrow(df), j = match(df$date, dates), x = 1))
  
  # condense the fitted value samples to each date
  date_fits = chain1$samples$fitted %*% mat
  
  # get the quantiles at the county level
  fitted_quants_C = t(apply(date_fits, 2, function(xx){
    quantile(xx, probs = quant_probs)
  }))
  fitted_quants_C = as.data.frame(fitted_quants_C)
  colnames(fitted_quants_C) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
  
  # aggregate the county results
  df_county = df %>% 
    dplyr::select(date, y_true, y_CARBayes_ST) %>% 
    group_by(date) %>%
    summarize(y_true = sum(y_true),
    y_CARBayes_ST = sum(y_CARBayes_ST)) %>% 
    arrange(date)
  
  if(!identical(df_county$date, dates)){
    browser()
  }

  df_county = cbind(df_county, fitted_quants_C)
  
  if(return_type == 'data.frame'){
    return(df)
  }else if(return_type == 'model'){
    return(chain1)
  }else if(return_type == 'all'){
    lst = list(facility_df = df, county_df = df_county, model_chain = chain1)
    return(lst)
  }
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
      WF_imputation(df_miss, col = col, family = harmonic_family, group = group)
    }, error = function(e){
      print(sprintf('error for periodic with p = %s', p))
      browser()
    })
    
    if(!is.na(bayes_harmonic_family)){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group)
      }, error = function(e){
        print(sprintf('error for bayes periodic with p = %s', p))
        browser()
      })
    }
    
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
      }, error = function(e){
        print(sprintf('error for bayes periodic beta prior with p = %s', p))
        browser()
      })
    }
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, df_OG = df, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
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

##### freqGLMepi Imputation Functions #####
model.mean <- function(D, params){
  mu = exp(params[1])*D$y.AR1 + # auto-regressive
    exp(params[2])*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

model.mean.exp <- function(D, params){
  mu = params[1]*D$y.AR1 + # auto-regressive
    params[2]*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

# likelihood function
ll.wrapper = function(params, D, target_col){
  mu = model.mean(D, params)
  logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
  return(-logL)
}

ll.wrapper.exp = function(params, D, target_col){
  mu = model.mean.exp(D, params)
  logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
  return(-logL)
}

freqGLMepi_variance <- function(param_vec, D){
  
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }
  
  # get the observed information matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  if(det(obs_I) < 0){
    return(list(variance = NULL, pos_def = F))
  }else{
    pos_def = T
  }
  
  # get the variance of the parameters
  V = tryCatch({solve(obs_I)}, error = function(e) {print(det(obs_I)); browser()})
  # V = solve(obs_I)
  
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
  #res_lst = list(variance = V, pos_def = pos_def)
  
  return(V)
}

freqGLMepi_variance2 <- function(param_vec, D){
  
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }
  
  ll.wrapper.exp = function(params, target_col = 'y_imp'){
    mu = model.mean.exp(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(-logL)
  }
  
  # get the observed information matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  # if(det(obs_I) < 0){
  #   return(list(variance = NULL, pos_def = F))
  # }else{
  #   pos_def = T
  # }
  
  # get the variance of the parameters
  #V = tryCatch({solve(obs_I)}, error = function(e) {print(det(obs_I)); browser()})
  V = solve(obs_I)
  
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
  #res_lst = list(variance = V, pos_def = pos_def)
  
  return(V)
}

freqGLMepi_variance_testing <- function(param_vec, D){
  
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }
  
  ll.wrapper.exp = function(params, target_col = 'y_imp'){
    mu = model.mean.exp(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(-logL)
  }
  
  
  # get the observed information matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  obs_2 = -pracma::hessian(ll.wrapper, param_vec)
  
  obs_grad = numDeriv::grad(ll.wrapper, param_vec)
  
  # for numerical stability of the variance components
  param_vec_exp = param_vec
  param_vec_exp[1:2] = exp(param_vec_exp[1:2])

  obs_test = -numDeriv::hessian(ll.wrapper.exp, param_vec_exp)
  obs_test_2 =  -pracma::hessian(ll.wrapper.exp, param_vec_exp)
  
  # need to do the delta method anyway
  
  ## Analytically, lezz do it.
  
  X = cbind(1, D[,c('year','cos1','sin1','cos2','sin2','cos3','sin3')])
 
  # ugly way - these's def some matrix algebra to speed this up but ignoring this for now
  logL_hess <- matrix(0, nrow = length(param_vec), ncol = length(param_vec))
  logL_grad <- matrix(0, nrow = length(param_vec), ncol = 1)
  mu = model.mean(D, param_vec)
  for(i in 2:nrow(D)){
    y_i = D$y_imp[i]
    x_i = as.matrix(t(X[i,]))
    
    mu_grad_i = as.matrix(unlist(c(exp(param_vec[1])*D$y.AR1[i],
                  exp(param_vec[2])*D$y.neighbors[i],
                  x_i*exp(t(x_i)%*%param_vec[3:10])[1,1])))
    
    mu_grad2_i = Matrix::bdiag(list(exp(param_vec[1])*D$y.AR1[i],
                           exp(param_vec[2])*D$y.neighbors[i],
                           x_i%*%t(x_i)*exp(t(x_i)%*%param_vec[3:10])[1,1]))

    logL_grad <- logL_grad + (y_i/mu[i] - 1)*mu_grad_i
    logL_hess <- logL_hess + (y_i/mu[i] - 1)*mu_grad2_i - y_i/(mu[i]^2)*mu_grad_i %*% t(mu_grad_i)
  }
  
  
  if(norm(obs_I + logL_hess)/norm(logL_hess) > .01){
    browser()
  }

  if(det(obs_I) < 0){
    if(det(logL_hess) < 0){
      browser()
    }else{
      obs_I = -logL_hess
    }
    
  }
  
  # get the variance of the parameters
  #V = tryCatch({solve(obs_I}, error = function(e) {browser()})})
  V = solve(obs_I)
  
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
  return(V)
}

### Fits the freqGLMepi model given the data
#   df: data
#   num_inits: number of different initial parameter tries
#   BFGS: Also try BFGS along with Nelder-Mead
#   verbose: Printing output or not
fit_freqGLMepi <- function(df, num_inits = 10, BFGS = T, verbose = T, target_col = 'y_imp', init = NULL){
  t0 = Sys.time()
  
  # set up initialization
  if(is.null(init)){
    init = rep(0,10)
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  # NM_control = list(maxit = 10000, reltol = 1e-12)
  NM_control = list(maxit = 10000)
  # BFGS_control = list(maxit = 10000, factr = 1e-11))
  BFGS_control = list(maxit = 10000)
  
  # fit using Nelder-Mead and L-BFGS-B and pick the better one
  tryCatch({
    params = optim(par = init, fn = ll.wrapper, D = df, target_col = target_col, control = NM_control)
  }, error = function(e){
    print(e)
    browser()
  })
  
  # L-BFGS
  if(BFGS){
    tryCatch({
      params2 = optim(init, ll.wrapper,  D = df, target_col = target_col, method = 'L-BFGS-B',  control = BFGS_control)
      if(params2$value < params$value & params2$convergence == 0){
        params = params2
      }
    }, error = function(e){
      #print(sprintf('skipping this round of L-BFGS-B because of error: %s', e))
    })
  }
  
  # try different initialization values to compare convergence
  if(num_inits > 1){
    for(i in 2:num_inits){
      # init = rnorm(10,0,10*i/num_inits)
      
      init = init_OG + rnorm(10,0,5*i/num_inits)
      
      # nelder-mead
      tryCatch({
        params2 = optim(init, ll.wrapper, D = df, target_col = target_col, control = NM_control)
      }, error = function(e){
        browser()
      })
      
      if(params2$value < params$value & params2$convergence == 0){
        if(verbose){print('using another initialization')}
        params = params2
      }
      
      # L-BFGS
      if(BFGS){
        tryCatch({
          params2 = optim(init, ll.wrapper,  D = df, target_col = target_col, method = 'L-BFGS-B',  control = BFGS_control)
          if(params2$value < params$value & params2$convergence == 0){
            if(verbose){print('using another initialization')}
            params = params2
          }
        }, error = function(e){
          if(verbose){
            print(sprintf('skipping this round of L-BFGS-B because of error: %s', e)) 
          }
        })
      }
    }
  }
  
  # for error checking
  if(params$convergence != 0){
    if(verbose){print('didnt converge for one iteration')}
    #browser()
  }
  
  if(verbose){
    print(sprintf('freqGLMepi fit in %s seconds', Sys.time() - t0))
  }
  
  names(params$par) = parnames
  return(params)
}

### the NNLS freqGLMepi model (so lambda and phi are constrained rather than exponentiated)
fit_freqGLMepi_nnls <- function(df, num_inits = 10, verbose = T, target_col = 'y_imp', init = NULL, BFGS = NULL){
  t0 = Sys.time()
  
  # set up initialization
  if(is.null(init)){
    init = c(0.1, 0.1, rep(0,8))
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  # params = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
  
  tryCatch({
    params = nlminb(start = init, objective = ll.wrapper.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
  }, error = function(e){
    browser()
  })
  
  # try different initialization values to compare convergence
  if(num_inits > 1){
    for(i in 2:num_inits){
      # init = rnorm(10,0,10*i/num_inits)
      
      init = init_OG + c(0.01*i, 0.01*i, rnorm(8,0,5*i/num_inits))
      
      # nelder-mead
      tryCatch({
        params2 = nlminb(start = init, objective = ll.wrapper.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
      }, error = function(e){
        browser()
      })
      
      if(params2$objective < params$objective & params2$convergence == 0){
        if(verbose){print('using another initialization')}
        params = params2
      }
    }
  }
  
  # for error checking
  if(params$convergence != 0){
    if(verbose){print('didnt converge for one iteration')}
    #browser()
  }
  
  if(verbose){
    print(sprintf('freqGLMepi fit in %s seconds', Sys.time() - t0))
  }
  
  names(params$par) = parnames
  return(params)
}

### Given data with missing values, fits the freqGLMepi model on all data points
#   df: data
#   max_iter: number of iterations of imputation
#   tol: tolerance of the converage of imputed values
#   individual_facility_models: fit each facility separately
#   prediction_intervals: whether to do prediction intervals and how to do them
#   R_PI: number of bootstrap iterations if doing so
#   quant_probs: the quantiles of the bootstrap to store in the data frame
#   verbose: printing updates
freqGLMepi_imputation = function(df, max_iter = 1, tol = 1e-4, individual_facility_models = T,  prediction_intervals= c('none','parametric_bootstrap','bootstrap','stationary_bootstrap'), R_PI = 100, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), refit_boot_outliers = F, verbose = F, optim_init = NULL, scale_by_num_neighbors = T, blocksize = 10, smart_boot_init = T, nnls = F){
  # check that we have the right columns
  if(!('y' %in% colnames(df) & 'y_true' %in% colnames(df))){
    stop('make sure the data has y (with NAs) and y_true')
  }
  
  # set fitting function to be regular (exponentiated spatial and temporal parameters) or non-negative least squares 
  if(nnls){
    fit_function <- fit_freqGLMepi_nnls
    model_function <- model.mean.exp
  }else{
    fit_function <- fit_freqGLMepi
    model_function <- model.mean
  }
  
  y_pred_list = list()
  
  ### Do initial filling of y
  # setting up the formula
  formula_col = as.formula(sprintf("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"))
  
  # the unique facility groups
  uni_group = unique(df$facility)
  
  #browser()
  # run the individual model for each group.
  tmp <- lapply(uni_group, function(xx) {
    tt <- df %>% filter(facility == xx)
    
    # run the model
    mod_col <- MASS::glm.nb(formula_col, data = tt)
    
    # update predictions
    tt$y_pred_freqGLMepi = predict(mod_col, tt, type = 'response')
    
    return(tt)
  })
  
  # combine into one data frame
  df = do.call('rbind',tmp)
  
  # filling in missing values by randomly sampling mean prediction from Poisson
  df$y_imp = df$y
  na.ind = which(is.na(df$y))
  df$y_imp[na.ind] <- rpois(n = length(na.ind), df$y_pred_freqGLMepi[na.ind])
  
  # add the neighbors and auto-regressive
  df = add_autoregressive(df, 'y_imp') %>%
    add_neighbors(., 'y_imp', scale_by_num_neighbors = scale_by_num_neighbors)
  
  ### Run freqGLM_epidemic model iteratively
  iter = 1
  y_pred_list[[1]] = df$y_pred_freqGLMepi
  prop_diffs = c(1)
  while(prop_diffs[length(prop_diffs)] > tol & iter <= max_iter){
    iter = iter + 1

    if(individual_facility_models){
      # run the individual model for each group.
      tmp <- lapply(uni_group, function(xx) {
        # subset data
        tt <- df %>% filter(facility == xx)
        
        # fit the model
        params = fit_function(tt, verbose = verbose, init = optim_init[[xx]])
        
        # update y_pred
        tt$y_pred_freqGLMepi = model_function(tt, params$par)
        
        return(list(df = tt, params = params))
      })
      
      # get matrix of the parameter estimates
      parmat = sapply(tmp, function(xx) xx[[2]]$par)
      
      # naming the parameter columns and rows
      rownames(parmat) = names(tmp[[1]]$params$par)
      colnames(parmat) = uni_group
      
      # get the convergence of each facility
      convergence = sapply(tmp, function(xx) (xx[[2]]$convergence == 0))
      
      # combine into one data frame
      df = do.call('rbind', lapply(tmp, '[[', 1)) %>% arrange(facility, date)
    }else{
      # fit the model
      parmat = fit_function(df)$par
      
      # update y_pred 
      df$y_pred_freqGLMepi = model_function(df, parmat)
      
    }
    
    # store the predictions for this iteration
    y_pred_list[[iter]] = df$y_pred_freqGLMepi
    
    if(length(na.ind) == 0){
      if(verbose){print('only running one iteration because there is no missingness')}
      break
    }
    
    # update y_imp
    na.ind.2 = intersect(na.ind, which(!is.na(df$y_pred_freqGLMepi)))
    df$y_imp[na.ind.2] <- rpois(n = length(na.ind.2), df$y_pred_freqGLMepi[na.ind.2])
    
    # compare y_imps
    prop_diffs = c(prop_diffs, mean(abs(y_pred_list[[iter]][na.ind] - y_pred_list[[iter-1]][na.ind])/y_pred_list[[iter-1]][na.ind], na.rm = T))
    
    # update the neighbors and auto-regressive
    df = add_autoregressive(df, 'y_imp') %>%
      add_neighbors(., 'y_imp')
    
    # update
    if(iter %% 10 == 0 & verbose){
      print(iter)
      if(individual_facility_models){
        print(parmat)
      }else{
        print(params$par)
      }
      
    }
  }
  
  if(prop_diffs[length(prop_diffs)] < tol){
    print(sprintf('convergence reached in %s iterations', iter))
  }
  
  num_errors = 0
  if(prediction_intervals == 'parametric_bootstrap'){
    warning('havent solved issue with non-positive-definite estimation of covariance')
    # estimate the variance of the MLEs
    V_list = list()
    
    # parametric bootstrap for prediction intervals
    # 1) For each facility,
    tmp <- lapply(1:length(uni_group), function(i) {
      # subset data
      tt <- df %>% filter(facility == uni_group[i])
      
      # get the coefficients and compute the variance
      coef_hat = parmat[,i]

      coef_vcov = tryCatch({freqGLMepi_variance(param_vec = parmat[,i], tt)
      }, error = function(e){
        return(NULL)
      })
      
      if(is.null(coef_vcov)){
        print('NULL!')
        return(NULL)
      }
      
      print('CALLING HAND-CALCULATED VARIANCE!')
      freqGLMepi_variance_testing(param_vec = parmat[,i], tt)
      #vcov_lst = freqGLMepi_variance(param_vec = parmat[,i], tt)
      # coef_vcov = vcov_lst$variance
      # 
      # # if not positive definite return an error
      # if(!vcov_lst$pos_def){
      #   return(NULL)
      # }
      
      # bootstrap 
      sim.boot <- tryCatch({
        sapply(1:R_PI, function(r){
        
        # sample coefficients and get the mean
        coef_boot <- MASS::mvrnorm(1,coef_hat,coef_vcov)
        pred_boot <- model_function(tt, coef_boot)
        
        pred_boot[which(is.na(pred_boot))] <- 1
        x = rpois(n = nrow(tt), pred_boot)
        
        return(x)})
        },
        error = function(e){
          return(-1)
        })

      if(length(sim.boot) == 1 | any(is.na(sim.boot))){
        return(NULL)
      }
      #browser()
      # get the quantiles and store them
      fitted_quants = t(apply(sim.boot, 1, function(xx){
        quantile(xx, probs = quant_probs)
      }))
      fitted_quants = as.data.frame(fitted_quants)
      colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)
      
      # combine the final results and return
      tt = cbind(tt, fitted_quants)
      tmp_lst = list(tt, sim.boot, coef_vcov)
      return(tmp_lst)
    })
    
    num_errors = sum(sapply(tmp, function(xx) is.null(xx[[1]])))
    vcov_list = lapply(tmp, function(xx) xx[[3]])
    df <- do.call('rbind',lapply(tmp,'[[',1))

  }else if(prediction_intervals == 'bootstrap'){
    warning('for bootstrap, only doing 1 initial value for optimization and only doing Nelder Mead')
    warning('for bootstrap, not testing for outliers')
    # 1) For each facility,
    tmp <- lapply(1:length(uni_group), function(i) {
      # subset data
      tt <- df %>% filter(facility == uni_group[i])
      
      # bootstrap 
      sapply(1:R_PI, function(r){
        
        # resample the data
        tt_boot = sample_n(tt, size = nrow(tt), replace = T)
        
        # fit the model
        params = fit_function(tt_boot, BFGS = F, num_inits = 1, verbose = verbose)
        
        # update y_pred sampled from the full data frame tt
        x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
        
        # test for outliers and then re-fit params (sometimes convergence can lead to crazy results)
        if(refit_boot_outliers){
          fit_iter = 1
          while(any(x > 10*median(x, na.rm = T), na.rm = T) & fit_iter <= 3){
            
            if(verbose){
              print(sort(x))
              print('rerunning fitting')
            }
            params = fit_function(tt_boot, BFGS = F, num_inits = 10^(fit_iter), verbose = verbose)
            x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
            
            fit_iter = fit_iter + 1
          }
          
          # check if we still have outliers
          if(any(x > 10*median(x, na.rm = T), na.rm = T)){
            print('outlier found for one bootstrap iteration')
          }
        }
        
        return(x)
        
      }) -> sim.boot
      
      # get the quantiles and store them
      fitted_quants = t(apply(sim.boot, 1, function(xx){
        quantile(xx, probs = quant_probs, na.rm = T)
      }))
      fitted_quants = as.data.frame(fitted_quants)
      colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)
      
      # combine the final results and return
      tt = cbind(tt, fitted_quants)
      tmp_lst = list(tt, sim.boot)
      return(tmp_lst)
    })
    
    df <- do.call('rbind',lapply(tmp,'[[',1))
  }else if(prediction_intervals == 'stationary_bootstrap'){
    
    tmp <- lapply(1:length(uni_group), function(i) {
      # subset data
      tt <- df %>% filter(facility == uni_group[i])
      
      ## function wrapper to predict values from the model
      if(smart_boot_init){
        # start the initialization at the values from the main model
        predict_function <- function(xx){
          # fit model on bootstrapped data set
          params = fit_function(xx, BFGS = F, num_inits = 1, verbose = verbose, init = parmat[,1])
          
          # predict model on original data set
          x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
          return(x)
        }
      }else{
        predict_function <- function(xx){
          # fit model on bootstrapped data set
          params = fit_function(xx, BFGS = F, num_inits = 1, verbose = verbose)
          
          # predict model on original data set
          x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
          return(x)
        }
      }

      
      # run the stationary bootstrap
      sim.boot <- boot::tsboot(tt, statistic = predict_function, R = R_PI, sim = 'geom', l = blocksize)$t
      
      # get the quantiles and store them
      fitted_quants = t(apply(sim.boot, 2, function(xx){
        quantile(xx, probs = quant_probs, na.rm = T)
      }))
      fitted_quants = as.data.frame(fitted_quants)
      colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)
      
      # combine the final results and return
      tt = cbind(tt, fitted_quants)
      tmp_lst = list(tt, sim.boot)
      return(tmp_lst)
    })
    
    # combining the final data frame
    df <- do.call('rbind',lapply(tmp,'[[',1))
    
    # for returning
    vcov_list = NULL
  }
  
  # prep data to return
  return_lst = list(df = df, params = parmat, convergence = convergence, y_pred_list = y_pred_list, prop_diffs = prop_diffs, num_errors = num_errors, vcov_list = vcov_list)
  return(return_lst)
}

#
##### Plotting + Metric Functions #####

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

plot_county_imputations <- function(df, imp_vec, color_vec, title = 'Aggregated County Imputations', labels = NULL){
  # set the labels for the plot
  if(is.null(labels)){
    labels = imp_vec
  }
  
  # rename columns of df to make it easier
  for(col in imp_vec){
    ind = grep(col, colnames(df))
    if(length(ind) != 1){browser()}
    colnames(df)[ind] = col
  }
  
  # initialize the data frame to store final results
  df_f = NULL
  
  for(j in 1:length(imp_vec)){
    col = imp_vec[j]
    
    # create the imputed outcome
    df$imputed_outcome = df$indicator_count_ari_total
    df$imputed_outcome[is.na(df$imputed_outcome)] = df[is.na(df$imputed_outcome), col]
    
    # aggregate by date  
    tmp = df %>% 
      group_by(date) %>%
      summarize(y = sum(imputed_outcome)) %>% mutate(method = col)
    
    # store results for this method
    df_f = rbind(df_f, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_f$method = factor(df_f$method, levels = imp_vec)
  
  # plot them all!
  p1 <- ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
    geom_line() +
    ggtitle(title) + 
    scale_color_manual(values = color_vec, labels = labels) + theme(legend.position = 'bottom', legend.text=element_text(size=20))
  
  return(p1)
}

plot_county_fits_from_fac <- function(df, imp_vec, color_vec, title = 'Aggregated County Fits', labels = NULL){
  
  # set the labels for the plot
  if(is.null(labels)){
    labels = imp_vec
  }
  
  # # rename columns of df to make it easier
  # for(col in imp_vec){
  #   ind = grep(col, colnames(df))
  #   if(length(ind) != 1){browser()}
  #   colnames(df)[ind] = col
  # }
  
  # initialize the data frame to store final results
  df_f = NULL

  for(j in 1:length(imp_vec)){
    col = imp_vec[j]

    # aggregate by date
    tmp = df %>%
      group_by(date) %>%
      summarize(y = sum(get(col))) %>% mutate(method = col)

    # store results for this method
    df_f = rbind(df_f, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_f$method = factor(df_f$method, levels = imp_vec)
  
  # plot them all!
  p1 <- ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
    geom_line() +
    ggtitle(title) + 
    scale_color_manual(values = color_vec, labels = labels) + theme(legend.position = 'bottom', legend.text=element_text(size=20))
  
  return(p1)
}

plot_imputations <- function(df, imp_vec, color_vec, fac_list = NULL){
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # rename columns of df to make it easier
  for(col in imp_vec){
    ind = grep(col, colnames(df))
    if(length(ind) != 1){browser()}
    colnames(df)[ind] = col
  }
  
  # initialize plotting
  plot_list = list()
  iter = 0
  
  # go through each facility
  for(f in fac_list){
    iter = iter + 1
    tmp = df %>% filter(facility == f)
    
    # store the true outcome for plotting
    df_f = tmp %>% dplyr::select(date, y = indicator_count_ari_total) %>% mutate(method = 'ari_true')
    
    #  # baseline plot for this
    # # p1 <- suppressWarnings(ggplot(tmp, aes(x = date, y = indicator_count_ari_total)) + 
    #                           geom_line() +
    #                           ggtitle(sprintf('%s (%0.1f %% M)', f, 100*mean(is.na(tmp$indicator_count_ari_total)))))
    for(j in 1:length(imp_vec)){
      col = imp_vec[j]
      tmp[!is.na(tmp$indicator_count_ari_total),col] = NA
      
      # filling in surrounding points to imputations to connect lines
      # I have to do multiple if statements because of R's handling of calling a value not in an array
      ind_list = c()
      for(i in 1:nrow(tmp)){
        if(i > 1){
          if(is.na(tmp[i,col]) & !is.na(tmp[i-1,col])){
            ind_list = c(ind_list, i)
          }
        }
        if(i < nrow(tmp)){
          if(is.na(tmp[i,col]) & !is.na(tmp[i+1,col])){
            ind_list = c(ind_list, i)
          }
        }
      }
      # replace surrounding values for plotting
      tmp[ind_list, col] = tmp$indicator_count_ari_total[ind_list]
      
      tmp2 = tmp %>% dplyr::select(date, y = imp_vec[j]) %>% mutate(method = col)
      df_f = rbind(df_f, tmp2)
      
      #scale_colour_manual(name="Error Bars",values=cols)
      #p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = imp_vec[j]), color = color_vec[j]))
      # p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = imp_vec[j], color = color_vec[j])))
    }
    
    #browser()
    
    p1 <- suppressWarnings(ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
                             geom_line() +
                             ggtitle(sprintf('%s (%0.1f%% M)', f, 100*mean(is.na(tmp$indicator_count_ari_total)))) + 
                             scale_color_manual(values = c('black', color_vec)))
    
    # store the legend for later
    legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
    
    # remove the legend position on this plot
    p1 = p1 + theme(legend.position = 'none') 
    
    # store the plot for this facility in the list
    plot_list[[iter]] = p1
  }
  
  #browser()
  plot_grid(plot_grid(plotlist = plot_list),  legend, ncol = 1, rel_heights = c(10,1))
  
}

# HERE! FOR COMMENTING PURPOSES. HAVENT COMMENTED THE FUNCTIONS ABOVE

plot_county_fits <- function(df, imp_vec, color_vec, imp_names = NULL, PIs = T, title = 'County-Level Predictions'){
  df = as.data.frame(df)
  
  # if no imputation names given, use the ones in the imputation vector
  if(is.null(imp_names)){
    imp_names = imp_vec
  }
  
  # initialize data frame for this county
  df_c = NULL
  
  # pull estimates for each method
  for(j in 1:length(imp_vec)){
    col = imp_vec[j]
    
    # get the lower and upper bounds
    tmp = df[,c('date',paste0(col,'_0.5'),paste0(col, '_0.025'),paste0(col, '_0.975'))] 
    colnames(tmp) = c('date', 'y', 'y_lower', 'y_upper')
    tmp$method = imp_names[j]
    
    df_c = rbind(df_c, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_c$method = factor(df_c$method, levels = imp_names)
  
  # make the plot!
  p1 <- ggplot() +
    geom_line(data = df, aes(x = date, y = y_true), size = 1) +
    geom_line(data = df_c, aes(x = date, y = y, group = method, color = method)) +
    geom_ribbon(data = df_c, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.3) +
    scale_color_manual(values = c(color_vec)) + 
    scale_fill_manual(values = c(color_vec)) + 
    ggtitle(title) + 
    ylab('y') +
    theme_bw() + 
    theme(text = element_text(size = 15))
  
  # store the legend for later
  legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
  
  # remove the legend position on this plot
  #p1 = p1 + theme(legend.position = 'none') 
  return(p1)
}

# plot the fits of the models for a group of facilities. Also plots the prediction intervals.
plot_facility_fits <- function(df, imp_vec, imp_names = NULL, color_vec, PIs = T, fac_list = NULL, plot_missing_points = T){
  df = as.data.frame(df)
  
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # if no imputation names given, use the ones in the imputation vector
  if(is.null(imp_names)){
    imp_names = imp_vec
  }
  
  # initialize plotting
  plot_list = list()
  iter = 0
  
  # go through each facility
  for(f in fac_list){
    iter = iter + 1
    tmp = df %>% filter(facility == f)
    
    if(!is.null(imp_vec)){
      # initialize data frame for this facility
      df_f = NULL
      
      for(j in 1:length(imp_vec)){
        col = imp_vec[j]
        
        # get the lower and upper bounds
        tmp2 = tmp[,c('date',col,paste0(col, '_0.025'),paste0(col, '_0.975'))] 
        colnames(tmp2) = c('date', 'y', 'y_lower', 'y_upper')
        tmp2$method = imp_names[j]
        
        df_f = rbind(df_f, tmp2)
      }
      
      # ordering the method to be consistent and for the labeling
      df_f$method = factor(df_f$method, levels = imp_names)
      
      
      # make the plot!
      p1 <- ggplot() +
        geom_line(data = tmp, aes(x = date, y = y_true), size = 1) +
        geom_line(data = df_f, aes(x = date, y = y, group = method, color = method)) +
        geom_ribbon(data = df_f, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.1) +
        scale_color_manual(values = c(color_vec)) + 
        scale_fill_manual(values = c(color_vec)) + 
        ylim(c(0,1.5*max(tmp$y_true))) + 
        ggtitle(sprintf('facility %s', f)) + 
        ylab('y') +
        theme_bw() +
        theme(text = element_text(size = 10))
      
      if(plot_missing_points){
        tmp2 <- tmp %>%
          filter(is.na(y)) %>%
          select(date)
        
        p1 <- p1 + 
          geom_point(data = tmp2, aes(x = date, y = 0),  color = 'red', size = 3)
      }
      
      # store the legend for later
      legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
      
      # remove the legend position on this plot
      p1 = p1 + theme(legend.position = 'none') 
      
      # store the plot for this facility in the list
      plot_list[[iter]] = p1
    }else{
      print('plotting baseline counts because no imputation methods are provided')
      
      p1 <- ggplot() +
        geom_line(data = tmp, aes(x = date, y = y), size = 1) + 
        ylim(c(0,1.5*max(tmp$y))) + 
        ggtitle(sprintf('facility %s', f)) + 
        ylab('y') +
        theme_bw() +
        theme(text = element_text(size = 10))
        
      plot_list[[iter]] = p1
    }
    
  }
  
  plot_grid(plot_grid(plotlist = plot_list),legend, ncol = 1, rel_heights = c(10,1))
}

# calculate the desired metrics for a set of imputation methods
calculate_metrics <- function(df, imp_vec, imputed_only = T, median_estimate = F){
  
  warning('this function incorrectly aggregates by simulation, not by data point, so dont trust these results')
  
  df = as.data.frame(df)
  
  # if imputed only, compute metrics only on the missing values
  if(imputed_only){
    df = df %>% filter(is.na(y))
  }
  
  
  # if using median as the point estimate for each method
  if(median_estimate){
    # for each method, compute metrics
    tmp_lst <- lapply(imp_vec, function(xx){
      tmp = data.frame(metric = c('coverage95','coverage50','bias','absolute_bias','RMSE', 'MAPE'),
                       value = c(mean(df$y_true >= df[,paste0(xx, '_0.025')] & df$y_true <= df[,paste0(xx, '_0.975')]),
                                 mean(df$y_true >= df[,paste0(xx, '_0.25')] & df$y_true <= df[,paste0(xx, '_0.75')]),
                                 mean(df[,paste0(xx, '_0.5')] - df$y_true),
                                 mean(abs(df[,paste0(xx, '_0.5')] - df$y_true)),
                                 sqrt(mean((df[,paste0(xx, '_0.5')] - df$y_true)^2)),
                                 mean(abs(df[,paste0(xx, '_0.5')] - df$y_true)/df$y_true, na.rm = T)))
      colnames(tmp)[2] = xx
      return(tmp)
    })
      
  }else{
    # for each method, compute metrics
    tmp_lst <- lapply(imp_vec, function(xx){
      tmp = data.frame(metric = c('coverage95','coverage50','bias','absolute_bias','RMSE', 'MAPE'),
                       value = c(mean(df$y_true >= df[,paste0(xx, '_0.025')] & df$y_true <= df[,paste0(xx, '_0.975')]),
                                 mean(df$y_true >= df[,paste0(xx, '_0.25')] & df$y_true <= df[,paste0(xx, '_0.75')]),
                                 mean(df[,xx] - df$y_true),
                                 mean(abs(df[,xx] - df$y_true)),
                                 sqrt(mean((df[,xx] - df$y_true)^2)),
                                 mean(abs(df[,paste0(xx, '_0.5')] - df$y_true)/df$y_true, na.rm = T)))
      colnames(tmp)[2] = xx
      
      return(tmp)
    })
  }
  
  # combine them all
  res = do.call(merge, tmp_lst)
  
  return(res)
}

# plot a set of metrics for a list of simulation runs
plot_metrics_bysim <- function(imputed_list, imp_vec, rename_vec = NULL, color_vec = c('red','blue'), imputed_only = T){
  
  warning('this function incorrectly aggregates by simulation, not by data point, so dont trust these results')
  
  # get the metrics from each simulation run
  res_full = NULL
  for(r in 1:length(imputed_list)){
    tmp <- calculate_metrics(imputed_list[[r]], imp_vec, imputed_only = imputed_only)
    tmp$r = r
    res_full <- rbind(res_full, tmp)
  }
  
  # convert to long form
  long = tidyr::gather(res_full, method, value, y_pred_harmonic:y_CARBayes_ST)
  
  # replace the names
  if(!is.null(rename_vec)){
    for(i in 1:length(imp_vec)){
      long$method = gsub(imp_vec[i], rename_vec[i],long$method)
    }
  }
  
  plot_list = list()
  # plot each metric
  for(i in 1:length(unique(long$metric))){
    
    m = unique(long$metric)[i]
    
    tmp2 = long %>% filter(metric == m)
    
    p1 <- ggplot(tmp2, aes(x = method, y = value, fill = method)) +  #, color = method)
      geom_violin(position="dodge", alpha=0.5) + 
      geom_jitter(position = position_jitter(0.2)) + 
      scale_color_manual(values = color_vec) +
      scale_fill_manual(values = color_vec) + 
      theme_bw() + 
      theme(legend.position = 'none') +
      ggtitle(m)
    
    plot_list[[i]] = p1
  }
  
  # make the final plot
  final_plot <- plot_grid(plotlist = plot_list)
  
  return(final_plot)
}

# calculate the metrics for individual data points across simulated imputations
calculate_metrics_by_point <- function(imputed_list, imp_vec = c("y_pred_WF", "y_CARBayes_ST"), imputed_only = T, rm_ARna = F, use_point_est = F){
  
  #y_true = imputed_list[[1]]$y_true
  
  # removing the starting points with NA AR1 values, since these
  if(rm_ARna){
    print('removing the starting points because of NA autoregressive term')
    imputed_list = lapply(imputed_list, function(xx) xx[!is.na(xx$y.AR1),])
  }
  
  # make sure the dates and facilities match
  date_facility_check(imputed_list)
  
  # get the true values everywhere and at the deleted time points
  y_true = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y_true']))
  y_missing = do.call('cbind', lapply(imputed_list, function(xx) {
    y_true = xx[,'y_true'];
    y_true[!is.na(xx[,'y'])] = NA
    y_true
  }))
  # ^above, an NA means the value was not missing, and a number means the value was missing
  
  # numeric missing matrix
  missing_mat <- apply(y_missing, 2, function(xx) 1 - as.numeric(is.na(xx)))
  missing_mat_NA <- missing_mat; missing_mat_NA[missing_mat_NA == 0] <- NA
  
  # get the number of times each data point was missing across simulations
  num_missing = apply(y_missing, 1, function(xx) sum(!is.na(xx)))
  
  if(!imputed_only){
    df = NULL
    # imputed only here
    for(method in imp_vec){
      
      tmp = imputed_list[[1]][,c('date','facility','district')]
      tmp$method = method
      tmp$num_missing = num_missing
      
      point_est = do.call('cbind', lapply(imputed_list, function(xx) xx[,method]))
      lower_025 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.025')]))
      upper_975 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.975')]))
      lower_25 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.25')]))
      upper_75 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.75')]))
      median = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.5')]))
      
      if(use_point_est){
        outcome = point_est
      }else{
        outcome = median
      }
      
      # full dataset
      tmp$y_missing <- sapply(1:nrow(y_missing), function(ii) mean(y_missing[ii,], na.rm = T))
      tmp$y_true <- sapply(1:nrow(y_true), function(ii) mean(y_true[ii,], na.rm = T))
      tmp$median <- sapply(1:nrow(median), function(ii) mean(median[ii,]))
      tmp$point_est <- sapply(1:nrow(point_est), function(ii) mean(point_est[ii,]))
      tmp$bias = rowMeans(sapply(1:ncol(outcome), function(ii) {outcome[,ii] - y_true[,ii]}))
      tmp$absolute_bias = rowMeans(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])}))
      tmp$MAPE = rowMeans(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])/y_true[,ii]}))
      tmp$RMSE = sqrt(rowMeans(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_true[,ii])^2})))
      
      tmp$coverage50 = rowMeans(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] > lower_25[,ii] & y_true[,ii] < upper_75[,ii])))
      tmp$coverage95 = rowMeans(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] > lower_025[,ii] & y_true[,ii] < upper_975[,ii])))
      
      tmp$lower_025 <- sapply(1:nrow(lower_025), function(ii) median(lower_025[ii,], na.rm = T))
      tmp$upper_975 <- sapply(1:nrow(upper_975), function(ii) median(upper_975[ii,], na.rm = T))
      
      # measure of how wide the 95% prediction intervals are
      tmp$interval_width = rowMeans(upper_975 - lower_025) 
      tmp$prop_interval_width = rowMeans((upper_975 - lower_025)/y_true)
      tmp$point_sd = apply(outcome, 1, sd)
      
      # update the results
      df = rbind(df, tmp)
    }
    
  # imputed only here
  }else{
    df = NULL
    
    for(method in imp_vec){
      
      tmp = imputed_list[[1]][,c('date','facility','district')]
      tmp$method = method
      tmp$num_missing = num_missing
      
      point_est = do.call('cbind', lapply(imputed_list, function(xx) xx[,method]))
      lower_025 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.025')]))
      upper_975 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.975')]))
      lower_25 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.25')]))
      upper_75 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.75')]))
      median = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.5')]))
      
      # make points only for missing data set
      point_est <- point_est*missing_mat_NA
      median <- median*missing_mat_NA
      lower_025 <- lower_025*missing_mat_NA
      lower_25 <- lower_25*missing_mat_NA
      upper_75 <- upper_75*missing_mat_NA
      upper_975 <- upper_975*missing_mat_NA
      
      if(use_point_est){
        outcome = point_est
      }else{
        outcome = median
      }
      
      # missing data set
      tmp$y_missing <- sapply(1:nrow(y_missing), function(ii) mean(y_missing[ii,], na.rm = T))
      tmp$y_true <- sapply(1:nrow(y_true), function(ii) mean(y_true[ii,], na.rm = T))
      tmp$median <- sapply(1:nrow(median), function(ii) mean(median[ii,], na.rm = T))
      tmp$point_est <- sapply(1:nrow(point_est), function(ii) mean(point_est[ii,], na.rm = T))
      tmp$bias = rowMeans(sapply(1:ncol(outcome), function(ii) {outcome[,ii] - y_missing[,ii]}), na.rm = T)
      tmp$absolute_bias = rowMeans(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_missing[,ii])}), na.rm = T)
      tmp$MAPE = rowMeans(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_missing[,ii])/y_missing[,ii]}), na.rm = T)
      tmp$RMSE = sqrt(rowMeans(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_missing[,ii])^2}), na.rm = T))
      
      tmp$coverage50 = rowMeans(sapply(1:ncol(lower_25), function(ii) (y_missing[,ii] > lower_25[,ii] & y_missing[,ii] < upper_75[,ii])), na.rm = T)
      tmp$coverage95 = rowMeans(sapply(1:ncol(lower_25), function(ii) (y_missing[,ii] > lower_025[,ii] & y_missing[,ii] < upper_975[,ii])), na.rm = T)
      
      tmp$lower_025 <- sapply(1:nrow(lower_025), function(ii) median(lower_025[ii,], na.rm = T))
      tmp$upper_975 <- sapply(1:nrow(upper_975), function(ii) median(upper_975[ii,], na.rm = T))
      
      # get the interval width and standard deviation only at points that are missing
      tmp$interval_width = rowMeans(do.call('cbind', lapply(imputed_list, function(xx){
        interval = xx[,paste0(method, '_0.975')] - xx[,paste0(method, '_0.025')];
        interval[!is.na(xx[,'y'])] = NA
        interval
      })), na.rm = T)
      
      tmp$prop_interval_width = rowMeans(do.call('cbind', lapply(imputed_list, function(xx){
        interval = xx[,paste0(method, '_0.975')] - xx[,paste0(method, '_0.025')];
        interval[!is.na(xx[,'y'])] = NA
        interval/xx$y_true
      })), na.rm = T)
      
      # if(method == "y_CB_facility"){
      #   browser()
      # }
 
      tmp$point_sd = apply(do.call('cbind', lapply(imputed_list, function(xx){
        median = xx[,paste0(method, '_0.5')];
        median[!is.na(xx[,'y'])] = NA
        median
      })), 1, function(xx) sd(xx, na.rm = T))
      
      # update the results
      df = rbind(df, tmp)
    }
  }
  
  #res_lst = list(df = df, num_missing = num_missing)
  return(df)
}

# plot the metrics by individual data points across simulated imputations
plot_metrics_by_point <- function(imputed_list, imp_vec = c('y_pred_WF', 'y_CARBayes_ST'), rename_vec = NULL, color_vec = c('red','blue'), imputed_only = T, outcomes = NULL, min_missing = 5, rm_ARna = F, use_point_est = F, violin_points = F, max_intW_lim = NULL, metric_list = c('bias','absolute_bias','MAPE','RMSE','coverage50','coverage95','interval_width','prop_interval_width')){
  
  if(is.null(outcomes)){
    outcomes = c("bias", "absolute_bias", "MAPE", "RMSE", "coverage50", "coverage95", "interval_width", "point_sd")
  }
  
  # get the metrics from each simulation run
  res = calculate_metrics_by_point(imputed_list, imp_vec = imp_vec, imputed_only = imputed_only, rm_ARna = rm_ARna, use_point_est = use_point_est)

  if(all(res$num_missing < min_missing) & imputed_only){
    stop('all points are not missing enough (check min_missing)')
  }else if(any(res$num_missing < min_missing) & imputed_only){
    # remove the data points that werent missing enough
    ind = which(res$num_missing < min_missing)
    res = res[-ind,]
    
    print(sprintf('removing %s data points because they werent missing %s times', length(ind), min_missing))
  }
  
  # convert to long form
  long = tidyr::gather(res, metric, value, bias:point_sd)
  
  # replace the names
  if(!is.null(rename_vec)){
    for(i in 1:length(imp_vec)){
      long$method = gsub(imp_vec[i], rename_vec[i],long$method)
    }
  }
  
  # make a factor so the ordering in the plots stays
  long$method=  factor(long$method, levels = rename_vec)
  
  long = long %>% 
    filter(metric %in% metric_list)
  
  plot_list = list()
  # plot each metric
  for(i in 1:length(unique(long$metric))){
    
    m = unique(long$metric)[i]
    
    tmp2 = long %>% filter(metric == m)
    
    p1 <- ggplot(tmp2, aes(x = method, y = value, fill = method)) +  #, color = method)
      geom_violin(position="dodge", alpha=0.5) + 
      scale_color_manual(values = color_vec) +
      scale_fill_manual(values = color_vec) + 
      theme_bw() + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank()) +
      ggtitle(m)
    
    # stash the legend
    if(i == 1){
      legend = get_legend(p1 + theme(legend.position = 'bottom'))
    }
    
    # remove the legend for this plot
    p1 <- p1 + theme(legend.position = 'none')
    
    if(violin_points){
      p1 <- p1 + geom_jitter(position = position_jitter(0.1))
    }
    
    if(m == 'bias'){
      p1 <- p1 + geom_hline(yintercept = 0)
    }else if(m == 'coverage50'){
      p1 <- p1 + geom_hline(yintercept = 0.5) + ylim(0,1)
    }else if(m == 'coverage95'){
      p1 <- p1 + geom_hline(yintercept = 0.95) + ylim(0,1)
    }else if(m == 'interval_width'){
      tt = long %>% filter(method != 'glmFreq_epi', metric == 'interval_width')
      if(is.null(max_intW_lim)){
        max_intW_lim = round(1.3*max(tt$value, na.rm = T))
      }
      p1 = p1 + ylim(c(0, max_intW_lim))
    }
    
    plot_list[[i]] = p1
  }
  
  # make the final plot
  final_plot <- plot_grid(plot_grid(plotlist = plot_list, nrow = 2), legend, ncol = 1, rel_heights = c(10,1))
  
  return(final_plot)
}

#
##### Simulation functions #####
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

initialize_df <- function(district_sizes, start_date = '2016-01-01', end_date = '2019-12-01'){
  facilities = unlist(lapply(1:length(district_sizes), function(xx) {paste0(toupper(letters[xx]), 1:district_sizes[xx])}))
  
  dates = seq(as.Date(start_date), as.Date(end_date), by = 'month')
  
  df = expand.grid(facilities, dates, stringsAsFactors = F)
  colnames(df) = c('facility','date')
  
  df$district = substr(df$facility, 1, 1) 
  
  df = df %>%
    dplyr::select(date, facility, district) %>%
    arrange(facility, date)
  return(df)
}

sample_real_betas <- function(facilities, file = 'results/all_facility_betas_filtered_09052022.csv'){
  tmp <- read.csv(file)
  
  ind <- sample(nrow(tmp), length(facilities))
  betas <- tmp[ind,3:10]
  rownames(betas) = facilities
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
  return(betas)
}

sample_betas = function(facilities, b0_mean = 4.3, b1_mean = -0.25, b1_sd = 0.26){
  betas = matrix(0, nrow = length(facilities), ncol = 8)
  
  betas[,1] = rnorm(b0_mean, 1, n = nrow(betas))
  betas[,2] = rnorm(b1_mean, b1_sd, n = nrow(betas))
  
  for(j in 3:8){
    betas[,j] = rnorm(0, 0.15, n = nrow(betas))
  }
  
  rownames(betas) = facilities
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
  #paste0('B',0:7)
  #betas = cbind(facilities, as.data.frame(betas))
  return(betas)
}

simulate_data <- function(district_sizes, R = 1, empirical_betas = F, ...){
  # set up data frame
  df = initialize_df(district_sizes, ...)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility names
  facilities = unique(df$facility)
  
  # set random seed and sample betas
  set.seed(10)
  if(empirical_betas){
    betas = sample_real_betas(facilities)
  }else{
    betas = sample_betas(facilities)  
  }
  
  # initialize list of data frames
  df_lst = list()
  
  # make n sampled sets of data
  for(i in 1:R){
    
    # simulate values given the betas
    tmp_lst = lapply(facilities, function(xx){
      tmp = df %>% filter(facility == xx)
      
      # keep the 1 for intercepts
      X = tmp %>% 
        mutate(intercept = 1) %>%
        dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
      
      # error checking
      if(!identical(colnames(betas), colnames(X))){
        browser()
      }
      
      # make the 8x1 beta vector for this facility
      beta_f = t(betas[xx,,drop = F])
      
      # get mean prediction from linear model
      mu = as.matrix(X)%*%beta_f
      
      # simluate random values
      tmp$y = rpois(length(mu), exp(mu))
      
      return(tmp)
    })
    
    # combine values into one data frame
    df_lst[[i]] = do.call('rbind', tmp_lst)
    
  }
  
  # make list of values to return
  res_lst = list(df_list = df_lst, betas = betas)
  return(res_lst)
}

simulate_data_spatiotemporal <- function(district_sizes, R = 1, rho = 0.3, alpha = 0.3, tau = 1, scale_by_num_neighbors = T, b0_mean = 7, b1_mean = -0.2){
  # set up data frame
  df = initialize_df(district_sizes)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility names
  facilities = unique(df$facility) %>% sort()
  
  # get all dates
  dates = unique(df$date) %>% sort()
  #n = length(dates)
  
  # set random seed and sample betas
  set.seed(10)
  betas = sample_betas(facilities, b0_mean = b0_mean, b1_mean = b1_mean)
  
  ### get the mean effects (since these don't change across the simulation samples)
  df = do.call('rbind', lapply(facilities, function(xx){
    tmp = df %>% filter(facility == xx)
    
    # keep the 1 for intercepts
    X = tmp %>% 
      mutate(intercept = 1) %>%
      dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
    
    # error checking
    if(!identical(colnames(betas), colnames(X))){
      browser()
    }
    
    # make the 8x1 beta vector for this facility
    beta_f = t(betas[xx,,drop = F])
    
    # get mean prediction from linear model
    tmp$mu = (as.matrix(X)%*%beta_f)[,1]
    
    return(tmp)
  }))
  
  # make the spatio-temporal precision matrix
  Q = make_precision_mat(df, rho = rho)
  
  tryCatch({
    V = tau^2*solve(Q)
  }, error = function(e){
    print(e)
    print('havent dealt with non-invertible precision matrices yet')
    browser()
  })
  
  # checking ordering of facilities matches
  if(!identical(colnames(V), facilities)){
    stop('names of covariances and facilities dont match')
  }
  
  # make R sampled sets of data
  df_lst = lapply(1:R, function(i){
    ### get the spatio-temporal random effects
    # initialize phi
    phi = matrix(0, nrow = length(dates), ncol = length(facilities))
    colnames(phi) = facilities
    
    # first time step
    phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
    
    # cycle through other time steps, using auto-correlated priors
    for(i in 2:length(dates)){
      phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
    }
    
    # convert to matching format
    phi = as.data.frame(phi)
    phi$date = dates
    phi_df = tidyr::gather(phi, facility, phi, setdiff(colnames(phi), c('facility','date')))
    
    # merge the phi values into the original data frame
    df = merge(df, phi_df, by = c('date','facility'))
    
    # simulate the observed values
    df$y = rpois(nrow(df), exp(df$mu + df$phi))
    
    return(df)
  })
  
  # make list of values to return
  res_lst = list(df_list = df_lst, betas = betas, V = V, rho = rho, alpha = alpha, tau = tau, b0_mean = b0_mean, b1_mean = b1_mean)
  return(res_lst)
}


#### Function to simulate data under the freqGLM_epi framework
##
## district_sizes = number of facilities in the districtss
## R = number of simulated data sets
## lambda = autoregressive term
## phi = neighbor term
simulate_data_freqGLM_epi <- function(district_sizes, R = 1, lambda = -2, phi = -2, num_iters = 10, scale_by_num_neighbors = F, seed = 10, start_date = '2016-01-01', end_date = '2019-12-01', b0_mean = 7, b1_mean = -0.2, b1_sd = 0.2){
  
  warning('counting all districts as neighbors')
  print(sprintf('lambda: exp(%0.2f) = %0.2f; phi: exp(%0.2f) = %0.2f', lambda, exp(lambda), phi, exp(phi)))
  
  # set up data frame
  df = initialize_df(district_sizes, start_date = start_date, end_date = end_date)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility names
  facilities = unique(df$facility) %>% sort()
  
  # get all dates
  dates = unique(df$date) %>% sort()
  #n = length(dates)
  
  # set random seed and sample betas
  set.seed(seed)
  betas = sample_betas(facilities, b0_mean = b0_mean, b1_mean = b1_mean, b1_sd = b1_sd)
  
  ### get the seasonal effects (since these don't change across the simulation samples)
  df = do.call('rbind', lapply(facilities, function(xx){
    tmp = df %>% filter(facility == xx)
    
    # keep the 1 for intercepts
    X = tmp %>% 
      mutate(intercept = 1) %>%
      dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
    
    # error checking
    if(!identical(colnames(betas), colnames(X))){
      browser()
    }
    
    # make the 8x1 beta vector for this facility
    beta_f = t(betas[xx,,drop = F])
    
    # get mean prediction from linear model
    tmp$mu_seasonal = (as.matrix(X)%*%beta_f)[,1]
    
    return(tmp)
  }))
  
  # initial sampling
  df$y_seasonal = rpois(n = nrow(df), lambda = exp(df$mu_seasonal))
  
  # initialize list of final data frame
  df_lst = list()
  
  for(r in 1:R){
    
    # to store the predicted values at successive iterations
    y_pred_list = list()
    
    tmp = df %>%
      mutate(y = y_seasonal)
    
    for(i in 1:num_iters){
      # add the neighbors and auto-regressive
      tmp = add_autoregressive(tmp, 'y') %>%
        add_neighbors(., 'y', scale_by_num_neighbors = scale_by_num_neighbors)
      
      # only resampling after the first time point
      ind = which(!is.na(tmp$y.AR1))
      
      # getting the mean according to the model and resampling
      mean_y = exp(lambda)*tmp$y.AR1 + exp(phi)*tmp$y.neighbors + exp(tmp$mu_seasonal)
      tmp$y[ind] = rpois(n = length(ind), lambda = mean_y[ind])
      
      # storing the results
      y_pred_list[[i]] = tmp$y
      
    }
    
    if(mean(y_pred_list[[num_iters]]) > 5*mean(y_pred_list[[1]])){
      print('there might be some divergence of the estimated values')
      browser()
    }
    
    df_lst[[r]] = tmp
  }

  params_true = as.data.frame(t(betas))
  rownames(params_true) = paste0('B', rownames(params_true))
  params_true = rbind(t(data.frame(By.AR1 = rep(lambda, ncol(params_true)), By.neighbors = rep(phi, ncol(params_true)), row.names = colnames(params_true))), params_true)
  
  
  # return it!
  res_lst = list(df_list = df_lst, betas = betas, lambda = lambda, phi = phi, params = params_true)
  return(res_lst)
}

make_precision_mat <- function(df, rho){
  # check the rho value
  if(rho < 0 | rho >= 1){
    stop('please input a rho in [0,1)')
  }
  
  # create the list of matching facilities
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  pairs = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district)
  
  # get unique facilities
  facilities = unique(df$facility) %>% sort()
  
  # initialize the W2 matrix (note that what I am calling W here is what the original paper calls diag(W1) - W)
  W2 = matrix(0, nrow = length(facilities), ncol = length(facilities))
  colnames(W2) = facilities
  rownames(W2) = facilities
  
  for(i in 1:length(facilities)){
    f1 = facilities[i]
    
    # get the pairs with this facility
    tmp = pairs %>% filter(facility.x == f1)
    
    if(nrow(tmp) == 0){
      print('havent done single-district facilities yet')
      browser()
    }else{
      # put -1 where there are pairs of facilities
      matid = match(tmp$facility.y, facilities)
      W2[i,matid] = -1
    }
    
    # get the number of neighbors this facility has
    W2[i,i] = pairs %>% filter(facility.x == f1) %>% nrow()
    # print(sprintf('%s: NUM neighbors = %s', f1, W2[i,i]))
  }
    
  # make the final Q matrix
  Q = rho*W2 + (1-rho)*diag(rep(1, length(facilities)))
  
  return(Q)
}

MCAR_sim <- function(df, p, by_facility = F, max_missing_date = '2019-12-01'){
  # save the true y value
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)

  # delete y values to add in missingness
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  # combine the hold-out/non-missing data with the training data with missingness
  df <- rbind(df, df_test)
  
  return(df)
}

MNAR_sim <- function(df, p, direction = NULL, gamma = 1.5, by_facility = T){
  df$y_true = df$y
  if(by_facility){
    # get the number of points to impute
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    
    # if removing low and high points
    if(is.null(direction)){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = median(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    # if removing high points
    }else if(direction == 'upper'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = min(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    # if removing low points
    }else if(direction == 'lower'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = max(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    }
    
  }else{
    print('havent coded this part - is it necessary?')
    browser()
    #num_impute = round(p*nrow(df))
    #df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  return(df)
}

MAR_spatiotemporal_sim <- function(df, p, rho = 0.3, alpha = 0.3, tau = 1, by_facility = T){
  # make phi for df
  # for all, or for each facility,
  # sample according to expit(phi)
  df$y_true = df$y
  
  # get all facility names
  facilities = unique(df$facility) %>% sort()
  
  # get all dates
  dates = unique(df$date) %>% sort()
  
  # make the spatio-temporal precision matrix
  Q = make_precision_mat(df, rho = rho)
  
  tryCatch({
    V = tau^2*solve(Q)
  }, error = function(e){
    print(e)
    print('havent dealt with non-invertible precision matrices yet')
    browser()
  })
  
  # checking ordering of facilities matches
  if(!identical(colnames(V), facilities)){
    stop('names of covariances and facilities dont match')
  }
  
  phi = matrix(0, nrow = length(dates), ncol = length(facilities))
  colnames(phi) = facilities
  
  # first time step
  phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
  
  # cycle through other time steps, using auto-correlated priors
  for(i in 2:length(dates)){
    phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
  }
  
  # convert to matching format
  phi = as.data.frame(phi)
  phi$date = dates
  phi_df = tidyr::gather(phi, facility, phiM, setdiff(colnames(phi), c('facility','date')))
  
  # convert to expit for probability
  phi_df$prob.sample = exp(phi_df$phiM)/(1 + exp(phi_df$phiM))
  
  # merge the phi values into the original data frame
  df = merge(df, phi_df, by = c('date','facility'))
  
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute, prob = tmp$prob.sample)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute, prob = df$prob.sample)] <- NA
  }
  
  return(df)
}
