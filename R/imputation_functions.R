### The functions used for imputation in Liberia

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
periodic_imputation <- function(df, col, group = 'facility', family = 'NB', period = 12, R_PI = 500){
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
      tt[,paste0(col, '_pred_harmonic')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)
    
  }else if(family %in% c('NB','negative binomial','neg_bin')){
    warning('havent coded PIs in quasipoisson or NB yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      #mod_col <- MASS::glm.nb(formula_col, data = tt, control = glm.control(maxit=200,trace = 3))
      
      tt[,paste0(col, '_pred_harmonic')] <- tryCatch({
        mod_col <- MASS::glm.nb(formula_col, data = tt)
        predict(mod_col, tt, type = 'response')
      }, error = function(e){
        print(sprintf('glm.nb failed for %s', xx))
        rep(NA, nrow(tt))
      })
      #tt[,paste0(col, '_pred_harmonic')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)
  }else if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=poisson)
      tt[,paste0(col, '_pred_harmonic')] = predict(mod_col, tt, type = 'response')
      
      if(R_PI > 0){
        # extract information from model fit
        beta_hat <- mod_col$coefficients
        beta_vcov <- vcov(mod_col)
        
        # bootstrap 
        sapply(1:R_PI, function(r){
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (tt %>% 
                          mutate(intercept=1) %>%
                          dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot) 
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
          x = rpois(n = nrow(tt), pred_boot_exp)
          
          return(x)
          
        }) -> sim.boot
        
        # get the quantiles and store them
        quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
        fitted_quants = t(apply(sim.boot, 1, function(xx){
          quantile(xx, probs = quant_probs)
        }))
        fitted_quants = as.data.frame(fitted_quants)
        colnames(fitted_quants) = paste0(paste0(col, '_pred_harmonic_'), quant_probs)
        
        tt = cbind(tt, fitted_quants)
        
      }
      
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
        tt[, paste0(col, '_pred_bayes_harmonic')] = apply(posterior_predict(mod_bayes, tmp, type = 'response'), 2, median)
        
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

CARBayes_imputation <- function(df, col, AR = 1, return_type = 'data.frame', facility_intercept = T){
  df = add_periodic_cov(df)
  
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
  if(facility_intercept){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + facility", col))
  }else{
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  }
  
  
  # run CAR Bayes
  chain1 <- ST.CARar(formula = formula_col, family = "poisson",
                     data = df, W = W, burnin = 20000, n.sample = 40000,
                     thin = 10, AR = AR)
  
  df[,paste0(col, '_CARBayes_ST')] = chain1$fitted.values
  
  # get the quantiles of fitted values 
  quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  fitted_quants = t(apply(chain1$samples$fitted, 2, function(xx){
    quantile(xx, probs = quant_probs)
  }))
  fitted_quants = as.data.frame(fitted_quants)
  colnames(fitted_quants) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
  
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

plot_county_fits <- function(df, imp_vec, color_vec, title = 'Aggregated County Fits', labels = NULL){
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

# plot the fits of the models for a group of facilities. Also plots the prediction intervals.
plot_facility_fits <- function(df, imp_vec, color_vec, PIs = T, fac_list = NULL){
  df = as.data.frame(df)
  
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # initialize plotting
  plot_list = list()
  iter = 0
  
  # go through each facility
  for(f in fac_list){
    iter = iter + 1
    tmp = df %>% filter(facility == f)
    
    # initialize data frame for this facility
    df_f = NULL
    
    for(j in 1:length(imp_vec)){
      col = imp_vec[j]
      
      # get the lower and upper bounds
      tmp2 = tmp[,c('date',col,paste0(col, '_0.025'),paste0(col, '_0.975'))] 
      colnames(tmp2) = c('date', 'y', 'y_lower', 'y_upper')
      tmp2$method = col

      df_f = rbind(df_f, tmp2)
    }
    
    # ordering the method to be consistent and for the labeling
    df_f$method = factor(df_f$method, levels = imp_vec)
    
    # make the plot!
    p1 <- ggplot() +
      geom_line(data = tmp, aes(x = date, y = y_true), size = 1) +
      geom_line(data = df_f, aes(x = date, y = y, group = method, color = method)) +
      geom_ribbon(data = df_f, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.3) +
      scale_color_manual(values = c(color_vec)) + 
      scale_fill_manual(values = c(color_vec)) + 
      ggtitle(sprintf('%s', f))
    
    # store the legend for later
    legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
    
    # remove the legend position on this plot
    p1 = p1 + theme(legend.position = 'none') 
    
    # store the plot for this facility in the list
    plot_list[[iter]] = p1
  }
  
  #browser()
  plot_grid(plot_grid(plotlist = plot_list),legend, ncol = 1, rel_heights = c(10,1))
}

# calculate the desired metrics for a set of imputation methods
calculate_metrics <- function(df, imp_vec, imputed_only = T){
  df = as.data.frame(df)
  
  # if imputed only, compute metrics only on the missing values
  if(imputed_only){
    df = df %>% filter(is.na(y))
  }
  
  # for each method, compute metrics
  tmp_lst <- lapply(imp_vec, function(xx){
    tmp = data.frame(metric = c('coverage95','coverage50','bias','absolute_bias','RMSE'),
                     value = c(mean(df$y_true >= df[,paste0(xx, '_0.025')] & df$y_true <= df[,paste0(xx, '_0.975')]),
                               mean(df$y_true >= df[,paste0(xx, '_0.25')] & df$y_true <= df[,paste0(xx, '_0.75')]),
                               mean(df[,xx] - df$y_true),
                               mean(abs(df[,xx] - df$y_true)),
                               sqrt(mean((df[,xx] - df$y_true)^2))))
    colnames(tmp)[2] = xx
    
    return(tmp)
  })
  
  # combine them all
  res = do.call(merge, tmp_lst)
  
  return(res)
}
