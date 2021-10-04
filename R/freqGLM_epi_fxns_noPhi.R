
model.mean <- function(D, params){
  mu = exp(params[1])*D$y.AR1 + # auto-regressive
    exp(params[2] + params[3]*D$year + params[4]*D$cos1 + params[5]*D$sin1 + params[6]*D$cos2 + params[7]*D$sin2 + params[8]*D$cos3 + params[9]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

# likelihood function
ll.wrapper = function(params, D, target_col){
  mu = model.mean(D, params)
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
  
  # get the observed infromation matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  obs_2 = -pracma::hessian(ll.wrapper, param_vec)
  
  #browser()
  
  if(det(obs_I) < 0){
    browser()
  }
  
  # get the variance of the parameters
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
    init = rep(0,9)
  }
  
  parnames = c('By.AR1', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  # init = c(-3.98630374, -1.72730745, 7.04425288, -0.13527676, -0.34012130, -0.07970224, -0.17471972, -0.10980406, -0.22505447, -0.04618626)
  
  # init = c(-2, -2, log(mean(df$y_imp)), rep(0,7))
  
  # NM_control = list(maxit = 10000, reltol = 1e-12)
  NM_control = list(maxit = 10000)
  # BFGS_control = list(maxit = 10000, factr = 1e-11))
  BFGS_control = list(maxit = 10000)
  
  # fit using Nelder-Mead and L-BFGS-B and pick the better one
  tryCatch({
    params = optim(par = init, fn = ll.wrapper, D = df, target_col = target_col, control = NM_control)
  }, error = function(e){
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
      
      init = init_OG + rnorm(9,0,10*i/num_inits)
      
      # nelder-mead
      tryCatch({
        params2 = optim(init, ll.wrapper, D = df, target_col = target_col, control = NM_control)
      }, error = function(e){
        browser()
      })
      
      if(params2$value < params$value & params2$convergence == 0){
        print('using another initialization')
        params = params2
      }
      
      # L-BFGS
      if(BFGS){
        tryCatch({
          params2 = optim(init, ll.wrapper,  D = df, target_col = target_col, method = 'L-BFGS-B',  control = BFGS_control)
          if(params2$value < params$value & params2$convergence == 0){
            print('using another initialization')
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
    print('didnt converge for one iteration')
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
freqGLMepi_imputation = function(df, max_iter = 1, tol = 1e-4, individual_facility_models = T,  prediction_intervals= c('none','parametric_bootstrap','bootstrap'), R_PI = 100, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), refit_boot_outliers = F, verbose = T, optim_init = NULL){
  # check that we have the right columns
  if(!('y' %in% colnames(df) & 'y_true' %in% colnames(df))){
    stop('make sure the data has y (with NAs) and y_true')
  }
  
  y_pred_list = list()
  
  ### Do initial filling of y
  # setting up the formula
  formula_col = as.formula(sprintf("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"))
  
  # the unique facility groups
  uni_group = unique(df$facility)
  
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
    add_neighbors(., 'y_imp')
  
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
        params = fit_freqGLMepi(tt, verbose = F, init = optim_init[[xx]])
        
        # update y_pred
        tt$y_pred_freqGLMepi = model.mean(tt, params$par)
        
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
      parmat = fit_freqGLMepi(df)$par
      
      # update y_pred 
      df$y_pred_freqGLMepi = model.mean(df, parmat)
      
    }
    
    # store the predictions for this iteration
    y_pred_list[[iter]] = df$y_pred_freqGLMepi
    
    if(length(na.ind) == 0){
      print('only running one iteration because there is no missingness')
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
  
  if(prediction_intervals == 'parametric_bootstrap'){
    warning('havent solved issue with non-positive-definite estimation of covariance')
    # estimate the variance of the MLEs
    V_list = list()
    
    # parametric bootstrap for prediction intervals
    # 1) For each facility,
    tmp <- lapply(1:length(uni_group), function(i) {
      print(i)
      # subset data
      tt <- df %>% filter(facility == uni_group[i])
      
      # get the coefficients and compute the variance
      coef_hat = parmat[,i]
      coef_vcov = freqGLMepi_variance(param_vec = parmat[,i], tt)
      
      print(det(coef_vcov))
      
      # if(i == 3){
      #   browser()
      # }
      
      # bootstrap 
      sapply(1:R_PI, function(r){
        
        # sample coefficients and get the mean
        coef_boot <- MASS::mvrnorm(1,coef_hat,coef_vcov)
        pred_boot <- model.mean(tt, coef_boot)
        
        pred_boot[which(is.na(pred_boot))] <- 1
        x = rpois(n = nrow(tt), pred_boot)
        
        return(x)
        
      }) -> sim.boot
      
      # get the quantiles and store them
      fitted_quants = t(apply(sim.boot, 1, function(xx){
        quantile(xx, probs = quant_probs)
      }))
      fitted_quants = as.data.frame(fitted_quants)
      colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)
      
      # combine the final results and return
      tt = cbind(tt, fitted_quants)
      tmp_lst = list(tt, sim.boot)
      return(tmp_lst)
    })
    
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
        params = fit_freqGLMepi(tt_boot, BFGS = F, num_inits = 1, verbose = verbose)
        
        # update y_pred sampled from the full data frame tt
        x = suppressWarnings(rpois(n = nrow(tt), model.mean(tt, params$par)))
        
        # test for outliers and then re-fit params (sometimes convergence can lead to crazy results)
        if(refit_boot_outliers){
          fit_iter = 1
          while(any(x > 10*median(x, na.rm = T), na.rm = T) & fit_iter <= 3){
            
            if(verbose){
              print(sort(x))
              print('rerunning fitting')
            }
            params = fit_freqGLMepi(tt_boot, BFGS = F, num_inits = 10^(fit_iter), verbose = verbose)
            x = suppressWarnings(rpois(n = nrow(tt), model.mean(tt, params$par)))
            
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
  }
  
  # prep data to return
  return_lst = list(df = df, params = parmat, convergence = convergence, y_pred_list = y_pred_list, prop_diffs = prop_diffs)
  return(return_lst)
}

#

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
      mean_y = exp(lambda)*tmp$y.AR1 + exp(tmp$mu_seasonal)
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
  params_true = rbind(t(data.frame(By.AR1 = rep(lambda, ncol(params_true)), row.names = colnames(params_true))), params_true)
  
  
  # return it!
  res_lst = list(df_list = df_lst, betas = betas, lambda = lambda, phi = phi, params = params_true)
  return(res_lst)
}