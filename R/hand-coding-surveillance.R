
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')

library(dplyr)
library(lubridate)

# returns the predicted mean at each data point
# params is a vector with the first value being alpha for the autoregressive term, the second value being beta for the spatial term, and the last 8 values being the gamma values for the seasonal terms
model.mean <- function(D, params){
  mu = exp(params[1])*D$y.AR1 + # auto-regressive
    exp(params[2])*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component

  return(mu)
}



# given the predicted means and y's, returns the likelihood (excluding the [log y!] portion)
# ll = function(mu, y){
#   ll = sum(-mu + y*log(mu), na.rm = T)
#   return(ll)
# }

ll.wrapper = function(params, target_col = 'y_imp'){
  mu = model.mean(D, params)
  #logL = ll(mu, D[,target_col])
  logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
  return(-logL)
}

freqGLMepi_variance <- function(param_vec, D){
  
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

### This function imputes missing values based on the surveillance method.

# to deal with missing values this fills in missing values iteratively: starting with the freqGLM model, fill in the missing values, run this model, then fill in again, then run, etc..., until the missing values are close enough.

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 1, rho = 0.3, alpha = 0.5, tau = 0.5)

df = lst$df_list[[1]]
df = MCAR_sim(df, p = 0.2, by_facility = T)

### Fits the freqGLMepi model given the data
#   df: data
#   num_inits: number of different initial parameter tries
#   verbose: Printing output or not
fit_freqGLMepi <- function(df, num_inits = 10, BFGS = T, verbose = T){
  t0 = Sys.time()
  
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(df, params)
    logL = sum(-mu + df[,target_col]*log(mu), na.rm = T)
    return(-logL)
  }
  
  # set up initialization
  init = rep(0,10)
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  
  # NM_control = list(maxit = 10000, reltol = 1e-12)
  NM_control = list(maxit = 5000)
  # BFGS_control = list(maxit = 10000, factr = 1e-11))
  BFGS_control = list(maxit = 5000)
  # fit using Nelder-Mead and L-BFGS-B and pick the better one
  params = optim(init, ll.wrapper, control = NM_control)
  
  # L-BFGS
  if(BFGS){
    tryCatch({
      params2 = optim(init, ll.wrapper, method = 'L-BFGS-B',  control = BFGS_control)
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
      init = rnorm(10,0,i/3)
      # nelder-mead
      params2 = optim(init, ll.wrapper,  control = NM_control)
      if(params2$value < params$value & params2$convergence == 0){
        params = params2
      }
      
      # L-BFGS
      if(BFGS){
        tryCatch({
          params2 = optim(init, ll.wrapper, method = 'L-BFGS-B',  control = BFGS_control)
          if(params2$value < params$value & params2$convergence == 0){
            params = params2
          }
        }, error = function(e){
          print(sprintf('skipping this round of L-BFGS-B because of error: %s', e))
        })
      }
    }
  }
  
  # for error checking
  if(params$convergence != 0){
    browser()
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

freqGLMepi_imputation = function(df, max_iter = 1, tol = 1e-4, individual_facility_models = T,  prediction_intervals= c('none','parametric_bootstrap','bootstrap'), R_PI = 100, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),  verbose = T){
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
    tt$y_pred = predict(mod_col, tt, type = 'response')
    
    return(tt)
  })
  
  # combine into one data frame
  df = do.call('rbind',tmp)
  
  # filling in missing values by randomly sampling mean prediction from Poisson
  df$y_imp = df$y
  na.ind = which(is.na(df$y))
  df$y_imp[na.ind] <- rpois(n = length(na.ind), df$y_pred[na.ind])

  # add the neighbors and auto-regressive
  df = add_autoregressive(df, 'y_imp') %>%
    add_neighbors(., 'y_imp')
  
  ### Run freqGLM_epidemic model iteratively
  iter = 1
  y_pred_list[[1]] = df$y_pred
  prop_diffs = c(1)
  while(prop_diffs[length(prop_diffs)] > tol & iter <= max_iter){
    iter = iter + 1

    if(individual_facility_models){
      # run the individual model for each group.
      tmp <- lapply(uni_group, function(xx) {
        # subset data
        tt <- df %>% filter(facility == xx)
        
        # fit the model
        params = fit_freqGLMepi(tt, verbose = F)
        
        # update y_pred
        tt$y_pred = model.mean(tt, params$par)
        
        return(list(df = tt, params = params))
      })
      
      # get matrix of the parameter estimates
      parmat = sapply(tmp, function(xx) xx[[2]]$par)
      
      # naming the parameter columns and rows
      rownames(parmat) = names(tmp[[1]]$params$par)
      colnames(parmat) = uni_group
      
      # combine into one data frame
      df = do.call('rbind', lapply(tmp, '[[', 1)) %>% arrange(facility, date)
    }else{
      # fit the model
      parmat = fit_freqGLMepi(df)$par
      
      # update y_pred 
      df$y_pred =model.mean(df, parmat)
      
    }

    # store the predictions for this iteration
    y_pred_list[[iter]] = df$y_pred
    
    # update y_imp
    na.ind.2 = intersect(na.ind, which(!is.na(df$y_pred)))
    df$y_imp[na.ind.2] <- rpois(n = length(na.ind.2), df$y_pred[na.ind.2])

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
    browser()
  }else if(prediction_intervals == 'bootstrap'){
    warning('for bootstrap, only doing 1 initial value for optimization and only doing Nelder Mead')
    # 1) For each facility,
    tmp <- lapply(1:length(uni_group), function(i) {
      # subset data
      tt <- df %>% filter(facility == uni_group[i])

      # bootstrap 
      sapply(1:R_PI, function(r){
        
        # resample the data
        tt_boot = sample_n(tt, size = nrow(tt), replace = T)
        
        # fit the model
        params = fit_freqGLMepi(tt_boot, BFGS = T, num_inits = 1, verbose = F)
        
        # update y_pred sampled from the full data frame tt
        x = suppressWarnings(rpois(n = nrow(tt), model.mean(tt, params$par)))
        
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
  return_lst = list(df = df, params = parmat, y_pred_list = y_pred_list, prop_diffs = prop_diffs)
  return(return_lst)
}

test <- freqGLMepi_imputation(df, prediction_intervals = 'bootstrap', R_PI = 100) 


####################

test = test$df

par(mfrow = c(2,2))
for(fac in unique(test$facility)){
  tmp = test %>% filter(facility == fac) %>% arrange(date)
  plot(tmp$date, tmp$y, type = 'l', main = fac, col = 'black')
  #lines(tmp$date, tmp$y_imp, type = 'l', col = 'blue')
  lines(tmp$date, tmp$y_pred, col = 'red')
}

# so it's reasonable but it seems to be over-predicting the auto-correlation


### TO DO

# (1) Test that these results match the surveillance package, roughly speaking
head(df)
# so I'll just use y_true for y
df$y = df$y_true
# my_results <- surveillance_imputation(df, max_iter = 1000) <- doesnt work without any missingness


D2 = df %>% dplyr::select(district, facility) %>% distinct()
W = full_join(D2, D2, by = 'district') %>%
  filter(facility.x != facility.y) %>%
  dplyr::select(-district) %>%
  igraph::graph_from_data_frame() %>%
  igraph::as_adjacency_matrix() %>%
  as.matrix()

ARI.counts <- df %>% 
  dplyr::select(date, facility, y_true) %>%
  tidyr::spread(facility,y_true) %>% 
  arrange(date) %>%
  dplyr::select(-date) %>%
  as.matrix()

ARI <- sts(ARI.counts, start = c(2016, 1), frequency = 12, neighbourhood = W)

f.end <- addSeason2formula(f = ~ 1 + I(t/12), S = 3, period = 12)

model.1 <- list(ar = list(f = ~ 1),
                ne = list(f = ~ 1,
                          weights = neighbourhood(ARI),
                          normalize = T),
                end = list(f = f.end),
                family = 'Poisson', verbose = T,
                optimizer = list(variance = list(method = 'Nelder-Mead')))

result.1 <- hhh4(ARI, model.1)

# (3) Figure out how to get parameter standard error estimates





