

# returns the predicted mean at each data point
# params is a vector with the first value being alpha for the autoregressive term, the second value being beta for the spatial term, and the last 8 values being the gamma values for the seasonal terms
model.mean <- function(D, params){
  mu = exp(params[1])*D$y.AR1 + # auto-regressive
    exp(params[2])*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

# given the predicted means and y's, returns the likelihood (excluding the [log y!] portion)
ll = function(mu, y){
  ll = sum(-mu + y*log(mu), na.rm = T)
  return(ll)
}


ll.wrapper = function(params, target_col = 'y_imp'){
  mu = model.mean(D, params)
  logL = ll(mu, D[,target_col])
  return(-logL)
}

params

test = model.mean(D, params)

params = rep(0,10)

D = na.omit(D)


test = optim(params, ll.wrapper)


### This function imputes missing values based on the surveillance method.

# to deal with missing values this fills in missing values iteratively: starting with the freqGLM model, fill in the missing values, run this model, then fill in again, then run, etc..., until the missing values are close enough.

lst <- simulate_data_spatiotemporal(district_sizes = c(4), R = 1, rho = 0.3, alpha = 0.5, tau = 0.5)

df = lst$df_list[[1]]
df = MCAR_sim(df, p = 0.2, by_facility = T)

surveillance_imputation = function(df, max_iter = 50, tol = 1e-4, verbose = T){
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
    
    # fill in the missing values from the data frame
    tt$y_pred = predict(mod_col, tt, type = 'response')
    
    return(tt)
  })
  
  # combine into one data frame
  df = do.call('rbind',tmp)
  
  # filling in missing values
  df$y_imp = df$y
  na.ind = which(is.na(df$y))
  df$y_imp[na.ind] <- df$y_pred[na.ind] 
  warning('currently only using mean-level imputations, not random predictive values')
  
  # add the neighbors and auto-regressive
  df = add_autoregressive(df, 'y_imp') %>%
    add_neighbors(., 'y_imp')
  
  ### Run freqGLM_epidemic model iteratively
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(df, params)
    logL = ll(mu, df[,target_col])
    return(-logL)
  }
  
  iter = 1
  y_pred_list[[1]] = df$y_pred
  prop_diffs = c(1)
  while(prop_diffs[length(prop_diffs)] > tol & iter < max_iter){
    iter = iter + 1
    
    # set up initialization
    init = rep(0,10)
    names(init) = c('AR.1', 'NE.1', 'B0', 'B1', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
    
    # fit using Nelder-Mead and L-BFGS-B and pick the better one
    params = optim(init, ll.wrapper, control = list(maxit = 5000))
    
    # L-BFGS
    tryCatch({
      params2 = optim(init, ll.wrapper, method = 'L-BFGS-B', control = list(maxit = 5000))
      if(params2$value < params$value & params2$convergence == 0){
        params = params2
      }
    }, error = function(e){
      print(sprintf('skipping this round of L-BFGS-B because of error: %s', e))
    })
    
    # try different initialization values to compare convergence
    for(i in 1:10){
      init = rnorm(10,0,i/3)
      # nelder-mead
      params2 = optim(init, ll.wrapper, control = list(maxit = 5000))
      if(params2$value < params$value & params2$convergence == 0){
        params = params2
      }
      
      # L-BFGS
      tryCatch({
        params2 = optim(init, ll.wrapper, method = 'L-BFGS-B', control = list(maxit = 5000))
        if(params2$value < params$value & params2$convergence == 0){
          params = params2
        }
      }, error = function(e){
        print(sprintf('skipping this round of L-BFGS-B because of error: %s', e))
      })
    }
    
    # THIS IS FOR ALL FACILITIES AT ONCE DINGUS.
    # The question is, do we specify the model to run for all facilities? Or do we run it for each facility individually and then just use the primary imputations for the other facilities

    # for error checking
    if(params$convergence != 0){
      browser()
    }
    
    # update y_pred 
    y_pred_list[[iter]] = model.mean(df, params$par)
    df$y_pred = y_pred_list[[iter]]
    
    # update y_imp
    na.ind.2 = intersect(na.ind, which(!is.na(df$y_pred)))
    df$y_imp[na.ind.2] <- df$y_pred[na.ind.2] 

    # compare y_imps
    prop_diffs = c(prop_diffs, mean(abs(y_pred_list[[iter]][na.ind] - y_pred_list[[iter-1]][na.ind])/y_pred_list[[iter-1]][na.ind], na.rm = T))
    
    # update the neighbors and auto-regressive
    df = add_autoregressive(df, 'y_imp') %>%
      add_neighbors(., 'y_imp')
    
    # update
    if(iter %% 10 == 0 & verbose){
      print(iter)
      print(params$par)
    }
  }
  
  if(prop_diffs[length(prop_diffs)] < tol){
    print(sprintf('convergence reached in %s iterations', iter))
  }
  print(prop_diffs)
  
  # prep data to return
  return_lst = list(df = df, params = params, y_pred_list = y_pred_list, prop_diffs = prop_diffs)
  return(return_lst)
  
}

surveillance_imputation(df, max_iter = 100) -> test



test = test$df

par(mfrow = c(2,2))
for(fac in unique(test$facility)){
  tmp = test %>% filter(facility == fac) %>% arrange(date)
  plot(tmp$date, tmp$y, type = 'l', main = fac)
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

# (2) Try different inits

# (3) Figure out how to get parameter standard error estimates





