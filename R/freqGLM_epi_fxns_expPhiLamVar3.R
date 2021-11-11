## Now finding the variance of the beta estimates all together. Is this legal? I'm not sure but I'm going to try it


# Nvm this is dumb. This can't work because of the nonlinearity of the seasonal coefficients/because there are the X values in the exponent

model.mean <- function(D, params){
  mu = params[1]*D$y.AR1 + # auto-regressive
    params[2]*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

freqGLMepi_variance <- function(param_vec, D){
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }

 # param_vec[1:2] = c(exp(param_vec[1:2]), exp(param_vec[3] + param_vec[4]*D$year + param_vec[5]*D$cos1 + param_vec[6]*D$sin1 + param_vec[7]*D$cos2 + param_vec[8]*D$sin2 + param_vec[9]*D$cos3 + param_vec[10]*D$sin3))
  
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
  #browser()
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
 # res_lst = list(variance = V, pos_def = pos_def)
  
  return(V)
}

fit_freqGLMepi <- function(df, num_inits = 10, verbose = T, target_col = 'y_imp', init = NULL){
  t0 = Sys.time()
  
  # set up initialization
  if(is.null(init)){
    init = c(0.1, 0.1, rep(0,8))
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  params = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))

  # fit using Nelder-Mead and L-BFGS-B and pick the better one
  tryCatch({
    params = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
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
        params2 = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
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


