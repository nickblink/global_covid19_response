
model.mean <- function(D, params){
  mu = exp(params[1] + params[2]*D$year + params[3]*D$cos1 + params[4]*D$sin1 + params[5]*D$cos2 + params[6]*D$sin2 + params[7]*D$cos3 + params[8]*D$sin3) # yearly + seasonal component
  
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
    init = rep(0,8)
  }
  
  parnames = c('Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
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
      
      init = init_OG + rnorm(8,0,10*i/num_inits)
      
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
