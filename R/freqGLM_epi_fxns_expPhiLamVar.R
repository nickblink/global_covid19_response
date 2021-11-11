model.mean2 <- function(D, params){
  mu = params[1]*D$y.AR1 + # auto-regressive
    params[2]*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}


freqGLMepi_variance <- function(param_vec, D){
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean2(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }

  param_vec[1:2] = exp(param_vec[1:2])
  
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
