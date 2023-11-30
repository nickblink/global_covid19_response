## This script makes bash commands for given simulations

bash_command <- function(p = NULL, b0_mean = 6, b1_mean = 'n0.25', missingness = 'mcar',DGP = 'WF', family = NULL, R = 1000, num_jobs = 50, output_path = NULL, theta = NULL, rho_DGP = NULL, alpha_DGP = NULL, tau2_DGP = NULL, rho_MAR = NULL, alpha_MAR = NULL, tau2_MAR = NULL, gamma = NULL){
  
  if(tolower(DGP) == 'mcar'){
    DGP_name = 'mcar'
  }else if(tolower(DGP) == 'car'){
    browser()
  }else if(tolower(DGP) == 'freqglm'){
    browser()
  }
  
  # make the output folder
  if(is.null(output_path)){
    output_path <- sprintf('%s%s_%s_beta0%s_beta1%s_%s', missingness, p, DGP, b0_mean, b1_mean, gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }
  
  job_name = sprintf('%s%s_%s_beta0%s_beta1%s', missingness, p, DGP, b0_mean, b1_mean)
  job_name <- gsub('\\.', '', job_name)
  
  params <- list(p = p,
                 b0_mean = b0_mean, 
                 b1_mean = b1_mean,
                 missingness = missingness,
                 DGP = DGP,
                 R = R, 
                 num_jobs = num_jobs,
                 output_path = output_path)
  
  # check that the right params are supplied
  {
    if(tolower(missingness) == 'mnar'){
      if(is.null(gamma)){
        stop('supply gamma for mnar')
      }else{
        params <- c(params, list(gamma = gamma))
      }
    } 
      
    if(tolower(missingness) == 'mar'){
      if(is.null(rho_MAR) | is.null(alpha_MAR) | is.null(tau2_MAR)){
        stop('supply params for MAR missingness')
      }else{
        params <- c(params, 
                    list(rho_MAR = rho_MAR,
                         alpha_MAR = alpha_MAR,
                         tau2_MAR = tau2_MAR))
      }
    } 
    
    if(tolower(DGP) == 'car'){
      if(is.null(rho_DGP) | is.null(alpha_DGP) | is.null(tau2_DGP)){
        stop('please input params for CAR DGP')
      }else{
        params <- c(params, 
                    list(rho_DGP = rho_DGP,
                         alpha_DGP = alpha_DGP,
                         tau2_DGP = tau2_DGP))
      }
    }else if(tolower(DGP) == 'freqglm'){
      if(is.null(rho_DGP) | is.null(alpha_DGP)){
        stop('please input params for freqGLM DGP')
      }else{
        params <- c(params, 
                    list(rho_DGP = rho_DGP,
                         alpha_DGP = alpha_DGP))
      }
    }
    
    if(!is.null(family)){
      if(family == 'quasipoisson'){
        if(is.null(theta)){
          stop('please input theta value for quasipoisson')
        }else{
          params <- c(params,
                      list(theta = theta))
        }
      }
    }
  }
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':')
  
  command_str = sprintf('sbatch --array=1-50 -J %s run_sim.sh %s', job_name, param_str)
  
  return(command_str)
}

bash_wrapper <- function(p_vec = seq(0, 0.5, 0.1), bash_file = NULL, ...){
  cmds <- data.frame(sapply(p_vec, function(xx){ bash_command(p = xx, ...)}))
  
  if(!is.null(bash_file)){
    write.table(cmds, file = bash_file, row.names = F, col.names = F, quote = F)
  }
  
  return(cmds)
}



bash_wrapper(missingness = 'mar', rho_DGP = 0.7, )