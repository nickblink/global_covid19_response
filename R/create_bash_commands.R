## This script makes bash commands for given simulations

bash_command <- function(p, b0_mean = 6, b1_mean = 'n0.25', missingness = 'mcar', DGP = 'WF', family = NULL, R = 1000, num_jobs = 50, output_path = NULL, theta = NULL, rho_DGP = NULL, alpha_DGP = NULL, tau2_DGP = NULL, rho_MAR = NULL, alpha_MAR = NULL, tau2_MAR = NULL, gamma = NULL){
  
  if(tolower(DGP) == 'wf'){
    DGP_name = 'WF'
  }else if(tolower(DGP) == 'car'){
    DGP_name = gsub('\\.', '', sprintf('CAR%s%s%s', rho_DGP, alpha_DGP, tau2_DGP))
  }else if(tolower(DGP) == 'freqglm'){
    DGP_name = gsub('\\.', '', sprintf('freqGLM%s%s', rho_DGP, alpha_DGP))
  }
  
  # make the output folder
  if(is.null(output_path)){
    output_path <- sprintf('%s%s_%s_beta0%s_beta1%s_%s', missingness, p, DGP_name, b0_mean, b1_mean, gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }
  
  job_name = sprintf('%s%s_%s_beta0%s_beta1%s', missingness, p, DGP_name, b0_mean, b1_mean)
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
                      list(theta = theta,
                           family = family))
        }
      }
    }
  }
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':')
  
  command_str = sprintf('sbatch --array=1-50 -J %s run_sim.sh %s', job_name, param_str)
  
  return(command_str)
}

bash_wrapper <- function(p_vec = seq(0, 0.5, 0.1), bash_file = NULL, ...){
  cmds <- lapply(p_vec, function(xx){ bash_command(p = xx, ...)})
  
  if(!is.null(bash_file)){
    # if it already exists, update it
    if(file.exists(bash_file)){
      lapply(cmds, write, bash_file, append = T)
    # if it doesn't exist, create it
    }else{
      out_file <- file(bash_file, open='wb')
      lapply(cmds, write, out_file, append = T)
      close(out_file)
    }
    # check that there are no repeats in commands
    test <- read.table(bash_file)
    
    if(length(unique(test[,ncol(test)])) != length(test[,ncol(test)])){
      stop('there are repeating simulation commands')
    }
  }
  
  return(cmds)
}


## commands for MAR quasipoisson
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 100, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

## commands for MNAR quasipoisson
bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

## commands for MNAR QP with WF diff. params
bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

