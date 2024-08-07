## This script makes bash commands for given simulations

bash_command <- function(p, b0_mean = 6, b1_mean = 'n0.25', missingness = 'mcar', DGP = 'WF', family = 'negbin', R = 1000, num_jobs = 50, output_path_addition = NULL, theta = NULL, DGP_theta_rate = NULL, DGP_theta_shape = NULL, rho_DGP = NULL, alpha_DGP = NULL, tau2_DGP = NULL, rho_MAR = NULL, alpha_MAR = NULL, tau2_MAR = NULL, gamma = NULL, CARburnin = 5000, CARnsample = 10000, R_PI = 200, models = c(2,4,7)){
  
  if(tolower(DGP) == 'wf'){
    DGP_name = 'WF'
  }else if(tolower(DGP) == 'car'){
    DGP_name = gsub('\\.', '', sprintf('CAR%s%s%s', rho_DGP, alpha_DGP, tau2_DGP))
  }else if(tolower(DGP) == 'freqglm'){
    if(b0_mean == 6){
      warning('are you sure you want this mean for freqglm?')
    }
    DGP_name = gsub('\\.', '', sprintf('freqGLM%s%s', rho_DGP, alpha_DGP))
  }
  
  if(!is.null(family)){
    if(family=='quasipoisson'){
      DGP_name = paste0(DGP_name, sprintf('_QPtheta%s', theta))
    }else if(!family %in% c('poisson', 'negbin')){
      stop('input a proper family type')
    }
  }else{
    family = 'poisson'
  }
  
  # make the output folder
  if(is.null(output_path_addition)){
    output_path <- sprintf('%s%s_%s_beta0%s_beta1%s_ID%s_%s', missingness, p, DGP_name, b0_mean, b1_mean, sample(1e6, size = 1), gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }else{
    output_path <- sprintf('%s%s_%s_beta0%s_beta1%s_ID%s_%s_%s', missingness, p, DGP_name, b0_mean, b1_mean, sample(1e6, size = 1), output_path_addition, gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }
  
  job_name = sprintf('%s%s_%s_beta0%s_beta1%s', missingness, p, DGP_name, b0_mean, b1_mean)
  job_name <- gsub('\\.', '', job_name)
  
  params <- list(p = p,
                 b0_mean = b0_mean, 
                 b1_mean = b1_mean,
                 missingness = missingness,
                 family = family,
                 DGP = DGP,
                 R = R, 
                 num_jobs = num_jobs,
                 output_path = output_path,
                 R_PI = R_PI,
                 models = models)
  
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
      }else if(family == 'negbin'){
        if(!is.null(DGP_theta_shape)){
          params <- c(params,
                      list(DGP_theta_shape = DGP_theta_shape))
        }
        if(!is.null(DGP_theta_rate)){
          params <- c(params,
                      list(DGP_theta_rate = DGP_theta_rate))
        }
      }
    }
  }

  # if adding in CAR burnin and nsample
  if(!is.null(CARburnin)){
    params[['CARburnin']] = CARburnin
  }
  if(!is.null(CARnsample)){
    params[['CARnsample']] = CARnsample
  }
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':') %>%
    gsub(' |c\\(|\\)','',.)
  
  command_str = sprintf('sbatch --array=1-%s -J %s run_sim.sh %s', num_jobs, job_name, param_str)
  
  return(command_str)
}

bash_wrapper <- function(p_vec = seq(0, 0.5, 0.1), bash_file = NULL, ...){
  cmds <- lapply(p_vec, function(xx){ bash_command(p = xx, ...)})
  
  if(!is.null(bash_file)){
    # if it already exists, update it
    if(file.exists(bash_file)){
      lapply(cmds, write, bash_file, append = T, sep = '')
    # if it doesn't exist, create it
    }else{
      out_file <- file(bash_file, open='wb')
      lapply(cmds, write, out_file, append = T, sep = '')
      close(out_file)
    }

    # check that there are no repeats in commands
    test <- read.table(bash_file)
    col <- lapply(test[,ncol(test)], function(str){
      tmp <- strsplit(str, ':')[[1]]
      tmp <- tmp[-grep('output_path', tmp)]
      paste(tmp, collapse = ':')
    })
    
    if(length(unique(col)) != length(col)){
      stop('there are repeating simulation commands')
    }
  }
  
  return(cmds)
}

#### Appendix simulations ####
bash_0716 <- 'cluster_code/cluster commands/bash_07162024.txt'

# All MCAR

# (1) WF beta0 = 6, beta1 = 0
bash_wrapper(DGP = 'WF', family = 'negbin', b1_mean = 0,
             bash_file = bash_0716)

# (2) freqGLM beta0 = 5.5, beta1 = 0
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, b1_mean = 0,
             bash_file = bash_0716)

# (3) CAR beta0 = 6, beta1 = 0
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, b1_mean = 0,
             bash_file = bash_0716)

# (4) WF beta0 = 2, beta1 = 0
bash_wrapper(DGP = 'WF', family = 'negbin', b1_mean = 0, b0_mean = 2,
             bash_file = bash_0716)

# (5) freqGLM beta0 = 1.5, beta1 = 0
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 1.5, b1_mean = 0,
             bash_file = bash_0716)

# (6) CAR beta0 = 2, beta1 = 0
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, b1_mean = 0, b0_mean = 2,
             bash_file = bash_0716)

# (7) CAR with higher spatial correlation (rho = 0.7 and tau2 = 1)
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.7, alpha_DGP = 0.3, tau2_DGP = 1, 
             bash_file = bash_0716)

# (8) freqGLM with higher spatial correlation (rho = 0.4, alpha = 0.2)
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.4, alpha_DGP = 0.2, b0_mean = 5.5, 
             bash_file = bash_0716)

# (9) WF with more overdispersion (higher theta rate)
bash_wrapper(DGP = 'WF', family = 'negbin', 
             DGP_theta_shape = 2.5, DGP_theta_rate = 0.667, output_path_addition = 'DGPthetarate0667',
             bash_file = bash_0716)

# (10) freqGLM with more overdispersion (higher theta rate)
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             DGP_theta_shape = 2.5, DGP_theta_rate = 0.667, output_path_addition = 'DGPthetarate0667',
             bash_file = bash_0716)

# (11) CAR with more overdispersion (higher theta rate)
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,
             DGP_theta_shape = 2.5, DGP_theta_rate = 0.667, output_path_addition = 'DGPthetarate0667',
             bash_file = bash_0716)

#### Running NB fit and DGP for MAR and MNAR - and re-doing freq ####

# (1) freqGLM MCAR
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             bash_file = bash_0716)

# (2) WF MAR
bash_wrapper(DGP = 'WF', family = 'negbin',
             missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4,
             bash_file = bash_0716)

# (3) freqGLM MAR
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4,
             bash_file = bash_0716)

# (4) CAR MAR
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, 
             missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4,
             bash_file = bash_0716)

# (5) WF MNAR
bash_wrapper(DGP = 'WF', family = 'negbin',
             missingness = 'mnar', gamma = 1, 
             bash_file = bash_0716)

# (6) freqGLM MNAR
bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             missingness = 'mnar', gamma = 1, 
             bash_file = bash_0716)

# (7) CAR MNAR
bash_wrapper(DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, 
             missingness = 'mnar', gamma = 1, 
             bash_file = bash_0716)

#### Running CAR NB DGP with all NB models ####
#bash_wrapper(missingness = 'mcar',DGP = 'CAR', family = 'negbin', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, CARburnin = 5000, CARnsample = 10000, output_path = 'mcar05_CAR0303025_beta_06_beta1n025_negbin_2024_07_11', bash_file = 'cluster_code/cluster commands/bash_07112024.txt')
DONT PUT SPECIFIC OUTPUT PATH SINCE IT WILL APPLY TO ALL p Values

#
#### Running WF NB, freqGLM NB with CAR NB fit ####
bash_wrapper(DGP = 'WF', family = 'negbin',
             bash_file = 'cluster_code/cluster commands/bash_07102024.txt')

bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             bash_file = 'cluster_code/cluster commands/bash_07102024.txt')

#
#### Running WF NB, freqGLM NB, and CAR DGP, with higher CAR burning and sample ####
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, bash_file = 'cluster_code/cluster commands/bash_07082024.txt')

bash_wrapper(DGP = 'WF', family = 'negbin',
             bash_file = 'cluster_code/cluster commands/bash_07082024.txt')

bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
             bash_file = 'cluster_code/cluster commands/bash_07082024.txt')

#### Running WF NB and freqGLM NB DGPs and CAR DGP with newer CAR method ####

bash_wrapper(DGP = 'WF', family = 'negbin',
             bash_file = 'cluster_code/cluster commands/bash_07072024.txt')

bash_wrapper(DGP = 'freqGLM', family = 'negbin', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, 
                          bash_file = 'cluster_code/cluster commands/bash_07072024.txt')

#### Creating the WF and freqGLM negbin tests ####
bash_command(p = 0.1, DGP = 'WF', family = 'negbin', R = 100, num_jobs = 5)

bash_command(p = 0.1, DGP = 'freqGLM', family = 'negbin', R = 100, num_jobs = 5, rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5)

#
#### 7/03/2024: Comparing CAR prediction on freqGLM and CAR DGP ####
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_07032024.txt', models = c(5,6))

bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_07032024.txt', models = c(5,6))
#
#### 6/24/2024: Running WF with MNAR and MAR - WITHOUT quasipoisson ####
#	CAR0303025 EB0 = 6, EB1 = -0.25: MAR
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, DGP = 'WF', b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06242024.txt')

#	CAR0303025 EB0 = 6, EB1 = -0.25: MNAR
bash_wrapper(missingness = 'mnar', gamma = 1, DGP = 'WF', b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06242024.txt')

#
#### 6/19/2024: Running all DGP with higher CAR nsample and running freqGLM DGP with MAR, freqGLM DGP with MNAR, CAR DGP with MAR, CAR DGP with MNAR ####

#	CAR0303025 EB0 = 6, EB1 = -0.25: MCAR
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	freqGLM EB0 = 5.5, EB1 = -0.25: MCAR
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	WF EB0 = 6, EB1 = -0.25: MCAR
bash_wrapper(missingness = 'mcar', b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	freqGLM EB0 = 5.5, EB1 = -0.25: MAR
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	freqGLM EB0 = 5.5, EB1 = -0.25: MNAR
bash_wrapper(missingness = 'mnar', gamma = 1, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	CAR0303025 EB0 = 6, EB1 = -0.25: MAR
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#	CAR0303025 EB0 = 6, EB1 = -0.25: MNAR
bash_wrapper(missingness = 'mnar', gamma = 1, DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06192024.txt')

#
#### 6/17/2024: Testing higher n.sample and burnin ####
#	CAR0303025 EB0 = 6, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, b1_mean = 0, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_06172024.txt')

bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_v2_06172024.txt')

#	WF EB0 = 6, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar', b0_mean = 6, b1_mean = -0.25, CARburnin = 5000, CARnsample = 10000, bash_file = 'cluster_code/cluster commands/bash_v2_06172024.txt')

#
#### 6/17/2024: Testing newer code ####
bash_command(p = 0.2, missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, R = 10, num_jobs = 5)


#### 12/14/2023: Results for the appendix! ####
# freqGLM0202 EB0 = 5.5, EB1 = -0.25: MCAR (not actually the appendix though)
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	freqGLM0202 EB0 = 5.5, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5.5, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 6, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar', b0_mean = 6, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	CAR0303025 EB0 = 6, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 6, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	freqGLM0202 EB0 = 1.5, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 1.5, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 2, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar', b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')


#	CAR0303025 EB0 = 2, EB1 = 0: MCAR
bash_wrapper(missingness = 'mcar',DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25,  b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 6, EB1 = 0: QP 4: MCAR
bash_wrapper(missingness = 'mcar', b1_mean = 0, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 6, EB1 = 0: QP 4: MAR
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

# WF EB0 = 6, EB1 = 0: QP 4: MNAR
bash_wrapper(missingness = 'mnar', b1_mean = 0, gamma = 1, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 6, EB1 = 0: QP 16: MCAR
bash_wrapper(missingness = 'mcar', b1_mean = 0, family = 'quasipoisson', theta = 16, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

#	WF EB0 = 6, EB1 = 0: QP 16: MAR
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, family = 'quasipoisson', theta = 16, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')

# WF EB0 = 6, EB1 = 0: QP 16: MNAR
bash_wrapper(missingness = 'mnar', b1_mean = 0, gamma = 1, family = 'quasipoisson', theta = 16, bash_file = 'cluster_code/cluster commands/bash_12142023.txt')



#### 12/10/2023: Take 2: Commands for results section of paper (not appendix) ####
# QP theta = 4 WF B0 = 6, b1 = -0.25 MCAR
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12102023_p2.txt')

# QP theta = 4 WF B0 = 6, b1 = -0.25 MAR
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12102023_p2.txt')

# QP theta = 4 WF B0 = 6, b1 = -0.25 MNAR
bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12102023_p2.txt')

#
#### 12/10/2023: Commands for appendix of paper ####
# First I will check that the results for the paper are as expected, but then I will do commands for:

MAKE THESE IN ORDER OF IMPORTANCE - LIKE IF ANY WILL BE USED FOR MY SUBMISSION TO BETHANY.

# WF, CAR, freqGLM DGP with beta = 2, 0
# Missing comparison MCAR, MAR, MNAR, with theta = 16
# Missing comparison MCAR, MAR, MNAR with Poisson variance
# Missing comparison MCAR, MAR, MNAR with QP theta = 4, beta = 6, 0


#### 12/9/2023: Commands for results section of paper (not appendix) ####
# freqGLM B0 = 6, b1 = -0.25, MCAR params 0202
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 5, bash_file = 'cluster_code/cluster commands/bash_12092023_p2.txt')

# QP theta = 4 WF B0 = 6, b1 = -0.25 MCAR
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12092023_p2.txt')

# QP theta = 4 WF B0 = 6, b1 = -0.25 MAR
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 4, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12092023_p2.txt')

# QP theta = 4 WF B0 = 6, b1 = -0.25 MNAR
bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 4, bash_file = 'cluster_code/cluster commands/bash_12092023_p2.txt')

#### On 12/8/2023, testing ish now, because ish is wack ####
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12092023.txt')
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12092023.txt')
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, bash_file = 'cluster_code/cluster commands/bash_12092023.txt')
#
#### On 12/5/2023 part 2 #####
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 100, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

## commands for MNAR quasipoisson
bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

## commands for MNAR QP with WF diff. params
bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MCAR QP9
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MCAR QP100
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MCAR QP9 EB1 = 0
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 9, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MCAR QP100 EB1 = 0
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 100, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MAR QP9 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MAR QP100 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# CAR MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# freqGLM MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')


#
#### On 12/5/2023 ####
# fixing the MAR QP B1_0 runs. Idk what happened with those.

# WF MAR QP9 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

# WF MAR QP100 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12052023.txt')

#
#### On 12/4/2023 ####
# filling in the two runs that were messed up last time
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12042023_fixed.txt')

bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12042023_fixed.txt')
# ^manually delete the ones I dont need

# WF MCAR QP9
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MCAR QP100
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MCAR QP9 EB1 = 0
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 9, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MCAR QP100 EB1 = 0
bash_wrapper(missingness = 'mcar', family = 'quasipoisson', theta = 100, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MAR QP9 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MAR QP100 EB1 = 0
bash_wrapper(missingness = 'mar', b1_mean = 0, rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# WF MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# CAR MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', DGP = 'CAR', rho_DGP = 0.3, alpha_DGP = 0.3, tau2_DGP = 0.25, b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

# freqGLM MCAR low mean EB0 = 2
bash_wrapper(missingness = 'mcar', DGP = 'freqglm', rho_DGP = 0.2, alpha_DGP = 0.2, b0_mean = 2, b1_mean = 0, bash_file = 'cluster_code/cluster commands/bash_12042023.txt')

#
#### On 11/30/2023 ####
## commands for MAR quasipoisson
bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 9, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(missingness = 'mar', rho_MAR = 0.7, alpha_MAR = 0.7, tau2_MAR = 100, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

## commands for MNAR quasipoisson
bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

## commands for MNAR QP with WF diff. params
bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 9, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

bash_wrapper(b1_mean = 0, missingness = 'mnar', gamma = 1, family = 'quasipoisson', theta = 100, bash_file = 'cluster_code/cluster commands/bash_11302023.txt')

