## This script makes bash commands for given simulations

TO DO:
  COMPLETE THE PARAMS input list in function
  complete the output path namer
  Complete the string creation
  make wrapper function to cycle through p and to write a file

sbatch --array=1-50 -J MNAR_WF_TEST run_sim.sh p=0.3:b0_mean=6:b1_mean=n0.25:missingness=mnar:gamma=1:DGP=WF:family=quasipoisson:theta=9:R=1000:num_jobs=50:output_path=mnar03_WF_QP9_TEST_beta6_n025_2023_11_28

bash_command <- function(p = NULL, b0_mean = 6, b1_mean = 'n0.25', missingness = 'mcar',DGP = 'WF', family = NULL, R = 1000, num_jobs = 50, output_path = NULL, theta = NULL, rho_DGP = NULL, alpha_DGP = NULL, tau2_DGP = NULL, rho_MAR = NULL, alpha_MAR = NULL, tau2_MAR = NULL, gamma = NULL){
  
  # check that the right params are supplied
  {
  if(tolower(missingness) == 'mnar' & is.null(gamma)){
    stop('supply gamma for mnar')
  }else if(towlower(missingness) == 'mar' & (is.null(rho_MAR) | is.null(alpha_MAR) | is.null(tau2_MAR))){
    stop('supply params for MAR missingness')
  }
  
  if(tolower(DGP) == 'car'){
    if(is.null(rho_DGP) | is.null(alpha_DGP) | is.null(tau2_DGP)){
      stop('please input params for CAR DGP')
    }
  }else if(tolower(DGP) == 'freqglm'){
    if(is.null(rho_DGP) | is.null(alpha_DGP)){
      stop('please input params for freqGLM DGP')
    }
  }
  
  if(!is.null(family)){
    if(family == 'quasipoisson' & is.null(theta)){
      stop('please input theta value for quasipoisson')
    }
  }
  }
  
  params <- list(p = p,
                 b0_mean = b0_mean, 
                 b1_mean = b1_mean,
                 missingness = missingness,
                 DGP)
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':')
}