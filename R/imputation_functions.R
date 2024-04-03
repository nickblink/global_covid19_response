### The functions used for imputation 
options(dplyr.summarise.inform = FALSE)

##### Helper Functions #####
# pull out the p value from a file name
get_p_from_name <- function(file){
  p <- c(stringr::str_match(file, 'mcar(.*?)_')[[2]],
         stringr::str_match(file, 'mnar(.*?)_')[[2]],
         stringr::str_match(file, 'mar(.*?)_')[[2]])
  p <- p[!is.na(p)]
  return(p)
}

rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

# make the adjacency matrix according to all facilities in a district being neighbors
make_district_adjacency <- function(df, scale_by_num_neighbors = F){
  
  # get the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% 
    distinct() %>%
    arrange(facility)
  W = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # scale the neighbor sum by the number of neighbors
  if(scale_by_num_neighbors){
    W = apply(W, 2, function(xx) xx/sum(xx))
  }
  
  return(W)
}

# The "W2" matrix, as I am calling it, is diag(W1) - W. Aka the negative adjacency matrix with the total number of neighbors for each facility on the diagonal
make_district_W2_matrix <- function(df){
  # create the list of matching facilities
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  pairs = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district)
  
  # get unique facilities
  facilities = unique(df$facility) %>% sort()
  
  # initialize the W2 matrix (note that what I am calling W2 here is what the original paper calls diag(W1) - W)
  W2 = matrix(0, nrow = length(facilities), ncol = length(facilities))
  colnames(W2) = facilities
  rownames(W2) = facilities
  
  for(i in 1:length(facilities)){
    f1 = facilities[i]
    
    # get the pairs with this facility
    tmp = pairs %>% filter(facility.x == f1)
    
    if(nrow(tmp) == 0){
      print('havent done single-district facilities yet')
      browser()
    }else{
      # put -1 where there are pairs of facilities
      matid = match(tmp$facility.y, facilities)
      W2[i,matid] = -1
    }
    
    # get the number of neighbors this facility has
    W2[i,i] = pairs %>% filter(facility.x == f1) %>% nrow()
    # print(sprintf('%s: NUM neighbors = %s', f1, W2[i,i]))
  }
  
  return(W2)
}

# function to add periodic terms
add_periodic_cov <- function(df, period = 12){
  df = df %>%
    dplyr::mutate(year = lubridate::year(date) - min(lubridate::year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  return(df)
}

# add autoregressive terms to the data
add_autoregressive <- function(df, target_col = 'y', num.terms = 1){
  
  if(!(target_col %in% colnames(df))){
    stop('need a target column to autoregress')
  }
  
  tmp <- lapply(unique(df$facility), function(xx) {
    tt <- df %>% filter(facility == xx) %>% arrange(date)
    #tt[,sprintf('%s.AR1', target_col)] = c(NA, tt[1:(nrow(tt) - 1), target_col, drop = T])
    tt[,'y.AR1'] = c(NA, tt[1:(nrow(tt) - 1), target_col, drop = T])
    return(tt)
  })
  
  df <- do.call('rbind',tmp)
  return(df)
}

# compute the sum of the values at all the neighbors to a given facility
add_neighbors <- function(df, target_col = 'y', lag = 1, scale_by_num_neighbors = F, W = NULL){
  if(lag == 0){
    print('WARNING: not doing a lag of 1 on the neighbors, so not using the same model as in the papers')
  }

  # remove the neighbor column
  if('y.neighbors' %in% colnames(df)){
    df$y.neighbors = NULL
  }
  
  # get the adjacency matrix
  if(is.null(W)){
    W <- make_district_adjacency(df, scale_by_num_neighbors)
  }

  # get the counts for each facility 
  y.counts <- df %>% 
    dplyr::select(date, facility, UQ(target_col)) %>%
    arrange(facility) %>%
    tidyr::spread(facility,get(target_col)) %>% 
    arrange(date)
  
  # check that the column names match
  if(!identical(colnames(W), colnames(y.counts)[-1])){
    # if the colnames do match but are in the wrong order, reorder them
    if(length(setdiff(colnames(W), colnames(y.counts)[-1])) == 0 &
       length(setdiff(colnames(y.counts)[-1], colnames(W))) == 0){
      
      matid = match(colnames(y.counts)[-1], colnames(W))
      y.counts = y.counts[,c(1, matid + 1)]
      
      # for error checking
      if(!identical(colnames(W), colnames(y.counts)[-1])){browser()}
    }else{
      stop('Adjacency and y matrix column names dont match')
    }
    ind = which(colnames(W) != colnames(y.counts)[-1])
  }
  
  # shift y.counts by lag
  if(lag > 0){
    tmp <- y.counts
    tmp[1:lag,-1] <- NA
    tmp[(lag+1):nrow(tmp),-1] <- y.counts[1:(nrow(y.counts) - lag),-1]
    y.counts <- tmp
  }
  
  # merge back into original data frame
  tmp = cbind(y.counts[,'date',drop = F], as.data.frame(as.matrix(y.counts[,-1])%*%W)) %>% 
    tidyr::gather(facility, y.neighbors, -date)
  if(is.factor(df$facility)){
    tmp$facility = factor(tmp$facility, levels = levels(df$facility))  
  }
  df = merge(df, tmp, by = c('date','facility'))
  
  return(df)
}

# prep the data for fitting the stan model
prep_stan_data_rushworth_sparse <- function(df, formula){
  N = nrow(df)
  N_T <- length(unique(df$date))
  N_F <- length(unique(df$facility))
  
  # order the data frame
  df <- df %>% arrange(date, facility)
  
  # W_star = make_district_W2_matrix(df)
  W = make_district_adjacency(df)
  W_n = sum(W)/2
  W2 = make_district_W2_matrix(df)
  
  # get eigenvalues for determinant calculation
  lambda = eigen(W2 - diag(1, nrow(W2)))$values
  
  # make the complete model matrix
  df2 <- df; df2$y[is.na(df2$y)] <- 0
  X = model.matrix(formula, data = df2)
  
  # get the outcome
  y = df$y
  
  # # comparing missingness.
  # if(!identical(as.integer(rownames(X_obs)), which(!is.na(df$y)))){
  #   stop('mismatch of model matrix and df missing rows')
  # }
  
  # missingness data
  N_miss = sum(is.na(y))
  N_obs = sum(!is.na(y))
  ind_miss = which(is.na(y))
  ind_obs = which(!is.na(y))
  y_obs = y[ind_obs]
  
  # compute the priors
  lm_fit <- glm(formula, family = 'poisson', data = df)
  coef_mat <- summary(lm_fit)$coefficients
  prior_mean_beta <- coef_mat[,1]
  sigma_beta = 10*vcov(lm_fit)
  
  # make the stan data frame
  stan_data <- list(
    N = N, # number of observations
    p = ncol(X), # number of variables
    N_F = N_F, # number of facilities
    N_T = N_T, # number of time points
    N_miss = N_miss,
    N_obs = N_obs,
    ind_miss = ind_miss,
    ind_obs = ind_obs,
    X = X, # design matrix
    y_obs = y_obs, # outcome variable 
    mu = prior_mean_beta, # prior mean
    Sigma = sigma_beta, # prior variance
    #W_star = W_star, 
    W = W,
    W_n = W_n,
    I = diag(1.0, N_F),
    lambda = lambda)
  
  return(stan_data)
}

# get the district and facility list from a data frame
get_district_facilities <- function(df){
  tt = unique(df[,c('district','facility')])
  res_lst <- NULL
  for(i in 1:nrow(tt)){
    res_lst[[tt[i,1]]] <- c(res_lst[[as.character(tt[i,1])]], as.character(tt[i,2]))
  }
  return(res_lst)
}

# check that the dates and facilities align in a list of data frames
date_facility_check <- function(imputed_list){
  date_fac_mat = imputed_list[[1]][, c("date", "facility")]
  date_fac = paste0(date_fac_mat[,1], '--', date_fac_mat[,2])
  
  for(i in 2:length(imputed_list)){
    date_fac_mat_new = imputed_list[[i]][, c("date", "facility")]
    date_fac_new = paste0(date_fac_mat_new[,1], '--', date_fac_mat[,2])
    if(!identical(date_fac_new, date_fac)){
      browser()
      stop('mismatch in dates and facilities in the imputed list run')
    }
  }
}

# check that the dates and districts align in a list of data frames
date_district_check <- function(imputed_list){
  date_fac_mat = imputed_list[[1]][, c("date", "district")]
  date_fac = paste0(date_fac_mat[,1], '--', date_fac_mat[,2])
  
  for(i in 2:length(imputed_list)){
    date_fac_mat_new = imputed_list[[i]][, c("date", "district")]
    date_fac_new = paste0(date_fac_mat_new[,1], '--', date_fac_mat[,2])
    if(!identical(date_fac_new, date_fac)){
      stop('mismatch in dates and districts in the imputed list run')
    }
  }
}

# fit the basic WF model
fit_WF_model <- function(data, outcome = 'indicator_count_ari_total', facilities = NULL, family = 'nb', max_date = '2019-12-01'){
  
  print(sprintf('filtering data to be less than %s', max_date))
  data = data %>%
    filter(date <= max_date)
  
  # use all facilities if not provided
  if(is.null(facilities)){
    facilities <- unique(data$facility)
  }
  
  data <- add_periodic_cov(data)
  data$y <- as.data.frame(data)[,outcome]

  formula_col = as.formula(sprintf("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"))
  
  print(sprintf('fitting models for %s family', family))
  if(family %in% c('nb','NB','negative binomial')){
    tmp <- lapply(facilities, function(xx) {
      tt <- data %>% filter(facility == xx)
      
      # run the model
      mod_col <- tryCatch({MASS::glm.nb(formula_col, data = tt)},
                          error = function(e){
                            NA
                          })
      
      return(mod_col)
    })
  }else{
    stop('havent coded for other families')
  }
  
  # combine the results
  df = do.call('rbind', lapply(tmp, function(xx) {
    if(all(is.na(xx))){
      return(rep(NA, 8))
    }else{
      if(!xx$converged){
        return(rep(NA, 8))
      }else{
        return(xx$coefficients)
      }
    }
  }))
    
  df <- as.data.frame(df) %>%
    mutate(facility = facilities)
  
  facility_miss = as.data.frame(data) %>% 
    group_by(facility) %>%
    summarize(num_missing = sum(is.na(y)))
  
  df <- merge(facility_miss, df) %>%
    arrange(num_missing)
  
  return(df)
}

# Check names from methods and output list
method_name_check <- function(lst, methods, give_method_name_err = T){
  err = F
  ind_rm = c()
  for(m in methods){
    chk <- sapply(lst, function(xx) any(grepl(m, colnames(xx))))
    if(any(chk == FALSE)){
      err = T
      print(sprintf('The method %s is not found in %s of the results list', m, sum(!chk)))
      print(head(which(chk == F)))
      ind_rm <- c(ind_rm, which(chk == F))
    }
  }
  if(err){
    if(give_method_name_err){
      print(colnames(lst[[1]]))
      stop('incorrect methods specified')
    }else{
      lst <- lst[-ind_rm]
      lst <- method_name_check(lst, methods, give_method_name_err = T)
    }
  }
  return(lst)
}

### Check that the names of the simulation outcomes are there. This is to avoid running a super long simulation and then find out I didn't return all the results
outcome_name_checker <- function(res, models = c('WF','CAR','freqGLM'), quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)){
  cols <- c('district','date','y','y_true','y_exp','y_var')
  for(m in models){
    cols <- c(cols, paste0('y_pred_',m,'_',quant_probs))
  }
  
  missing_cols = setdiff(c('facility', cols), colnames(res$df_miss))
  missing_dist_cols = setdiff(cols, colnames(res$district_df))
  
  if(length(missing_cols) > 0){
    print(sprintf('missing cols from facility df: %s', paste0(missing_cols, collapse = ', ')))
  }
  
  if(length(missing_dist_cols) > 0){
    print(sprintf('missing cols from district df: %s', paste0(missing_dist_cols, collapse = ', ')))
  }
}

### combine saved results from batch runs into one file
combine_results <- function(input_folder, results_file = NULL, return_lst = T, return_raw_list = F, ignore_size_err = F, district_results = F, subset_results = NULL){
  files <- dir(input_folder, full.names = T)
  files <- grep('sim_results', files, value = T)
  
  load(files[1])

  lst_full <- imputed_list
  params_list <- list()
  tmp <- list(params); names(tmp) = files[1]
  params_list <- c(params_list, tmp)

  rm(imputed_list)
  for(i in 2:length(files)){
    load(files[i])
    lst_full <- c(lst_full, imputed_list)
    tmp <- list(params); names(tmp) = files[i]
    params_list <- c(params_list, tmp)
  }
  
  if(length(lst_full) != 1000){
    stop(sprintf('lst_full only has %s simulations. There should be 1000.', length(lst_full)))
  }
  
  # combine the parameters
  param_mat = data.table::rbindlist(params_list, idcol = 'file')
  
  if(!is.null(subset_results)){
    lst_full = lst_full[subset_results]
  }
  
  if(return_raw_list){
    return(lst_full)
  }
  
  if(class(lst_full[[1]]) == 'data.frame'){
    df_lst <- lst_full
    WF_lst <- error_lst <- true_betas <- NULL
  }else{
    df_lst <- lapply(lst_full, function(tmp) {tmp$df_miss})
    #if('WF_betas' %in% names(lst_full[[1]])){
    WF_lst <- lapply(lst_full, function(tmp) {tmp$WF_betas})
    #}
    CARstan_summary <- lapply(lst_full, function(tmp) {tmp$CARstan_summary})
    
    freqGLM_lst <- lapply(lst_full, function(tmp) {tmp$freqGLM_params})
    
    if(district_results){
      district_lst <- lapply(lst_full, function(tmp) {tmp$district_df})
    }
    
    models = names(lst_full[[1]][[2]])
    error_lst <- NULL
    for(m in models){
      error_lst[[m]] <- do.call('rbind',lapply(lst_full, function(tmp) tmp$errors[[m]]))
    }
  }

  # checking the size of the objects - can tell if there's an error
  if(district_results == T){
    check_lst <- district_lst
  }else{
    check_lst <- df_lst
  }
  object_sizes <- c()
  for(i in 1:length(check_lst)){
    obj <- object.size(check_lst[[i]])
    object_sizes <- c(object_sizes, obj)
    if(obj == 0 & !ignore_size_err){
      stop(sprintf('object size is 0 at list value %i', i))
    }
  }
  if(any(object_sizes < .9*mean(object_sizes))){
    ind <- which(object_sizes < .9*mean(object_sizes))
    warning(sprintf('there are %s out of %s object sizes less than half the mean size. That shouldnt be. Removing them', length(ind), length(check_lst)))
    print(head(ind))
    if((length(ind) > 0.05*length(check_lst)) & !ignore_size_err){
      stop('too many files of incomplete size.')
    }else{
      check_lst = check_lst[-ind]
    }
  }
  if(length(files) != 50 | length(check_lst) != length(lst_full)){
    print(sprintf('num files: %i. %s out of %s simulation results kept', length(files), length(check_lst), length(lst_full)))
  }
  
  if(district_results == T){
    district_lst <- check_lst
  }else{
    df_lst <- check_lst
  }
  
  # store the results
  res_lst <- list(error_lst = error_lst, 
                  WF_lst = WF_lst,
                  CARstan_summary = CARstan_summary,
                  freqGLM_lst = freqGLM_lst,
                  params = param_mat)
  
  if(district_results){
    res_lst[['district_lst']] = district_lst 
  }else{
    res_lst[['df_lst']] = df_lst
  }
  
  # doing this because deprecated results didnt have true betas
  if(exists('true_betas')){
    res_lst$true_betas = true_betas
  }else{
    print('no true betas found in these results')
  }
  
  if(!is.null(results_file)){
    save(res_lst, file = results_file)
  }
  
  if(return_lst){
    return(res_lst)
  }
}

### Takes in a list of files and/or directories. For each of these, pulls in the data using "combine_results" and then computes the metrics for these results together.
combine_results_wrapper <- function(files, district_results = F, methods = c("y_pred_CCA_WF", "y_pred_CCA_CAR", "y_pred_CCA_freqGLMepi"), rename_vec = c('WF','CAR','freqGLM'), metrics = c('bias', 'relative_bias', 'RMSE', 'coverage95','specificity', 'interval_width','outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'),  results_by_point = F, give_method_name_err = T, return_unprocessed = F,  QP_variance_adjustment = NULL, ...){ 
  
  # initialize results catchers
  res <- NULL
  params <- NULL
  
  if(district_results){
    print('getting district level results')
    output = 'district_lst'
  }else{
    print('getting facility level results')
    output = 'df_lst'
  }
  
  for(file in files){
    p <- get_p_from_name(file)
    if(length(p) == 0){
      print('cant find p. Artificially set equal to 0')
      p = 0
    }
    # }else{
    #   print(sprintf('p = %s', p))
    # }
    
    print(file)
    
    # if a directory, combine results. If already combined, load them
    if(dir.exists(file)){
      lst_full <- combine_results(input_folder = file, return_lst = T, results_file = NULL, district_results = district_results, ...)#, subset_results = subset_results, ignore_size)
    }else{
      load(file)
    }
    tmp = lst_full$params
    tmp$dir = file
    params = tryCatch(rbind(params, tmp),
             error = function(e){
               print(sprintf('ERROR WITH PARAMS IN FILE %s', file))
               params = plyr::rbind.fill(params, tmp)
               return(params)
             })
    #browser()
    
    # Check that the names match
    lst_full[[output]] <- method_name_check(lst_full[[output]], methods, give_method_name_err = give_method_name_err)
    
    if(return_unprocessed){
      res[[p]] <- lst_full
      next
    }
  
    if(district_results){
      date_district_check(lst_full$district_lst)
    }else{
      # make sure the dates and facilities match
      date_facility_check(lst_full$df_lst)
    }
    
    if(!any(grepl('outbreak_detection', metrics))){
      warning('no outbreak detection. Artificially creating the y_exp value')
      lst_full[[output]] <- lapply(lst_full[[output]], function(df){
        df$y_exp = df$y_true + 1
        df
      })
    }

    # Calculate the metrics
    tmp <- calculate_metrics(lst_full[[output]], results_by_point = results_by_point, methods = methods, imputed_only = F, rm_ARna = F, use_point_est = F, date = '2020-01-01', district_results = district_results, QP_variance_adjustment = QP_variance_adjustment) 
    #tmp$method = paste0(tmp$method, sprintf('_p%s_', p))
    tmp$prop_missing = as.numeric(p)/10
    
    #tmp$b0_mean <- paste(as.character(lst_full$params$b0_mean), collapse = '/')
    #tmp$b1_mean <- paste(as.character(lst_full$params$b1_mean), collapse = '/')
    
    res <- rbind(res, tmp)
  }
  
  # replace the names
  if(!is.null(rename_vec)){
    for(i in 1:length(methods)){
      res$method = gsub(methods[i], rename_vec[i], res$method)
    }
    
    # make a factor so the ordering in the plots stays
    res$method =  factor(res$method, levels = rename_vec)
  }
  
  res_lst = list(results = res, params = params)
  return(res_lst)
}

# pulls out the CAR parameter estimates from an imputed list and organizes them
process_CAR_params <- function(imputed_list, all_betas = F, CAR_name = 'CAR_summary'){
  # get row indices of betas
  tt <- imputed_list[[1]][[CAR_name]]
  ind = which(!(rownames(tt) %in% c('tau2','rho.S','rho.T')))
  
  # get all row names
  rn = rownames(tt)
  
  # Organize betas and pull out mean and n.effective
  tmp <- lapply(imputed_list, function(xx){
    tt = xx[[CAR_name]]
    if(!identical(rownames(tt), rn)){
      stop('unmatched row names')
    }
    tt2 <- rbind(tt, matrix(colMeans(tt[ind,]), 
                            nrow = 1, 
                            dimnames = list('betas', NULL)))
    return(list(tt2[,1], tt2[,'n.effective'], tt2[,'% accept']))
  })
  
  # create matrix of means and n effective numbers
  param_means <- do.call('cbind', lapply(tmp, '[[', 1))
  param_n <- do.call('cbind', lapply(tmp, '[[', 2))
  param_accept <- do.call('cbind', lapply(tmp, '[[', 3))
  
  # combine into one df
  df <- data.frame(param = rownames(param_means),
                   mean_acceptance = apply(param_accept, 1, mean),
                   mean_param = apply(param_means, 1, mean),
                   median_param = apply(param_means, 1, median),
                   param_05 = apply(param_means, 1, function(x){ quantile(x, probs = 0.05)}),
                   param_95 = apply(param_means, 1, function(x){ quantile(x, probs = 0.95)}),
                   mean_n_eff = apply(param_n, 1, mean),
                   median_n_eff = apply(param_n, 1, median),
                   n_eff_05 = apply(param_n, 1, function(x){ quantile(x, probs = 0.05)}),
                   n_eff_95 = apply(param_n, 1, function(x){ quantile(x, probs = 0.95)}))
  
  if(!all_betas){
    df <- df %>% 
      filter(param %in% c('betas','tau2', 'rho.S','rho.T'))
  }
  
  return(df)
}

# pulls out CAR estimates from a set of files where results are located and splits up results by the proportion missing
process_CAR_params_wrapper <- function(files, rename_params = T, all_betas = F, use_prop_missing = F, CAR_names = c('CAR_summary')){
  # expected = data.frame(param = c('tau2','rho.S','alpha'), expected = c(1, 0.3,0.3))
  res_df <- NULL
  for(f in files){
    lst <- combine_results(f, return_raw_list = T)
    for(CN in CAR_names){
      res <- process_CAR_params(lst, all_betas, CAR_name = CN)
      if(use_prop_missing){
        p <- get_p_from_name(f)
        res$prop_missing = as.numeric(p)/10
      }
      res$CAR_name = CN
      res_df <- rbind(res_df, res)
    }
  }
  
  if(rename_params){
    res_df$param <- gsub('rho.T','alpha', res_df$param)
  }
  
  return(res_df)
}

#
##### Imputation Functions #####

# Weingberger-Fulcher imputation method
WF_imputation <- function(df, col, group = 'facility', family = 'NB', period = 12, R_PI = 500, bias_correction_chen = F, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)){
  
  # print(sprintf('family is %s', family))
  
  # check if this method has already been run
  if(any(grepl('y_pred_WF', colnames(df)))){
    print('previous WF predictions found. Removing them')
    df[,grep('y_pred_WF', colnames(df))] <- NULL
  }
  
  # prep the data with the harmonic functions
  df <- add_periodic_cov(df, period = period)
  
  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  
  if(family == 'quasipoisson'){
    warning('havent coded PIs in quasipoisson or NB yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=quasipoisson)
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    #df <- data.table::rbindlist(tmp)
    df <- do.call('rbind',tmp)
    
  }else if(family %in% c('NB','negative binomial','neg_bin')){
    warning('havent coded PIs in quasipoisson or NB yet')
    
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      #mod_col <- MASS::glm.nb(formula_col, data = tt, control = glm.control(maxit=200,trace = 3))
      
      tt[,paste0(col, '_pred_WF')] <- tryCatch({
        mod_col <- MASS::glm.nb(formula_col, data = tt)
        predict(mod_col, tt, type = 'response')
      }, error = function(e){
        print(sprintf('glm.nb failed for %s', xx))
        rep(NA, nrow(tt))
      })
      #tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      return(tt)
    })
    
    df <- data.table::rbindlist(tmp)
  }else if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(get(group) == xx)
      mod_col <- glm(formula_col, data = tt, family=poisson)
      
      tt[,paste0(col, '_pred_WF')] = predict(mod_col, tt, type = 'response')
      
      if(R_PI > 0){
        # extract information from model fit
        beta_hat <- mod_col$coefficients
        beta_vcov <- vcov(mod_col)

        # if correcting the bias, compute and correct it
        if(bias_correction_chen){
          X = tt %>%
            mutate(intercept = 1) %>%
            dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3) %>%
            as.matrix()
          
          # the not good lin alg way
          tmp = lapply(1:nrow(X), function(i){
            x = X[i,]
            res = tt$y_exp[i]*x%*%t(x)
            return(res)
          })
          
          H1 = -Reduce("+", tmp) / length(tmp)
          Q = solve(H1)
          
          H2 = lapply(1:length(beta_hat), function(p){
            tmp = lapply(1:nrow(X), function(i){
              x = X[i,]
              res = x[p]*tt$y_exp[i]*x%*%t(x)
              return(res)
            })
            
            val = -Reduce("+", tmp) / length(tmp)
            return(val)
          })
          
          H2 = do.call('rbind', H2)
          bias_analytical = 1/(2*nrow(tt))*Q%*%t(H2)%*%c(Q)
          
          # update the beta hats to s
          beta_hat <- beta_hat - bias_analytical  
        }
        
        # store the model results to return
        model_res = list()
        model_res[[as.character(xx)]][['beta_hat']] <- beta_hat
        model_res[[as.character(xx)]][['beta_vcov']] <- beta_vcov
        
        # bootstrap 
        sapply(1:R_PI, function(r){
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (tt %>% 
                          mutate(intercept=1) %>%
                          #dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          dplyr::select(intercept,year,cos1,sin1,cos2,sin2,cos3,sin3) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot) 
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
          x = rpois(n = nrow(tt), pred_boot_exp)
          
          return(x)
          
        }) -> sim.boot
        
        # get the quantiles and store them
        fitted_quants = t(apply(sim.boot, 1, function(xx){
          quantile(xx, probs = quant_probs)
        }))
        fitted_quants = as.data.frame(fitted_quants)
        colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
        
        # merge the quantiles back into data frame
        tt = cbind(tt, fitted_quants)
        
      }
      tmp_lst = list(tt, sim.boot, model_res)
      return(tmp_lst)
    })
    names(tmp) = uni_group
    
    # initialize district results
    district_df = NULL
    
    # district-level analysis
    for(d in names(dist_fac)){
      tt = data.frame(district = d,
                      date = tmp[[1]][[1]]$date)
      facs = dist_fac[[d]]
      
      # get the sums by district: returns n x R data frame
      sum_district = Reduce('+', lapply(facs, function(f){
        tmp[[f]][[2]]
      }))
      
      # get the quantiles and store them
      fitted_quants = t(apply(sum_district, 1, function(xx){
        quantile(xx, probs = quant_probs)
      }))
      fitted_quants = as.data.frame(fitted_quants)
      colnames(fitted_quants) = paste0(paste0(col, '_pred_WF_'), quant_probs)
      
      # merge the quantiles back into data frame
      tt = cbind(tt, fitted_quants)
      
      district_df = rbind(district_df, tt)
    }
    
    betas <- do.call('rbind',lapply(1:length(tmp), function(ii){
      tmp_fac <- tmp[[ii]][[3]]
      tt <- matrix(NA, ncol = 8)
      tt[1,] <- tmp_fac[[1]][[1]]
      rownames(tt) <- names(tmp_fac)
      colnames(tt) <- names(tmp_fac[[1]][[1]])
      return(tt)
    }))
    
    beta_vcovs <- lapply(1:length(tmp), function(ii){
      return(tmp[[ii]][[3]])
    })
    
    # combine the individual facility results into a larger data frame
    df <- do.call('rbind',lapply(tmp,'[[',1))
    
    # take the sum of the bootstrap samples of each facility (this returns an n X R matrix itself)
    sim.full = Reduce('+', lapply(tmp, '[[', 2))

  }
  res_lst = list(df = df, district_df = district_df, betas = betas, beta_vcovs = beta_vcovs)
  return(res_lst)
}

# run a complete case analysis for the WF method
WF_CCA <- function(df, district_df = NULL, train_end_date = '2019-12-01', ...){
  # replace values in the test set with missing ones
  tmp <- df
  tmp$y[tmp$date > train_end_date] <- NA
  
  res <- WF_imputation(tmp, ...)

  # store the results and return the original y values
  # This is because the y values in 2020 are returned as NA from WF_imputation
  tmp <- res$df
  tmp = merge(tmp %>% dplyr::select(-y),
              df %>% dplyr::select(date, facility, y),
              by = c('date','facility'))
  
  res$df <- tmp
  
  return(res)
}

# Run Weinberger-Fulcher (WF) model using all observed data with no missingness as a baseline comparison.
WF_baseline <- function(df, train_end_date = '2019-12-01', family = 'poisson', R_PI = 100){
  # replace missing points with their true values
  tmp <- df
  tmp$y <- tmp$y_true
  tmp$y[tmp$date > train_end_date] <- NA
  
  # temporary store of column names
  colnames(tmp) <- gsub('y_pred_WF','TEMPORARY',colnames(tmp))
  
  res <- WF_imputation(tmp, col = 'y', family = family, R_PI = R_PI)
  
  # store the results and return the original y values
  tmp <- res$df
  tmp = merge(tmp %>% select(-y),
               df %>% select(date, facility, y),
               by = c('date','facility'))

  
  # rename columns of the results
  colnames(tmp) <- gsub('y_pred_WF', 'y_pred_baseline_WF', colnames(tmp)) 
  colnames(tmp) <- gsub('TEMPORARY', 'y_pred_WF',colnames(tmp))
  
  return(tmp)
}

# Bayes imputatioFfacility_Fixedsdfsdfn method
bayes_WF_imputation <- function(df, df_OG = NULL, col, group = 'facility', family = 'NB', period = 12, iterations = 500, harmonic_priors = F){
  # prep the data with the harmonic functions
  df <- add_periodic_cov(df, period = period) %>% as.data.frame()
  
  # make the target an integer for the Bayesian approach to work.
  df[,col] = as.integer(df[,col])
  
  # pulling unique groups
  uni_group = df %>% pull(get(group)) %>% unique()
  
  # setting up the formula
  formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  
  if(harmonic_priors){
    
    if(!is.null(df_OG)){
      # setting the prior to use the un-imputed matrix models
      df_prior = add_periodic_cov(df_OG, period = period) %>% as.data.frame()
      pred_col_name = paste0(col, '_pred_bayes_WF_ideal')
    }else{
      # setting the prior to use the imputed matrix models
      df_prior = df
      pred_col_name = paste0(col, '_pred_bayes_WF_miss')
    }
    
    # get the max value
    max_val = max(df_prior[,col])
    
    # run the individual model for each group.
    tmp <- lapply(uni_group, function(xx) {
      tt <- df_prior %>% filter(get(group) == xx)
      
      # run the model
      mod_col <- MASS::glm.nb(formula_col, data = tt)
      
      # return the coefficients
      res_df = data.frame(coef(mod_col))
      colnames(res_df)[1] = xx
      
      return(data.frame(t(res_df)))
      
      # scale to get an idea of the intercept (ignoring)
      # tt[,col] = (max_val/max(tt[,col]))*tt[,col]
      
    })
    
    # get a matrix of the betas, including the intercepts.
    betas = do.call('rbind',tmp)
  }
  
  if(family == 'quasipoisson'){
    stop('havent coded Bayes quasipoisson')
    
    
  }else if(family %in% c('NB','negative binomial','neg_bin')){
    
    
    # train the model
    if(harmonic_priors){
      tmp <- lapply(uni_group, function(xx) {
        # get the data frame
        tt <- df %>% filter(get(group) == xx)
        
        # get the distribution of the other betas
        betas2 = betas[-which(rownames(betas) == xx),]
        beta_prior = normal(location = colMeans(betas2[,-1], na.rm = T), scale = apply(betas2[,-1], 2, function(xx) {sd(xx, na.rm = T)}))
        intercept_prior = normal(location = mean(betas2[,1]), scale = sd(betas2[,1]))
        
        # run the model with priors
        tt[, pred_col_name] <- tryCatch({
          # fit the model
          mod_bayes <- stan_glm(formula_col, data = tt, family = neg_binomial_2, refresh = 0, iter = iterations, prior = beta_prior, prior_intercept = intercept_prior)
          
          # need to fill in NAs, even if they aren't used
          tmp2 = tt; tmp2[is.na(tmp2)] = 1e6
          
          # predict the median output from a posterior draw
          tt[, pred_col_name] = apply(posterior_predict(mod_bayes, tmp2, type = 'response'), 2, median)
        }, error = function(e){
          print(sprinft('error for facility %s. Returning NAs', xx))
          rep(NA, nrow(tt))
        })
        
        
        # if(xx == "Family Diagnostic Clinic"){
        #   browser()
        # }
        
        return(tt)
      })
      
      # prior_summary(mod_bayes)
      # so the prior is adjusted based on scale of values. Ok ok.
      
    }else{
      tmp <- lapply(uni_group, function(xx) {
        # get the data frame
        tt <- df %>% filter(get(group) == xx)
        
        # run the models
        mod_bayes <- stan_glm(formula_col, data = tt, family = neg_binomial_2, refresh = 0, iter = iterations)
        
        # need to fill in NAs, even if they aren't used
        tmp = tt; tmp[is.na(tmp)] = 1e6
        
        # predict the median output from a posterior draw
        tt[, paste0(col, '_pred_bayes_WF')] = apply(posterior_predict(mod_bayes, tmp, type = 'response'), 2, median)
        
        return(tt)
      })
      
    }
    #browser()
    df <- data.table::rbindlist(tmp)
  }
  return(df)
}

cutoff_imputation <- function(df, df_spread = NULL, group = 'facility', method = 'number', N = NULL, cutoff = NULL){
  # if data is not wide yet, spread it!
  if(is.null(df_spread)){
    tmp = df %>% dplyr::select(date, UQ(group), indicator_denom)
    df_spread = tidyr::spread(tmp, UQ(group), indicator_denom)
  }
  
  # impute!
  df_imp = cutoff(df_spread, N = N, cutoff = cutoff, method = method)
  
  # convert data to long form for merging
  df_long = cbind(df_spread[,'date',drop = F], df_imp) %>% pivot_longer(!date, names_to = group, values_to = 'indicator_denom_pred_cutoff')
  
  # put back in original form
  df = merge(df, df_long, by = c('date',group))
  
  # return it!
  return(df)
}

CARBayes_fitting <- function(df, col, AR = 1, return_type = 'all', model = c('fixed','facility_intercept','facility_fixed'), burnin = 1000, n.sample = 2000, prediction_sample = T, thin = 10, prior = 'none', prior_var_scale = 1, prior_mean = NULL, prior_var = NULL, MALA = T, MCMC_sampler = 'stan'){

  # check if this method has already been run
  if(any(grepl('y_CARBayes_ST', colnames(df)))){
    print('previous CAR Bayes predictions found. Removing them')
    df[,grep('y_CARBayes_ST', colnames(df))] <- NULL
  }
  
  # checking that we don't have any single districts
  tt = df %>% filter(date == unique(date)[1]) %>% pull(district) %>% table
  if(any(tt == 1)){
    warning('Randomly choosing a district for a facility without that district (this should be fixed later)')
    # replace the single districts with the biggest ones
    for(nn in names(which(tt == 1))){
      df$district = gsub(nn, names(which.max(tt)), df$district)
    }
  }
  
  #districts = df %>% group_by(district) %>% summarize(n = length(unique(facility)))
  
  # create the adjacency matrix
  D2 = df %>% dplyr::select(district, facility) %>% distinct()
  W = full_join(D2, D2, by = 'district') %>%
    filter(facility.x != facility.y) %>%
    dplyr::select(-district) %>%
    igraph::graph_from_data_frame() %>%
    igraph::as_adjacency_matrix() %>%
    as.matrix()
  
  # model formula
  if(model == 'facility_intercept'){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + facility", col))
  }else if(model == 'facility_fixed'){
    formula_col = as.formula(sprintf("%s ~ facility + facility*year + facility*cos1 + facility*sin1 + facility*cos2 + facility*sin2 + facility*cos3 + facility*sin3", col))
  }else if(model == 'fixed'){
    formula_col = as.formula(sprintf("%s ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3", col))
  }else{
    print('Using user-defined model')
    formula_col = as.formula(model)
  }
  
  if(prior == 'WF'){
    # Get the WF estimates
    lm_fit <- glm(formula_col, family = 'poisson', data = df)
    coef_mat <- summary(lm_fit)$coefficients
    prior_mean_beta <- coef_mat[,1]
    prior_var_beta <- coef_mat[,2]^2*prior_var_scale
  }else if(prior == 'constant'){
    prior_mean_beta = prior_mean
    prior_var_beta = prior_var
  }else if(prior == 'none'){
    prior_mean_beta = NULL
    prior_var_beta = NULL
  }else{
    stop('please input a proper prior value (WF, constant, none)')
  }
  
  # run CAR Bayes
  if(MCMC_sampler == 'CARBayesST'){
    chain1 <- ST.CARar(formula = formula_col, 
                       family = "poisson",
                       data = df, W = W, 
                       prior.mean.beta = prior_mean_beta,
                       prior.var.beta = prior_var_beta,
                       burnin = burnin, 
                       n.sample = n.sample,
                       thin = thin, 
                       AR = AR, 
                       verbose = F,
                       MALA = MALA)
    
    # check that the prior names matched the CAR fitted names
    if(prior == 'WF'){
      if(!(identical(rownames(coef_mat), rownames(chain1$summary.results)[1:160]))){
        stop('error in WF prior and CAR formula name mismatch')
      }
    }
    
    beta_df = as.data.frame(chain1$samples$beta)
    colnames(beta_df) <- setdiff(gsub('\\(|\\)', '', row.names(chain1$summary.results)), c('tau2', 'rho.S','rho.T'))
    
    model_chain <- list(
      fitted_mean = chain1$fitted.values,
      fitted = chain1$samples$fitted,
      beta = beta_df,
      phi = chain1$samples$phi,
      rho = chain1$samples$rho[,'rho.S'],
      alpha = chain1$samples$rho[,'rho.T'],
      tau2 = chain1$samples$tau2,
      CARBayesST_summary = chain1$summary.results
    )
    
  }else if(MCMC_sampler == 'stan'){
    stan_data <- prep_stan_data_rushworth_sparse(df, formula_col)
    stan_fit <- stan(file = "R/regression_rushworth_sparse.stan",
               data = stan_data, 
               iter = n.sample, 
               warmup = burnin,
               chains = 1, 
               init = '0',
               cores = 1)
    
    # extract out the important features from the model
    stan_out <- extract(stan_fit)
    
    # pull out the beta params
    beta_df <- as.data.frame(stan_out$beta)
    model_col_names <- gsub('\\(|\\)', '', colnames(stan_data$X))
    colnames(beta_df) <- model_col_names
    
    # pull out the summary values
    stan_summary = summary(stan_fit, pars = c('tau2','rho','alpha','beta'))$summary
    rownames(stan_summary)[grep('beta', rownames(stan_summary))] <- model_col_names
    
    model_chain = list(
      fitted_mean = apply(stan_out$y_exp, 2, mean),
      fitted = stan_out$y_exp,
      beta = beta_df,
      phi = stan_out$phi,
      rho = stan_out$rho,
      alpha = stan_out$alpha,
      tau2 = stan_out$tau2,
      CARstan_summary = stan_summary
    )
  }else{
    stop('input a proper MCMC sampler')
  }
  
  df[,paste0(col, '_CARBayes_ST')] = model_chain$fitted_mean
  
  # Poisson sample the fitted values for the posterior predictive distribution
  if(prediction_sample){
    # pull the fitted values and randomly select prediction values based on the Poisson distribution
    tt = model_chain$fitted
    tt = apply(tt, 2, function(xx) rpois(n = length(xx), lambda = xx))
    
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(tt, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
    
  # ignore the sampling of the Poisson and take mean fitted value quantiles
  }else{
    # get the quantiles of fitted values 
    quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    fitted_quants = t(apply(model_chain$fitted, 2, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
  }
  
  df = as.data.frame(cbind(df, fitted_quants))

  # mapping matrix from each data point to a date
  dates = sort(unique(df$date))
  mat = as.matrix(sparseMatrix(i = 1:nrow(df), j = match(df$date, dates), x = 1))
  
  # condense the fitted value samples to each date
  date_fits = model_chain$fitted %*% mat
  
  # get the quantiles at the county level
  fitted_quants_C = t(apply(date_fits, 2, function(xx){
    quantile(xx, probs = quant_probs)
  }))
  fitted_quants_C = as.data.frame(fitted_quants_C)
  colnames(fitted_quants_C) = paste0(paste0(col, '_CARBayes_ST_'), quant_probs)
  
  # aggregate the county results
  df_county = df %>% 
    dplyr::select(date, y_true, y_CARBayes_ST) %>% 
    group_by(date) %>%
    summarize(y_true = sum(y_true),
    y_CARBayes_ST = sum(y_CARBayes_ST)) %>% 
    arrange(date)
  
  if(!identical(df_county$date, dates)){
    browser()
  }

  df_county = cbind(df_county, fitted_quants_C)
  
  if(return_type == 'data.frame'){
    return(df)
  }else if(return_type == 'model'){
    return(chain1)
  }else if(return_type == 'all'){
    lst = list(facility_df = df, county_df = df_county, model_chain = model_chain)
    return(lst)
  }
}

### Fit the model on all data using the CARBayes imputation function and return results
# df: data frame with facilities, dates, and outcomes
# R_posterior: A specified number of posterior simulations to run (if you want it to be smaller than the number of returned posterior samples from CARBayes)
# train_end_date: cutoff date for the training data to fit the model
# predict_start_date: the starting time point for where predictions should be run. If null, defaults to all dates after train_end_date
# col: outcome column
# quant_probs: quantiles to be returned from prediction samples
CARBayes_wrapper <- function(df, R_posterior = NULL, train_end_date = '2019-12-01', predict_start_date = NULL, col = 'y', quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), return_chain = F, return_raw_fit = F, ...){
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  facilities <- unique(df$facility)
  
  # get the max date
  max_date <- max(df$date)
  
  if(is.null(predict_start_date)){
    dates = df %>% 
      dplyr::filter(date > train_end_date) %>%
      dplry::select(date) %>%
      unique() %>%
      .$date
    predict_start_date = min(dates)
  }else{
    dates = df %>% 
      dplyr::filter(date >= predict_start_date) %>%
      dplyr::select(date) %>%
      unique() %>%
      .$date
  }
  
  if(predict_start_date > min(dates)){
    warning(sprintf('Setting phi to 0 at %s in the posterior prediction simulation', predict_start_date))
  }
  
  # make a train set to fit on
  train <- df %>%
    filter(date <= train_end_date) %>%
    arrange(facility, date)
  
  # fit the model!
  res <- CARBayes_fitting(train, col, ...)
  
  # return the fit if not doing post-processing
  if(return_raw_fit){
    return(res)
  }
  
  #### the rest is for future model prediction
  betas <- res$model_chain$beta
  
  phi <- res$model_chain$phi
  rho <- res$model_chain$rho
  alpha <- res$model_chain$alpha
  tau2 <- res$model_chain$tau2
  
  # set the R_posterior
  if(is.null(R_posterior)){
    R_posterior = nrow(betas)
  }else if(R_posterior > nrow(betas)){
    stop('cant sample more posterior predictions than the original model fit returns')
  }else{
    # havent implemented yet
    browser()
  }
  
  # group the betas into their separate facility values
  fac_beta_list <- list()
  beta_ref <- betas[,c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")]
  for(f in facilities){
    # if this is the reference facility (likely A1 in my simulations)
    if(sum(grepl(f, names(betas))) == 0){
      beta_f <- beta_ref
    }else{
      cols <- paste0('facility', f, c("",":year", ":cos1", ":sin1", ":cos2", ":sin2", ":cos3", ":sin3"))
      beta_f <- betas[,cols] + beta_ref 
      colnames(beta_f) <- c("Intercept", "year", "cos1", "sin1", "cos2", "sin2", "cos3", "sin3")
    }
    fac_beta_list[[f]] <- beta_f
  }
  
  # pull the fixed effects for the simulation
  fixed_effects <- lapply(facilities, function(f){
    tmp = df %>% 
      filter(facility == f,
             date >= predict_start_date) %>%
      arrange(date)
    
    # keep the 1 for intercepts
    X = tmp %>% 
      mutate(Intercept = 1) %>%
      dplyr::select(Intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
    
    rownames(X) = tmp$date
    
    # check that the ordering is correct
    if(!identical(names(X), names(fac_beta_list[[f]]))){
      browser()
    }
    
    betas <- as.matrix(fac_beta_list[[f]])
    
    # (# posterior fits x # params) x (# params x # of data points)
    mean_sims <- betas%*%t(X)
    return(mean_sims)
  })
  names(fixed_effects) <- facilities
  
  # make the "W2" matrix only once for repeated use
  W2 <- make_district_W2_matrix(df)
  
  # make the spatio-temporal precision matrices
  covar_mats <- lapply(1:R_posterior, function(ii){
    Q = make_precision_mat(df, rho = rho[ii], W2)
    
    tryCatch({
      V = tau2[ii]*solve(Q)
    }, error = function(e){
      print(e)
      print('havent dealt with non-invertible precision matrices yet')
      browser()
    })
    return(V)
  })
  
  # make R sampled sets of data
  phi_lst = lapply(1:R_posterior, function(i){
    ### get the spatio-temporal random effects
    # initialize phi
    phi = matrix(0, nrow = length(dates), ncol = length(facilities))
    colnames(phi) = facilities
    
    # first time step (at 2016-01-01, not 2020-01-01)
    phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, ncol(phi)), Sigma = covar_mats[[i]])
    
    # cycle through other time steps, using auto-correlated priors
    for(t in 2:length(dates)){
      phi[t,] = MASS::mvrnorm(n = 1, mu = alpha[i]*phi[t-1,], Sigma = covar_mats[[i]])
    }
    
    # convert to matching format
    phi = as.data.frame(phi)
    phi$date = dates
    phi_df = tidyr::gather(phi, facility, phi, setdiff(colnames(phi), c('facility','date')))
    
    return(phi_df)
  })
  
  # rearrange phi_lst to match the format of fixed effects
  phi_by_fac <- lapply(facilities, function(f){
    sapply(phi_lst, function(xx){
      xx %>% 
        filter(facility == f) %>%
        arrange(date) %>%
        pull(phi)
    })
  })
  names(phi_by_fac) <- facilities
  
  # get mean predictions at each point
  mean_pred <- lapply(facilities, function(f){
    # get mean fits
    res <- t(fixed_effects[[f]]) + phi_by_fac[[f]]
    
    # make sure dimensions line up
    if(nrow(res) != length(dates)){
      browser()
    }
    return(res)
  })
  names(mean_pred) <- facilities
  
  # predict new points using the fixed effects, phi values and poisson distribution
  predicted_vals <- lapply(facilities, function(f){
    # get mean fits
    res <- mean_pred[[f]]
    
    # get poisson sampled values from mean
    tt <- apply(res, 2, function(xx){
      rpois(length(dates), exp(xx))
    })
  })
  names(predicted_vals) <- facilities
  
  # get the quantiles of the predictions across the simulations
  fitted_quants = do.call('rbind', lapply(facilities, function(f){
    res <- as.data.frame(t(apply(predicted_vals[[f]], 1, function(xx){
      quantile(xx, probs = quant_probs)
    })))
    colnames(res) <- paste0(paste0(col, '_pred_CAR_'), quant_probs)
    res$facility = f
    res$date = dates
    return(res)
  }))
  
  fitted_quants$y_pred_CCA_CAR <- fitted_quants$y_pred_CCA_CAR_0.5
  
  # merge the results together
  df <- merge(df, fitted_quants, by = c('date', 'facility'), all = T)
  
  # initialize district results
  district_df = NULL
  
  # district-level analysis
  for(d in names(dist_fac)){
    tt = data.frame(district = d,
                    date = dates)
    facs = dist_fac[[d]]
    
    # get the sums by district: returns n x R data frame
    sum_district = Reduce('+', lapply(facs, function(f){
      predicted_vals[[f]]
    }))
    
    # get the quantiles and store them
    fitted_quants = t(apply(sum_district, 1, function(xx){
      quantile(xx, probs = quant_probs)
    }))
    fitted_quants = as.data.frame(fitted_quants)
    colnames(fitted_quants) = paste0(paste0(col, '_pred_CAR_'), quant_probs)
    
    # merge the quantiles back into data frame
    tt = cbind(tt, fitted_quants)
    
    district_df = rbind(district_df, tt)
  }
  
  res_lst <- list(df = df, district_df = district_df)
  
  if(return_chain){
    res_lst[['model_chain']] <- res$model_chain
  }
  
  if(list(...)$MCMC_sampler == 'stan'){
    res_lst[['CARstan_summary']] <- res$model_chain$CARstan_summary
  }
  
  if(list(...)$MCMC_sampler == 'CARBayesST'){
    res_lst[['CARBayesST_summary']] <- res$model_chain$CARBayesST_summary
  }
  
  return(res_lst)
}

impute_wrapper <- function(df, col = 'indicator_denom', p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5), N = 2, cutoff_cor = NULL, harmonic_family = 'NB', bayes_harmonic_family = 'NB', cutoff_method = 'number', group = 'facility', bayes_iterations = 500, bayes_beta_prior = T){
  
  if(col != 'indicator_denom'){
    warning('cutoff not coded to do imputation on other columns')
  }
  
  # create the res data frame
  res = NULL
  group_res = NULL
  
  for(i in (1:length(p_vec))){
    # simulate the missing data
    p = p_vec[i]
    print(p)
    df_miss <- simulate_imputation(df, col, p = p, group = group)
    
    # run two imputations
    df_imp <- tryCatch({
      WF_imputation(df_miss, col = col, family = harmonic_family, group = group)
    }, error = function(e){
      print(sprintf('error for periodic with p = %s', p))
      browser()
    })
    
    if(!is.na(bayes_harmonic_family)){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group)
      }, error = function(e){
        print(sprintf('error for bayes periodic with p = %s', p))
        browser()
      })
    }
    
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
      }, error = function(e){
        print(sprintf('error for bayes periodic beta prior with p = %s', p))
        browser()
      })
    }
    
    if(bayes_beta_prior){
      df_imp <- tryCatch({
        bayes_WF_imputation(df_imp, df_OG = df, col = col, family = harmonic_family, iterations = bayes_iterations, group = group, harmonic_priors = T)
      }, error = function(e){
        print(sprintf('error for bayes periodic beta prior with p = %s', p))
        browser()
      })
    }
    
    if(!is.na(cutoff_method)){
      df_imp <- tryCatch({
        cutoff_imputation(df_imp, cutoff = cutoff_cor, N = N, method = cutoff_method, group = group)
      }, error = function(e){
        print(sprintf('error for cutoff with p = %s', p))
        tmp <- df_imp
        tmp$indicator_denom_pred_cutoff = NA
        return(tmp)
      })
    }
    
    # get the "true" column names
    true_col = paste0(col,'_true')
    
    # get the eval period
    eval <- df_imp %>% filter(!is.na(get(true_col)))
    
    ### get the results over all districts
    for(imputation_method in grep(paste0(col,'_pred'), colnames(eval), value = T)){
      
      # create the res data frame and the comparison vectors
      y_true = eval[,get(true_col)]
      y_pred = eval %>% pull(imputation_method)
      
      #if(imputation_method == 'indicator_denom_pred_bayes_harmonic'){browser()}
      
      # calculate metrics
      res_row = data.frame(p = p, imputation_method = imputation_method, 
                           RMSE = sqrt(mean((y_true - y_pred)^2)), 
                           MAE = mean(abs(y_true - y_pred)), 
                           MAPE = mean(abs(y_true - y_pred)/y_true), 
                           N_imp = nrow(eval))
      
      # update full results
      res = rbind(res, res_row)
    }
    
    ### get the group level results
    # get the unique groups
    uni_group = eval %>% pull(get(group)) %>% unique()
    
    for(imputation_method in grep(paste0(col,'_pred'), colnames(eval), value = T)){
      
      # cycle through each group/facility
      for(g in uni_group){
        eval2 = eval %>% filter(get(group) == g)
        y_true = eval2 %>% pull(UQ(paste0(col,'_true')))
        y_pred = eval2 %>% pull(imputation_method)
        
        # calculate metrics
        res_row = data.frame(p = p, imputation_method = imputation_method, group = g,
                             RMSE = sqrt(mean((y_true - y_pred)^2)), 
                             MAE = mean(abs(y_true - y_pred)), 
                             MAPE = mean(abs(y_true - y_pred)/y_true), 
                             N_imp = nrow(eval2))
        
        # update full results
        group_res = rbind(group_res, res_row)
      }
      
    }
  }
  
  # return the results
  res_lst = list(full_results = res, group_results = group_res)
  return(res_lst)
}

##### freqGLMepi Imputation Functions #####
model.mean <- function(D, params){
  mu = exp(params[1])*D$y.AR1 + # auto-regressive
    exp(params[2])*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

model.mean.exp <- function(D, params){
  mu = params[1]*D$y.AR1 + # auto-regressive
    params[2]*D$y.neighbors + # neighbors
    exp(params[3] + params[4]*D$year + params[5]*D$cos1 + params[6]*D$sin1 + params[7]*D$cos2 + params[8]*D$sin2 + params[9]*D$cos3 + params[10]*D$sin3) # yearly + seasonal component
  
  return(mu)
}

# likelihood function
ll.wrapper = function(params, D, target_col){
  mu = model.mean(D, params)
  logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
  return(-logL)
}

ll.wrapper.exp = function(params, D, target_col){
  mu = model.mean.exp(D, params)
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
  
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
  #res_lst = list(variance = V, pos_def = pos_def)
  
  return(V)
}

freqGLMepi_variance2 <- function(param_vec, D){
  
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }
  
  ll.wrapper.exp = function(params, target_col = 'y_imp'){
    mu = model.mean.exp(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(-logL)
  }
  
  # get the observed information matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  # if(det(obs_I) < 0){
  #   return(list(variance = NULL, pos_def = F))
  # }else{
  #   pos_def = T
  # }
  
  # get the variance of the parameters
  #V = tryCatch({solve(obs_I)}, error = function(e) {print(det(obs_I)); browser()})
  V = solve(obs_I)
  
  # update names
  colnames(V) = names(param_vec)
  rownames(V) = names(param_vec)
  
  #res_lst = list(variance = V, pos_def = pos_def)
  
  return(V)
}

freqGLMepi_variance_testing <- function(param_vec, D){
  
  warning('weird scoping in variance function. Try to make it like the freqGLM_epi fitting procedure')
  # make likelihood function
  ll.wrapper = function(params, target_col = 'y_imp'){
    mu = model.mean(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(logL)
  }
  
  ll.wrapper.exp = function(params, target_col = 'y_imp'){
    mu = model.mean.exp(D, params)
    logL = sum(-mu + D[,target_col]*log(mu), na.rm = T)
    return(-logL)
  }
  
  
  # get the observed information matrix
  obs_I = -numDeriv::hessian(ll.wrapper, param_vec)
  
  obs_2 = -pracma::hessian(ll.wrapper, param_vec)
  
  obs_grad = numDeriv::grad(ll.wrapper, param_vec)
  
  # for numerical stability of the variance components
  param_vec_exp = param_vec
  param_vec_exp[1:2] = exp(param_vec_exp[1:2])

  obs_test = -numDeriv::hessian(ll.wrapper.exp, param_vec_exp)
  obs_test_2 =  -pracma::hessian(ll.wrapper.exp, param_vec_exp)
  
  # need to do the delta method anyway
  
  ## Analytically, lezz do it.
  
  X = cbind(1, D[,c('year','cos1','sin1','cos2','sin2','cos3','sin3')])
 
  # ugly way - these's def some matrix algebra to speed this up but ignoring this for now
  logL_hess <- matrix(0, nrow = length(param_vec), ncol = length(param_vec))
  logL_grad <- matrix(0, nrow = length(param_vec), ncol = 1)
  mu = model.mean(D, param_vec)
  for(i in 2:nrow(D)){
    y_i = D$y_imp[i]
    x_i = as.matrix(t(X[i,]))
    
    mu_grad_i = as.matrix(unlist(c(exp(param_vec[1])*D$y.AR1[i],
                  exp(param_vec[2])*D$y.neighbors[i],
                  x_i*exp(t(x_i)%*%param_vec[3:10])[1,1])))
    
    mu_grad2_i = Matrix::bdiag(list(exp(param_vec[1])*D$y.AR1[i],
                           exp(param_vec[2])*D$y.neighbors[i],
                           x_i%*%t(x_i)*exp(t(x_i)%*%param_vec[3:10])[1,1]))

    logL_grad <- logL_grad + (y_i/mu[i] - 1)*mu_grad_i
    logL_hess <- logL_hess + (y_i/mu[i] - 1)*mu_grad2_i - y_i/(mu[i]^2)*mu_grad_i %*% t(mu_grad_i)
  }
  
  
  if(norm(obs_I + logL_hess)/norm(logL_hess) > .01){
    browser()
  }

  if(det(obs_I) < 0){
    if(det(logL_hess) < 0){
      browser()
    }else{
      obs_I = -logL_hess
    }
    
  }
  
  # get the variance of the parameters
  #V = tryCatch({solve(obs_I}, error = function(e) {browser()})})
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
    init = rep(0,10)
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  # NM_control = list(maxit = 10000, reltol = 1e-12)
  NM_control = list(maxit = 10000)
  # BFGS_control = list(maxit = 10000, factr = 1e-11))
  BFGS_control = list(maxit = 10000)
  
  # fit using Nelder-Mead and L-BFGS-B and pick the better one
  tryCatch({
    params = optim(par = init, fn = ll.wrapper, D = df, target_col = target_col, control = NM_control)
  }, error = function(e){
    print(e)
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
      
      init = init_OG + rnorm(10,0,5*i/num_inits)
      
      # nelder-mead
      tryCatch({
        params2 = optim(init, ll.wrapper, D = df, target_col = target_col, control = NM_control)
      }, error = function(e){
        browser()
      })
      
      if(params2$value < params$value & params2$convergence == 0){
        if(verbose){print('using another initialization')}
        params = params2
      }
      
      # L-BFGS
      if(BFGS){
        tryCatch({
          params2 = optim(init, ll.wrapper,  D = df, target_col = target_col, method = 'L-BFGS-B',  control = BFGS_control)
          if(params2$value < params$value & params2$convergence == 0){
            if(verbose){print('using another initialization')}
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
    if(verbose){print('didnt converge for one iteration')}
    #browser()
  }
  
  if(verbose){
    print(sprintf('freqGLMepi fit in %s seconds', Sys.time() - t0))
  }
  
  names(params$par) = parnames
  return(params)
}

### the NNLS freqGLMepi model (so lambda and phi are constrained rather than exponentiated)
fit_freqGLMepi_nnls <- function(df, num_inits = 10, verbose = T, target_col = 'y_imp', init = NULL, BFGS = NULL){
  t0 = Sys.time()
  
  # set up initialization
  if(is.null(init)){
    init = c(0.1, 0.1, rep(0,8))
  }
  
  parnames = c('By.AR1', 'By.neighbors', 'Bintercept', 'Byear', 'Bcos1', 'Bsin1', 'Bcos2', 'Bsin2', 'Bcos3', 'Bsin3')
  names(init) = parnames
  
  init_OG = init
  
  # params = nlminb(start = init, objective = ll.wrapper, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
  
  tryCatch({
    params = nlminb(start = init, objective = ll.wrapper.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
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
        params2 = nlminb(start = init, objective = ll.wrapper.exp, D = df, target_col = target_col, lower = c(0, 0, rep(-10, 10)), upper = c(1, 1, rep(10, 10)))
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

### Given data with missing values, fits the freqGLMepi model on all data points
#   df: data
#   max_iter: number of iterations of imputation
#   tol: tolerance of the converage of imputed values
#   individual_facility_models: fit each facility separately

#   R_PI: number of bootstrap iterations if doing so
#   quant_probs: the quantiles of the bootstrap to store in the data frame
#   verbose: printing updates
freqGLMepi_CCA = function(df, train_end_date = '2019-12-01', max_iter = 1, tol = 1e-4, individual_facility_models = T, family = 'poisson', R_PI = 100, quant_probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), verbose = F, optim_init = NULL, scale_by_num_neighbors = T, blocksize = 10, nnls = T){
  # check that we have the right columns
  if(!('y' %in% colnames(df) & 'y_true' %in% colnames(df))){
    stop('make sure the data has y (with NAs) and y_true')
  }
  
  # get districts and facilities
  dist_fac <- get_district_facilities(df)
  
  # convert to proper format
  train_end_date <- as.Date(train_end_date)
  
  # split the hold out set and train set
  df_test <- df %>%
    filter(date > train_end_date)
  df_test$y_imp <- NA
  df_test$y_pred_freqGLMepi <- NA
  df <- df %>%
    filter(date <= train_end_date)
  
  # get the adjacency matrix
  W <- make_district_adjacency(df, scale_by_num_neighbors)
  
  # set fitting function to be regular (exponentiated spatial and temporal parameters) or non-negative least squares 
  if(nnls){
    fit_function <- fit_freqGLMepi_nnls
    model_function <- model.mean.exp
  }else{
    fit_function <- fit_freqGLMepi
    model_function <- model.mean
  }
  
  y_pred_list = list()
  
  ### Do initial filling of y
  # setting up the formula
  formula_col = as.formula("y ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  
  # the unique facility groups
  uni_group = unique(df$facility)
  
  # run the individual model for each group.
  if(family == 'poisson'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(facility == xx)
      
      # run the model
      mod_col <- glm(formula_col, family = 'poisson', data = tt)
      
      # update predictions
      tt$y_pred_freqGLMepi = predict(mod_col, tt, type = 'response')
      
      return(tt)
    })
  }
  else if(family == 'negative_binomial'){
    tmp <- lapply(uni_group, function(xx) {
      tt <- df %>% filter(facility == xx)
      
      # run the model
      mod_col <- MASS::glm.nb(formula_col, data = tt)
      
      # update predictions
      tt$y_pred_freqGLMepi = predict(mod_col, tt, type = 'response')
      
      return(tt)
    })
  }
  # combine into one data frame
  df = do.call('rbind',tmp)
  
  # filling in missing values by randomly sampling mean prediction from Poisson
  df$y_imp = df$y
  na.ind = which(is.na(df$y))
  df$y_imp[na.ind] <- rpois(n = length(na.ind), df$y_pred_freqGLMepi[na.ind])
  
  # add the neighbors and auto-regressive
  df = add_autoregressive(df, 'y_imp') %>%
    add_neighbors(., 'y_imp', scale_by_num_neighbors = scale_by_num_neighbors, W = W)

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
        params = fit_function(tt, verbose = verbose, init = optim_init[[xx]])
        
        # update y_pred
        tt$y_pred_freqGLMepi = model_function(tt, params$par)
        
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
      parmat = fit_function(df)$par
      
      # update y_pred 
      df$y_pred_freqGLMepi = model_function(df, parmat)
      
    }
    
    # store the predictions for this iteration
    y_pred_list[[iter]] = df$y_pred_freqGLMepi
    
    if(length(na.ind) == 0){
      if(verbose){print('only running one iteration because there is no missingness')}
      break
    }
    
    # update y_imp
    na.ind.2 = intersect(na.ind, which(!is.na(df$y_pred_freqGLMepi)))
    df$y_imp[na.ind.2] <- rpois(n = length(na.ind.2), df$y_pred_freqGLMepi[na.ind.2])
    
    # compare y_imps
    prop_diffs = c(prop_diffs, mean(abs(y_pred_list[[iter]][na.ind] - y_pred_list[[iter-1]][na.ind])/y_pred_list[[iter-1]][na.ind], na.rm = T))

    # update the neighbors and auto-regressive
    df = add_autoregressive(df, 'y_imp') %>%
      add_neighbors(., 'y_imp', W = W)
    
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
  
  ### run the stationary bootstrap
  param_boots <- lapply(1:length(uni_group), function(i) {
    # subset data
    tt <- df %>% filter(facility == uni_group[i])
    tt_test <- df_test %>% filter(facility == uni_group[i])
    
    predict_function <- function(xx){
      # fit model on bootstrapped data set
      # start the initialization at the values from the main model
      params = fit_function(xx, BFGS = F, num_inits = 1, verbose = verbose, init = parmat[,1])
      
      # predict model on original training data set
      x = suppressWarnings(rpois(n = nrow(tt), model_function(tt, params$par)))
      
      return(params$par)
    }
    
    # run the stationary bootstrap and return the parameters from each fit
    sim_boot <- boot::tsboot(tt, statistic = predict_function, R = R_PI, sim = 'geom', l = blocksize)$t
    
    # match the parameter names
    colnames(sim_boot) <- names(fit_function(tt, BFGS = F, num_inits = 1, verbose = verbose, init = parmat[,1])$par)
    
    return(sim_boot)
  })
  names(param_boots) <- as.character(uni_group)
  
  # get the parameter results from the bootstrap
  param_results <- do.call('rbind', lapply(as.character(uni_group), function(fac){
    tt <- param_boots[[as.character(fac)]]
    tmp <- data.frame(t(apply(tt, 2, function(col) {
        #mean(col)
      vec = c(mean(col), quantile(col, probs = c(.025,0.25,0.5,0.75,0.975), na.rm = T))
      names(vec) <- c('mean', paste0('Q_', names(vec)[-1]))
      vec
    })))
    tmp$param = rownames(tmp)
    tmp$facility = fac
    tmp
  }))
  param_results <- param_results[,c(8,7,1:6)]
  rownames(param_results) <- NULL
  
  # create long form of estimated parameters
  tmp = parmat %>%
    as.data.frame(.) %>%
    mutate(param = rownames(.)) 
  par_long = NULL
  for(i in 1:ncol(tmp)){
    par_long = rbind(par_long,
                     data.frame(facility = colnames(tmp)[i],
                                param = tmp$param, 
                                full_estimate = tmp[,i]))
  }

  # merge them!
  param_results = merge(param_results, par_long, by = c('facility','param'))
  
  # rename columns appropriately
  param_results$param = gsub('y.neighbors','rho',
                             gsub('y.AR1', 'alpha',
                                  gsub('B','',param_results$param)))
  
  # combine the data sets and split by facility
  df_combined <- rbind(df[,colnames(df_test)], df_test) %>%
    add_autoregressive(., 'y_imp') %>%
    add_neighbors(., 'y_imp', scale_by_num_neighbors = scale_by_num_neighbors, W = W) %>%
    arrange(date)
  df_combined$y_pred <- NA
  
  # cycle through all the bootstrap iterations
  pred_boots <- lapply(1:nrow(param_boots[[1]]), function(i) {
    # resetting the data frame
    df_tmp <- df_combined 
    
    # do all the baseline predictions (+ 1)
    for(f in uni_group){
      tt <- df_tmp %>% filter(facility == f, date <= (train_end_date + 31))
      
      y = suppressWarnings(rpois(n = nrow(tt), model_function(tt, param_boots[[f]][1,]))) 
      
      # put the results back into the data frame
      df_tmp$y_pred[df_tmp$facility == f][1:length(y)] <- y
    }
    
    # update the neighbors and auto-regressive terms
    df_tmp <- add_autoregressive(df_tmp, 'y_pred') %>%
      add_neighbors(., 'y_pred', scale_by_num_neighbors = scale_by_num_neighbors, W = W)
    
    # get remaining time points
    time_points <- unique(df_tmp$date[is.na(df_tmp$y_pred) & df_tmp$date > train_end_date]) 
    
    # cycle through remaining time points
    for(t in time_points){
      for(f in uni_group){
        tt <- df_tmp %>% filter(facility == f, date == t)
        
        y = suppressWarnings(rpois(n = nrow(tt), model_function(tt, param_boots[[f]][i,]))) 
        
        # put the results back into the data frame
        df_tmp$y_pred[df_tmp$facility == f & df_tmp$date == t] <- y
      }
      
      # update the neighbors and auto-regressive terms
      df_tmp <- add_autoregressive(df_tmp, 'y_pred') %>%
        add_neighbors(., 'y_pred', scale_by_num_neighbors = scale_by_num_neighbors, W = W)
    }
    
    return(df_tmp$y_pred)
  })
  
  # combine into a data frame
  pred_boots <- do.call('cbind', pred_boots)
  
  # get the quantiles and store them
  fitted_quants = t(apply(pred_boots, 1, function(xx){
    quantile(xx, probs = quant_probs, na.rm = T)
  }))
  fitted_quants = as.data.frame(fitted_quants)
  colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)

  # combine the final results and return
  df_combined = cbind(df_combined, fitted_quants)
  
  # set the prediction to the median
  df_combined$y_pred_CCA_freqGLMepi <- df_combined$y_pred_CCA_freqGLMepi_0.5

  # remove the ambiguous "y_pred" column
  df_combined$y_pred <- NULL
  
  ### Make the district results
  district_df = data.frame(cbind(df_combined[,c('date','district')], pred_boots)) %>% 
    group_by(date, district) %>%
    summarize_all(sum)
  
  district_mat = district_df[,3:ncol(district_df)]
  # get the quantiles and store them
  fitted_quants = t(apply(district_df[,3:ncol(district_df)], 1, function(xx){
    quantile(xx, probs = quant_probs, na.rm = T)
  }))
  fitted_quants = as.data.frame(fitted_quants)
  colnames(fitted_quants) = paste0(paste0('y_pred_freqGLMepi_'), quant_probs)
  
  district_df = cbind(district_df[,c('date','district')], fitted_quants)
  
  # prep data to return
  return_lst = list(df = df_combined, district_df = district_df, param_results = param_results, convergence = convergence, y_pred_list = y_pred_list, prop_diffs = prop_diffs)

  return(return_lst)
}

#
##### Plotting + Metric Functions #####

plot_results <- function(res, subset = NULL){
  tt = res$group_results %>%
    mutate(p = as.factor(p),
           imputation_method = gsub('observed_count_pred_', '', imputation_method)) %>%
    mutate(imputation_method = factor(imputation_method, levels = c("harmonic", "bayes_harmonic", "bayes_harmonic_miss", "bayes_harmonic_ideal")))
  
  # subset data if that's what's asked for
  if(!is.null(subset)){
    tt = tt %>% filter(imputation_method %in% subset)
  }
  
  # RMSE plot
  p1 <- ggplot(tt, aes(x = p, y = RMSE, fill = imputation_method)) +
    geom_boxplot() + 
    ggtitle('RMSE of imputation predictions') + 
    ylim(c(0,1500))
  
  # MAPE plot
  p2 <- ggplot(tt, aes(x = p, y = MAPE, fill = imputation_method)) +
    geom_boxplot() + 
    ggtitle('Mean Absolute Prediction Error of imputation predictions') + 
    ylim(c(0,1))
  
  
  # final plot
  plot_final = cowplot::plot_grid(p1, p2, ncol = 1)
  return(plot_final)
}

plot_county_imputations <- function(df, methods, color_vec, title = 'Aggregated County Imputations', labels = NULL){
  # set the labels for the plot
  if(is.null(labels)){
    labels = methods
  }
  
  # rename columns of df to make it easier
  for(col in methods){
    ind = grep(col, colnames(df))
    if(length(ind) != 1){browser()}
    colnames(df)[ind] = col
  }
  
  # initialize the data frame to store final results
  df_f = NULL
  
  for(j in 1:length(methods)){
    col = methods[j]
    
    # create the imputed outcome
    df$imputed_outcome = df$indicator_count_ari_total
    df$imputed_outcome[is.na(df$imputed_outcome)] = df[is.na(df$imputed_outcome), col]
    
    # aggregate by date  
    tmp = df %>% 
      group_by(date) %>%
      summarize(y = sum(imputed_outcome)) %>% mutate(method = col)
    
    # store results for this method
    df_f = rbind(df_f, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_f$method = factor(df_f$method, levels = methods)
  
  # plot them all!
  p1 <- ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
    geom_line() +
    ggtitle(title) + 
    scale_color_manual(values = color_vec, labels = labels) + theme(legend.position = 'bottom', legend.text=element_text(size=20))
  
  return(p1)
}

plot_county_fits_from_fac <- function(df, methods, color_vec, title = 'Aggregated County Fits', labels = NULL){
  
  # set the labels for the plot
  if(is.null(labels)){
    labels = methods
  }
  
  # # rename columns of df to make it easier
  # for(col in methods){
  #   ind = grep(col, colnames(df))
  #   if(length(ind) != 1){browser()}
  #   colnames(df)[ind] = col
  # }
  
  # initialize the data frame to store final results
  df_f = NULL

  for(j in 1:length(methods)){
    col = methods[j]

    # aggregate by date
    tmp = df %>%
      group_by(date) %>%
      summarize(y = sum(get(col))) %>% mutate(method = col)

    # store results for this method
    df_f = rbind(df_f, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_f$method = factor(df_f$method, levels = methods)
  
  # plot them all!
  p1 <- ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
    geom_line() +
    ggtitle(title) + 
    scale_color_manual(values = color_vec, labels = labels) + theme(legend.position = 'bottom', legend.text=element_text(size=20))
  
  return(p1)
}

plot_imputations <- function(df, methods, color_vec, fac_list = NULL){
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # rename columns of df to make it easier
  for(col in methods){
    ind = grep(col, colnames(df))
    if(length(ind) != 1){browser()}
    colnames(df)[ind] = col
  }
  
  # initialize plotting
  plot_list = list()
  iter = 0
  
  # go through each facility
  for(f in fac_list){
    iter = iter + 1
    tmp = df %>% filter(facility == f)
    
    # store the true outcome for plotting
    df_f = tmp %>% dplyr::select(date, y = indicator_count_ari_total) %>% mutate(method = 'ari_true')
    
    #  # baseline plot for this
    # # p1 <- suppressWarnings(ggplot(tmp, aes(x = date, y = indicator_count_ari_total)) + 
    #                           geom_line() +
    #                           ggtitle(sprintf('%s (%0.1f %% M)', f, 100*mean(is.na(tmp$indicator_count_ari_total)))))
    for(j in 1:length(methods)){
      col = methods[j]
      tmp[!is.na(tmp$indicator_count_ari_total),col] = NA
      
      # filling in surrounding points to imputations to connect lines
      # I have to do multiple if statements because of R's handling of calling a value not in an array
      ind_list = c()
      for(i in 1:nrow(tmp)){
        if(i > 1){
          if(is.na(tmp[i,col]) & !is.na(tmp[i-1,col])){
            ind_list = c(ind_list, i)
          }
        }
        if(i < nrow(tmp)){
          if(is.na(tmp[i,col]) & !is.na(tmp[i+1,col])){
            ind_list = c(ind_list, i)
          }
        }
      }
      # replace surrounding values for plotting
      tmp[ind_list, col] = tmp$indicator_count_ari_total[ind_list]
      
      tmp2 = tmp %>% dplyr::select(date, y = methods[j]) %>% mutate(method = col)
      df_f = rbind(df_f, tmp2)
      
      #scale_colour_manual(name="Error Bars",values=cols)
      #p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = methods[j]), color = color_vec[j]))
      # p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = methods[j], color = color_vec[j])))
    }
    
    #browser()
    
    p1 <- suppressWarnings(ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
                             geom_line() +
                             ggtitle(sprintf('%s (%0.1f%% M)', f, 100*mean(is.na(tmp$indicator_count_ari_total)))) + 
                             scale_color_manual(values = c('black', color_vec)))
    
    # store the legend for later
    legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
    
    # remove the legend position on this plot
    p1 = p1 + theme(legend.position = 'none') 
    
    # store the plot for this facility in the list
    plot_list[[iter]] = p1
  }
  
  #browser()
  plot_grid(plot_grid(plotlist = plot_list),  legend, ncol = 1, rel_heights = c(10,1))
  
}

# HERE! FOR COMMENTING PURPOSES. HAVENT COMMENTED THE FUNCTIONS ABOVE

plot_county_fits <- function(df, methods, color_vec, imp_names = NULL, PIs = T, title = 'County-Level Predictions'){
  df = as.data.frame(df)
  
  # if no imputation names given, use the ones in the imputation vector
  if(is.null(imp_names)){
    imp_names = methods
  }
  
  # initialize data frame for this county
  df_c = NULL
  
  # pull estimates for each method
  for(j in 1:length(methods)){
    col = methods[j]
    
    # get the lower and upper bounds
    tmp = df[,c('date',paste0(col,'_0.5'),paste0(col, '_0.025'),paste0(col, '_0.975'))] 
    colnames(tmp) = c('date', 'y', 'y_lower', 'y_upper')
    tmp$method = imp_names[j]
    
    df_c = rbind(df_c, tmp)
  }
  
  # ordering the method to be consistent and for the labeling
  df_c$method = factor(df_c$method, levels = imp_names)
  
  # make the plot!
  p1 <- ggplot() +
    geom_line(data = df, aes(x = date, y = y_true), size = 1) +
    geom_line(data = df_c, aes(x = date, y = y, group = method, color = method)) +
    geom_ribbon(data = df_c, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.3) +
    scale_color_manual(values = c(color_vec)) + 
    scale_fill_manual(values = c(color_vec)) + 
    ggtitle(title) + 
    ylab('y') +
    theme_bw() + 
    theme(text = element_text(size = 15))
  
  # store the legend for later
  legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
  
  # remove the legend position on this plot
  #p1 = p1 + theme(legend.position = 'none') 
  return(p1)
}

# plot the fits of the models for a group of facilities. Also plots the prediction intervals.
# Also can plot the raw values if no imputations are provided (methods is NULL)
### Parameters
plot_facility_fits <- function(df, methods = NULL, imp_names = NULL, color_vec = NULL, PIs = T, fac_list = NULL, plot_missing_points = T, vertical_line = '2020-01-01', outbreak_points = NULL, include_legend = T, ...){
  df = as.data.frame(df)
  
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # if no imputation names given, use the ones in the imputation vector
  if(is.null(imp_names)){
    imp_names = methods
  }
  
  # initialize plotting
  plot_list = list()
  iter = 0
  
  # go through each facility
  for(f in fac_list){
    iter = iter + 1
    tmp = df %>% filter(facility == f)
    
    if(!is.null(methods)){
      # initialize data frame for this facility
      df_f = NULL
      
      for(j in 1:length(methods)){
        col = methods[j]
        
        # get the lower and upper bounds
        tmp2 = tmp[,c('date',col,paste0(col, '_0.025'),paste0(col, '_0.975'))] 
        colnames(tmp2) = c('date', 'y', 'y_lower', 'y_upper')
        tmp2$method = imp_names[j]
        
        df_f = rbind(df_f, tmp2)
      }
      
      # ordering the method to be consistent and for the labeling
      df_f$method = factor(df_f$method, levels = imp_names)
      
      # get the ymax value
      ymax <- min(1.1*max(tmp[,c(12:ncol(tmp))], na.rm = T), 
                  2*max(tmp$y_true))
      
      # make the plot!
      p1 <- ggplot() +
        geom_line(data = tmp, aes(x = date, y = y), size = 1) + 
        geom_line(data = df_f, aes(x = date, y = y, group = method, color = method), show.legend = include_legend) +
        #scale_color_manual(values = c(color_vec)) + 
        #scale_fill_manual(values = c(color_vec)) + 
        ylim(c(0,ymax)) +
        # ggtitle(sprintf('facility %s', f)) +
        ggtitle(f) + 
        ylab('y') +
        theme_bw() +
        theme(text = element_text(size = 10))
      
      if(PIs){
        p1 <- p1 + geom_ribbon(data = df_f, aes(x = date,ymin = y_lower, ymax = y_upper, fill = method, colour = method), alpha = 0.1, show.legend = F)
      }
      
      if(!is.null(color_vec)){
        p1 <- p1 + 
          scale_color_manual(values = c(color_vec)) + 
          scale_fill_manual(values = c(color_vec)) 
      }
      
      if(plot_missing_points){
        tmp2 <- tmp %>%
          filter(is.na(y)) %>%
          select(date)
        
        p1 <- p1 + 
          geom_point(data = tmp2, aes(x = date, y = 0),  color = 'red', size = 3)
      }
      
      # plot vertical line at outbreak time
      if(!is.null(vertical_line)){
        p1 <- p1 + geom_vline(xintercept = as.Date(vertical_line))
      }
      
      # plot outbreak points
      if(!is.null(outbreak_points)){
        ind = which(tmp$date == '2020-01-01')
        out_df = data.frame(date = rep(as.Date('2020-01-01'), length(outbreak_points)),
                            outbreak = tmp$y_exp[ind] + outbreak_points*sqrt(tmp$y_var[ind]),
                            k = factor(outbreak_points))
        p1 <- p1 + geom_point(data = out_df, aes(x = date, y = outbreak, shape = k), size = 3) + 
          guides(shape = guide_legend(title = "Outbreak size k:"))
      }
      
      # store the legend for later
      legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20), 
                                     legend.title = element_text(size = 20)))
      
      # remove the legend position on this plot
      p1 <- p1 + theme(legend.position = 'none') 
      
      # store the plot for this facility in the list
      plot_list[[iter]] = p1
    }else{
      # print('plotting baseline counts because no imputation methods are provided')
      p1 <- ggplot() +
        geom_line(data = tmp, aes(x = date, y=y_true), col = 'red') + 
        geom_line(data = tmp, aes(x = date, y = y), size = 1) + 
        ylim(c(0,1.5*max(tmp$y))) + 
        ggtitle(sprintf('facility %s', f)) + 
        ylab('y') +
        theme_bw() +
        theme(text = element_text(size = 10))
        
      # plot vertical line at outbreak point
      if(!is.null(vertical_line)){
        p1 <- p1 + geom_vline(xintercept = as.Date(vertical_line))
      }
      
      # plot outbreak points
      if(!is.null(outbreak_points)){
        ind = which(tmp$date == '2020-01-01')
        out_df = data.frame(date = rep(as.Date('2020-01-01'), length(outbreak_points)),
                            outbreak = tmp$y_exp[ind] + outbreak_points*sqrt(tmp$y_var[ind]),
                            k = factor(outbreak_points))
        p1 <- p1 + geom_point(data = out_df, aes(x = date, y = outbreak, shape = k), size = 3, col = 'red' )
      }
      
      # store the legend for later
      legend = get_legend(p1 + theme(legend.position = 'bottom', legend.text=element_text(size=20)))
      
      # remove the legend position on this plot
      p1 <- p1 + theme(legend.position = 'none') 
      
      plot_list[[iter]] = p1
    }
    
  }
  
  if(include_legend){
    final_plot <- cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_list, ...),legend, ncol = 1, rel_heights = c(10,1))
  }else{
    final_plot <- cowplot::plot_grid(plotlist = plot_list, ...)
  }
  return(final_plot)
}

# process the imputed list for metric calculations
clean_data_list <- function(imputed_list, dates = '2020-01-01',  min_date = NULL, rm_ARna = F, imputed_only = F){
  # filter to only be greater than the specified date
  if(!is.null(min_date)){
    # print(sprintf('only getting metrics with dates on or after %s', min_date))
    imputed_list <- lapply(imputed_list, function(xx){
      xx <- xx %>% dplyr::filter(date >= min_date)
    })
  }
  
  # filter to only be the specified date
  if(!is.null(dates)){
    # print(sprintf('only getting metrics with dates on  %s', dates))
    imputed_list <- lapply(imputed_list, function(xx){
      xx <- xx %>% dplyr::filter(date == dates)
    })
  }
  
  # removing the starting points with NA AR1 values, since these
  if(rm_ARna){
    print('removing the starting points because of NA autoregressive term')
    imputed_list = lapply(imputed_list, function(xx) xx[!is.na(xx$y.AR1),])
  }
  
  # remove all non-missing points
  if(imputed_only){
    imputed_list = lapply(imputed_list, function(xx) xx[is.na(xx$y),])
  }
  
  return(imputed_list)
}

# calculate the metrics across simulations
calculate_metrics <- function(imputed_list, methods = c("y_pred_WF", "y_CARBayes_ST"), results_by_point = F, date = '2020-01-01', min_date = NULL, rm_ARna = F, imputed_only = F,  use_point_est = F, k = NULL, district_results = F, QP_variance_adjustment){

  if(results_by_point){
    # getting the results for each point across all simulations
    avg_fxn = rowMeans
  }else{
    # getting the results for each simulation across all points
    avg_fxn = colMeans
  }
  
  # process the imputed_list
  imputed_list = clean_data_list(imputed_list, date, min_date, rm_ARna, imputed_only)
    
  ## Get the outcome values across all simulations
  {
  if(district_results){
    # updating the y_true and y missing since we don't need those
    y_true = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y']))
    y_missing = matrix(NA, nrow = nrow(y_true), ncol = ncol(y_true))
  }else{
    # get the true values everywhere and at the deleted time points
    y_true = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y_true']))
    y_missing = do.call('cbind', lapply(imputed_list, function(xx) {
      y_true = xx[,'y_true'];
      y_true[!is.na(xx[,'y'])] = NA
      y_true
    }))
  }
  
  # get the expected y values and the variance associated with them
  y_exp = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y_exp']))
  if(is.null(QP_variance_adjustment)){
    y_var = tryCatch({
      y_var = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y_var']))
    }, error = function(e){
      warning('making var(Y) = E(Y). This DOES NOT hold under the CAR DGP')
      y_var = y_exp
    })
  }else{
    y_var = do.call('cbind', lapply(imputed_list, function(xx) xx[,'y_var']))
    # check that it's equal to y_exp (as in QP not already adjusted for)
    if(sum(y_var != y_exp) > 0){
      stop('QP variance adjustment error: Has it already been adjusted?')
    }else{
      y_var = QP_variance_adjustment*y_var
    }
  }
  
  # calculate the outbreak values
  y_outbreak3 <- y_exp + 3*sqrt(y_var)
  y_outbreak5 <- y_exp + 5*sqrt(y_var)
  y_outbreak10 <- y_exp + 10*sqrt(y_var)
  
  # numeric missing matrix
  missing_mat <- apply(y_missing, 2, function(xx) 1 - as.numeric(is.na(xx)))
  missing_mat_NA <- missing_mat; missing_mat_NA[missing_mat_NA == 0] <- NA
  
  # get the number of times each data point was missing across simulations
  num_missing = apply(y_missing, 1, function(xx) sum(!is.na(xx)))
  }
  
  df = NULL
  for(method in methods){
    lower_025 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.025')]))
    upper_975 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.975')]))
    lower_25 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.25')]))
    upper_75 = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.75')]))
    
    if(use_point_est){
      point_est = do.call('cbind', lapply(imputed_list, function(xx) xx[,method]))
      #tmp$point_est <- sapply(1:nrow(point_est), function(ii) mean(point_est[ii,]))
      outcome = point_est
    }else{
      median = do.call('cbind', lapply(imputed_list, function(xx) xx[,paste0(method, '_0.5')]))
      #tmp$median <- sapply(1:nrow(median), function(ii) mean(median[ii,]))
      outcome = median
    }
    
    # prep the temporary results data frame
    if(results_by_point){
      if(district_results){
        tmp = imputed_list[[1]][,c('date','district')]
      }else{
        tmp = imputed_list[[1]][,c('date','facility','district')]
      }
      tmp$num_missing = num_missing

    }else{
      tmp = data.frame(r = 1:length(imputed_list))
    }
    
    tmp$method = method
    
    # point estimation metrics
    tmp$bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {outcome[,ii] - y_true[,ii]}))
    tmp$relative_bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_true[,ii])/y_exp[,ii]}))
    tmp$absolute_bias = avg_fxn(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])}))
    tmp$MAPE = avg_fxn(sapply(1:ncol(outcome), function(ii) {abs(outcome[,ii] - y_true[,ii])/y_true[,ii]}))
    tmp$RMSE = sqrt(avg_fxn(sapply(1:ncol(outcome), function(ii) {(outcome[,ii] - y_true[,ii])^2})))
    
    # coverage metrics and specificity
    tmp$coverage50 = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] >= lower_25[,ii] & y_true[,ii] <= upper_75[,ii])))
    tmp$coverage95 = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] >= lower_025[,ii] & y_true[,ii] <= upper_975[,ii])))
    tmp$specificity = avg_fxn(sapply(1:ncol(lower_25), function(ii) (y_true[,ii] <= upper_975[,ii])))
    
    # outbreak detection metrics
    tmp$outbreak_detection3 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak3[,ii] >= upper_975[,ii]))
    tmp$outbreak_detection5 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak5[,ii] >= upper_975[,ii]))
    tmp$outbreak_detection10 <- avg_fxn(sapply(1:ncol(y_exp), function(ii) y_outbreak10[,ii] >= upper_975[,ii]))
    
    # measure of how wide the 95% prediction intervals are
    tmp$interval_width = avg_fxn(upper_975 - lower_025) 
    tmp$prop_interval_width = avg_fxn((upper_975 - lower_025)/y_true)
    
    # update the results
    df = rbind(df, tmp)
  }
  
  #res_lst = list(df = df, num_missing = num_missing)
  return(df)
}

#### plot all methods against each other
### Inputs:
## files: Files to pull the results from. Calls these if "res" is null
## res: the results from the simulation runs
## fix_axis: a list/vector containing y limits. Can be of the form list(ylim(0,1), F, ylim(0.9,1), F)
## add_lines: a vector of where to add horizontal lines in plots (if F no lines are added)
## bar_quants: the outer quantiles of the metric results
## metrics: which metrics to plot
## metric_rename: a vector of what to ttile the metrics
## rows: number of rows of resulting plots
## title: title of overall plot
## ...: params to be passed into "combine_results_wrapper"

plot_all_methods <- function(files = NULL, res = NULL, fix_axis = F, add_lines = rep(F, 4), bar_quants = c(0.25, 0.75), metrics = c('specificity', 'outbreak_detection3', 'outbreak_detection5', 'outbreak_detection10'), metric_rename = c('specificity', 'sensitivity-3', 'sensitivity-5', 'sensitivity-10'), rows = 2, title = NULL, include_legend = T, ...){
  if(is.null(res)){
    if(!is.null(files)){
      tmp <- combine_results_wrapper(files,  ...)
      res = tmp$results
    }else{
      stop('need files if not providing results')
    }
  }
  
  if((length(metric_rename) != length(metrics)) & (length(metrics) >0) & (length(metric_rename) > 0)){
    stop('metric rename and metrics dont match in length')
  }
  
  # get color set up
  method_colors <- setNames(c('orange3', 'forestgreen', 'blue'), levels(res$method))
  
  # fixing the axis dimensions
  if(length(fix_axis) == 1){
    fix_axis = rep(fix_axis, length(metrics))
  }
  
  options(dplyr.summarise.inform = FALSE)
  
  # making the length of the alphas work for odd and even numbers
  stripes = as.factor(rep(c(0, 1), (length(unique(res$prop_missing)) + 1)/2)[1:length(unique(res$prop_missing))])
  rects = data.frame(xstart = seq(-0.05, max(res$prop_missing) - 0.05, by = 0.1),
                     xend = seq(0.05, max(res$prop_missing) + 0.05, by = 0.1),
                     stripe = stripes)
  
  # rename the metrics
  if(!is.null(metric_rename)){
    for(i in 1:length(metrics)){
      colnames(res)[which(colnames(res) == metrics[i])] <- metric_rename[i]
    }
    metrics <- metric_rename
  }
  
  plot_list = list()
  i = 0
  for(metric in metrics){
    i = i + 1
    
    # aggregate the results by the metric and quantile of plots.
    tmp <- res %>%
      group_by(method, prop_missing) %>% 
      summarize(median = median(get(metric)),
                lower = stats::quantile(get(metric), probs = bar_quants[1]),
                upper = stats::quantile(get(metric), probs = bar_quants[2])) 
    
    # remove parts of string not wanted in plot.
    tmp$method <- gsub('y_pred_|_MCAR', '', tmp$method)
    
    # refactor for ordering
    if(class(res$method) == 'factor'){
      tmp$method = factor(tmp$method, levels = levels(res$method))
    }
    
    # plot!
    p1 <- ggplot() + 
      # plot the background shading
      geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = stripe), alpha = 0.2,show.legend = F)  + scale_fill_manual(values = c('white', 'grey50')) + 
      # plot the median points
      geom_point(data = tmp, aes(x = prop_missing, y = median, color = method, shape = method), position = position_dodge(width = 0.1)) +
      # plot the error bars
      geom_errorbar(data = tmp, aes(x = prop_missing, ymin = lower, ymax = upper, color = method), position = position_dodge(width = 0.1)) + 
      scale_color_manual(values = method_colors) + 
      labs(x = 'proportion missing') + 
      guides(alpha = 'none') +
      theme_bw() + 
      ylab(metric) + guides(color = guide_legend(title = 'model fit', override.aes = list(size = 4)), shape = guide_legend(title = 'model fit'))
    
    # fix axis limits
    if(!is.logical(fix_axis[[i]])){
      p1 <- p1 + fix_axis[[i]]
    }
    
    if(!is.logical(add_lines[[i]])){
      p1 <- p1 + geom_hline(yintercept = add_lines[[i]])
    }

    legend = get_legend(p1 + theme(legend.position = 'bottom', 
                                   legend.text = element_text(size = 17),
                                   legend.title = element_text(size = 20),
                                   legend.key.size = unit(0.1, 'cm')))
    
    p1 <- p1 + theme(legend.position = 'none')
    
    plot_list[[i]] <- p1
  }
  
  ## ggpubr::ggarrange might be the way to go.
  if(include_legend){
    final_plot <- plot_grid(plot_grid(plotlist = plot_list, nrow = rows), legend, ncol = 1, rel_heights = c(10,1))
  }else{
    final_plot = plot_grid(plotlist = plot_list, nrow = rows)
  }
  
  if(!is.null(title)){
    title_plot <- ggplot() + 
      labs(title = title) + 
      theme_bw()
    final_plot <- plot_grid(title_plot, final_plot, ncol = 1, rel_heights = c(0.2 - 0.05*rows,1))
  }
  
  res_lst <- list(plot = final_plot, legend = legend)
  return(res_lst)
}

# plot the CAR parameter estimates from a set of simulations
# This function can either take in a list of files with results or a processed CAR results data set and plot them
plot_CAR_params <- function(files = NULL, res_df = NULL, expected = data.frame(param = c('tau2','rho.S','alpha'), expected_value = c(1, 0.3,0.3)), ...){
  if(!is.null(files)){
    res_df <- process_CAR_params_wrapper(files, ...)
  }
  
  tmp <- res_df %>% filter(param != 'betas')
  
  ## plot the parameter estimates
  if(!is.null(expected)){
    tmp = merge(tmp, expected)
    p1 <- ggplot(tmp) +
      facet_wrap(vars(param)) + 
      geom_point(aes(x = prop_missing, y = mean_param)) + 
      geom_errorbar(aes(x = prop_missing,  ymin = param_05, ymax = param_95), position = position_dodge(width = 0.1)) +
      geom_hline(aes(yintercept = expected_value), color = 'red') + 
      ylab('mean MCMC parameter estimate (5%, 95%) of simulations')
  }else{
    p1 <- ggplot(tmp) +
      facet_wrap(vars(param)) + 
      geom_point(aes(x = prop_missing, y = mean_param)) + 
      geom_errorbar(aes(x = prop_missing,  ymin = param_05, ymax = param_95), position = position_dodge(width = 0.1)) +
      ylab('mean param estimate (5%, 95%) of simulations')
  }
  
  ## plot the n effective results
  p2 <- ggplot(res_df) +
    facet_wrap(vars(param), nrow = 1) + 
    geom_point(aes(x = prop_missing, y = mean_n_eff)) + 
    geom_errorbar(aes(x = prop_missing,  ymin = n_eff_05, ymax = n_eff_95), position = position_dodge(width = 0.1)) +
    scale_y_continuous(trans='log2') + 
    ylab('mean n.effective (5%, 95%)')
  
  ## Combine into final plot
  final_plot <- cowplot::plot_grid(p1, p2, ncol = 1)
  
  return(final_plot)
}

#
##### Simulation functions #####
# function to simulate data
# simulate_imputation <- function(df, col, p = 0.1, group = NULL){
#   set.seed(1)
#   # if not worrying about the grouping of randomization
#   
#   # set the true value to null for the df
#   df[,paste0(col, '_true')] = as.integer(NA)
#   
#   if(is.null(group)){
#     # get the number of data points to impute
#     num_impute = round(p*nrow(df)) - sum(is.na(df[,col]))
#     
#     if(num_impute <= 0){
#       warning('skipping imputation. Already too many missing')
#       return(df)
#     }
#     
#     # sample the indices to impute
#     ind_sample = sample(setdiff(1:nrow(df), which(is.na(df[,col]))), num_impute)
#     
#     # store the true value of the values to replace
#     df[ind_sample, paste0(col, '_true')] <- df[ind_sample, col]
#     
#     # replace the values with NAs
#     df[ind_sample, col] = NA
#     
#   }else{
#     # doing grouping
#     uni_group = df %>% pull(UQ(group)) %>% unique()
#     
#     # apply to each group
#     tmp <- lapply(uni_group, function(xx) {
#       tt <- df %>% filter(get(group) == xx)
#       num_impute = round(p*nrow(tt)) - sum(is.na(tt[,col]))
#       if(num_impute <= 0){
#         print(sprintf('skipping imputation for %s', xx))
#         return(tt)
#       }
#       ind_sample = sample(setdiff(1:nrow(tt), which(is.na(tt[,col]))), num_impute) 
#       tt[ind_sample, paste0(col, '_true')] <- tt[ind_sample, col]
#       tt[ind_sample, col] = NA
#       return(tt)
#     })
#     
#     # collapse results
#     df <- data.table::rbindlist(tmp)
#   }
#   
#   return(df)
# }

initialize_df <- function(district_sizes, start_date = '2016-01-01', end_date = '2019-12-01', ...){
  facilities = unlist(lapply(1:length(district_sizes), function(xx) {paste0(toupper(letters[xx]), sprintf('%03d',1:district_sizes[xx]))}))
  
  dates = seq(as.Date(start_date), as.Date(end_date), by = 'month')
  
  df = expand.grid(facilities, dates, stringsAsFactors = T)
  colnames(df) = c('facility','date')
  
  df$district = substr(df$facility, 1, 1) 
  
  df = df %>%
    dplyr::select(date, facility, district) %>%
    arrange(facility, date)
  return(df)
}

sample_real_betas <- function(facilities, file = 'results/all_facility_betas_filtered_09052022.csv'){
  tmp <- read.csv(file)
  
  ind <- sample(nrow(tmp), length(facilities))
  betas <- tmp[ind,3:10]
  rownames(betas) = facilities
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
  return(betas)
}

sample_betas = function(facilities, b0_mean = 4.3, b1_mean = -0.25, b1_sd = 0.25, ...){
  betas = matrix(0, nrow = length(facilities), ncol = 8)
  
  if(length(b0_mean) == 1){
    betas[,1] = rnorm(mean = b0_mean, sd = 1, n = nrow(betas))
  }else{
    u = sample(b0_mean, nrow(betas), replace = T)
    betas[, 1] = rnorm(mean = u, sd = 1, n = nrow(betas))
  }
  
  betas[,2] = rnorm(b1_mean, b1_sd, n = nrow(betas))
  
  for(j in 3:8){
    betas[,j] = rnorm(0, 0.15, n = nrow(betas))
  }
  
  rownames(betas) = facilities
  colnames(betas) = c('intercept', 'year', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3')
  #paste0('B',0:7)
  #betas = cbind(facilities, as.data.frame(betas))
  return(betas)
}

### Function to simulate the data for a variety of situations
simulate_data <- function(district_sizes, R = 1, empirical_betas = F, seed = 10, type = 'WF', family = 'poisson', ...){
  
  # set seed so the betas are always the same (the seed input is used later for simulating the data on top of these betas)
  set.seed(10)
  
  # set up data frame
  df = initialize_df(district_sizes, ...)
  
  # make periodic covariates
  df = add_periodic_cov(df)
  
  # get all facility and district names
  facilities = unique(df$facility)
  districts = unique(df$district)
  
  # get all dates
  dates = unique(df$date) %>% sort()
  
  # sample betas
  if(empirical_betas){
    betas = sample_real_betas(facilities)
  }else{
    betas = sample_betas(facilities, ...)  
  }
  
  if(family == 'poisson'){
    DGP_function = rpois
  }else if(family == 'quasipoisson'){
    DGP_function <- function(n, mu){
      y <- rqpois(n, mu, theta = list(...)$theta)
    }
  }else if(family == 'negbin'){
    DGP_function <- function(n, mu){
      y <- rnbinom(n = n, mu = mu, size = list(...)$dispersion)
    }
  }
  
  # initialize list of data frames
  df_lst = list()
  district_lst = list()
  
  # set seed for the data generation
  set.seed(seed)
  
  #browser()
  # simulate the data according to the DGP
  if(type == 'WF'){
    # make R sampled sets of data
    for(i in 1:R){
      
      # simulate values given the betas
      tmp_lst = lapply(facilities, function(xx){
        tmp = df %>% filter(facility == xx)
        
        # keep the 1 for intercepts
        X = tmp %>% 
          mutate(intercept = 1) %>%
          dplyr::select(intercept, year, cos1, sin1, cos2, sin2, cos3, sin3)
        
        # error checking
        if(!identical(colnames(betas), colnames(X))){
          print(colnamaes(betas))
          print(colnames(X))
          stop(sprintf('colnames of betas not equal to X: %s, %s',paste(colnames(betas), collapse = ';'), paste(colnames(X), collapse = ';') ))
        }
        
        # make the 8x1 beta vector for this facility
        beta_f = t(betas[xx,,drop = F])
        
        # get mean prediction from linear model
        mu = as.matrix(X)%*%beta_f
        tmp$y_exp = exp(mu)[,1]
        
        # get Poisson or quasipoisson variance
        if(family == 'poisson'){
          tmp$y_var <- tmp$y_exp
        }else if(family == 'quasipoisson'){
          tmp$y_var <- list(...)$theta*tmp$y_exp
        }else if(family == 'negbin'){
          tmp$y_var <- tmp$y_exp + tmp$y_exp^2/list(...)$dispersion
        }else{
          stop('input a proper family')
        }
        
        # simluate random values
        tmp$y = DGP_function(length(mu), exp(mu))
        
        return(tmp)
      })
      
      # combine values into one data frame
      df = do.call('rbind', tmp_lst)
      df_lst[[i]] = df
      
      # group the results by district
      district <- df %>%
        group_by(district, date) %>%
        summarize(y_exp = sum(y_exp),
                  y_var = sum(y_var),
                  y = sum(y),
                  y_true = sum(y))
      district_lst[[i]] = district
      
    }

  }else if(type == 'CAR'){
    
    rho = list(...)$rho
    alpha = list(...)$alpha
    tau2 = list(...)$tau2
    
    # add in the mean effects because these are the same for all simulations
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
      tmp$mu = (as.matrix(X)%*%beta_f)[,1]
      
      return(tmp)
    }))
    
    # make the spatio-temporal precision matrix
    Q = make_precision_mat(df, rho = rho)
    
    tryCatch({
      V = tau2*solve(Q)
    }, error = function(e){
      print(e)
      print('havent dealt with non-invertible precision matrices yet')
      browser()
    })

    # checking ordering of facilities matches
    if(!identical(colnames(V), as.character(facilities))){
      browser()
      stop('names of covariances and facilities dont match')
    }
    
    # add in the marginal variance to original df
    dV = diag(V)
    matid = match(df$facility, names(dV))
    df$sigma2_marginal = dV[matid]
    
    # make R sampled sets of data
    df_lst = lapply(1:R, function(r){
      ### get the spatio-temporal random effects
      # initialize phi
      phi = matrix(0, nrow = length(dates), ncol = length(facilities))
      colnames(phi) = facilities
      
      # first time step
      phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
      
      # cycle through other time steps, using auto-correlated priors
      for(i in 2:length(dates)){
        phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
      }
      
      # convert to matching format
      phi = as.data.frame(phi)
      phi$date = dates
      phi_df = tidyr::gather(phi, facility, phi, setdiff(colnames(phi), c('facility','date'))) %>% 
        arrange(facility, date)
      
      # add in the previous phi values
      phi_df$phi_previous <- c(0,phi_df$phi[1:(nrow(phi_df) - 1)])
      phi_df$phi_previous[phi_df$date == min(phi_df$date)] <- 0
      
      # merge the phi values into the original data frame
      df = merge(df, phi_df, by = c('date','facility'))
      
      # calculate the expected value of the Poisson log-normal
      df$mu_marginal = df$mu + alpha*df$phi_previous
      df$y_exp = exp(df$mu_marginal + df$sigma2_marginal/2)
      
      # variance from Poisson or quasi-poisson log-normal
      if(family == 'poisson'){
        df$y_var = df$y_exp + (exp(df$sigma2_marginal) - 1)*exp(2*df$mu_marginal + df$sigma2_marginal)
      }else if(family == 'quasipoisson'){
        df$y_var = list(...)$theta*df$y_exp + (exp(df$sigma2_marginal) - 1)*exp(2*df$mu_marginal + df$sigma2_marginal)
      }else{
        stop('input a proper family')
      }
      
      # simulate the observed values
      df$y = DGP_function(nrow(df), exp(df$mu + df$phi))
      
      return(df)
    })
    
    # cycle through each created data frame
    district_lst = lapply(df_lst, function(df){
      # cycle through each district
      district <- do.call('rbind', lapply(districts, function(d){
        # filter data to this district
        df2 <- df %>% filter(district == d)
        facs = unique(df2$facility)
        
        # get covariance by this district
        V_d = V[facs, facs]
        V_exp = exp(V_d) - 1
        ind_cov = upper.tri(V_exp) + lower.tri(V_exp)
        
        # if(family == 'quasipoisson'){
        #   warning('cant do district level quasipoisson variance for CAR yet. Havent coded it.')
        # }
        
        district_df = do.call('rbind', lapply(1:length(dates), function(i_date){
          df3 <- df2 %>% filter(date == dates[i_date])
          covs = (df3$y_exp%*%t(df3$y_exp))*V_exp
          cov = sum(covs*ind_cov)
          
          tmp_df = data.frame(district = d,
                              date = dates[i_date],
                              y_exp = sum(df3$y_exp),
                              y_var_ind = sum(df3$y_var),
                              y_cov = cov,
                              y_var = sum(df3$y_var) + cov,
                              y = sum(df3$y),
                              y_true = sum(df3$y))
          return(tmp_df)
          }
        ))
        return(district_df)
      }))
      return(district)
    })
    
  }else if(type == 'freqGLM'){
    rho = list(...)$rho
    alpha = list(...)$alpha
    
    # get the adjacency matrix
    W <- make_district_adjacency(df, scale_by_num_neighbors = T)
    
    # add in the mean effects because these are the same for all simulations
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
      tmp$mu = (as.matrix(X)%*%beta_f)[,1]
      
      return(tmp)
    }))
    
    df$y_exp = df$y_var = df$y = NA
    
    # make R sampled sets of data
    df_lst = lapply(1:R, function(r){
      
      ### Get the first time point values
      ind = which(df$date == min(dates))
      
      # the adjustment accounts for the fact that there aren't additive auto-regressive and spatial terms at the first time point.
      # this adjustment comes from the sum of a geometric series (since this is roughly the effect that the spatial and autoregressive terms approach as we increase the time series)
      adjustment = 1 + rho/(1-rho) + alpha/(1-alpha)
      df$y_exp[ind] = df$y_var[ind] = adjustment*exp(df$mu[ind])
      df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
      
      # update the neighbors and auto-regressive terms
      df <- add_autoregressive(df, 'y') %>%
        add_neighbors(., 'y', scale_by_num_neighbors = T, W = W)
      
      ### remaining time points
      for(d in dates[-1]){
        # get the subset of dates
        ind = which(df$date == d)
        
        # get the mean estimates for these dates
        df$y_exp[ind] <- exp(df$mu[ind]) + 
          alpha*df$y.AR1[ind] + rho*df$y.neighbors[ind]
        
        # get Poisson or quasipoisson variance
        if(family == 'poisson'){
          df$y_var[ind] <-df$y_exp[ind]
        }else if(family == 'quasipoisson'){
          df$y_var[ind] <- list(...)$theta*df$y_exp[ind]
        }else{
          stop('input a proper family')
        }
        
        # predict!
        df$y[ind] = DGP_function(length(ind), df$y_exp[ind])
        
        # update the neighbors and auto-regressive terms
        df <- add_autoregressive(df, 'y') %>%
          add_neighbors(., 'y', scale_by_num_neighbors = T, W = W)
      }
      # ggplot(df, aes(x = date, y = y)) + 
      #   geom_line() + 
      #   facet_wrap(~facility)
      df
    })

    # group the results by district
    district_lst <- lapply(df_lst, function(df){
      district <- df %>%
        group_by(district, date) %>%
        summarize(y_exp = sum(y_exp),
                  y_var = sum(y_var),
                  y = sum(y),
                  y_true = sum(y))
      district
    })
  }else{
    stop('please input a proper type')
  }
  
  # make list of values to return
  res_lst = list(df_list = df_lst, district_list = district_lst, betas = betas)
  return(res_lst)
}

#### Function to simulate data under the freqGLM_epi framework
##
## district_sizes = number of facilities in the districts
## R = number of simulated data sets
## lambda = autoregressive term
## phi = neighbor term
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
      mean_y = exp(lambda)*tmp$y.AR1 + exp(phi)*tmp$y.neighbors + exp(tmp$mu_seasonal)
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
  params_true = rbind(t(data.frame(By.AR1 = rep(lambda, ncol(params_true)), By.neighbors = rep(phi, ncol(params_true)), row.names = colnames(params_true))), params_true)
  
  
  # return it!
  res_lst = list(df_list = df_lst, betas = betas, lambda = lambda, phi = phi, params = params_true)
  return(res_lst)
}

make_precision_mat <- function(df, rho, W2 = NULL){
  # check the rho value
  if(rho < 0 | rho >= 1){
    stop('please input a rho in [0,1)')
  }
  
  # get unique facilities
  facilities = unique(df$facility) %>% sort()
  
  # create the <W2 = diag(W1) - W> matrix
  if(is.null(W2)){
    W2 <- make_district_W2_matrix(df)
  }
  
  # make the final Q matrix
  Q = rho*W2 + (1-rho)*diag(rep(1, length(facilities)))
  
  return(Q)
}

MCAR_sim <- function(df, p, by_facility = F, max_missing_date = '2019-12-01'){
  # save the true y value
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)

  # delete y values to add in missingness
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  # combine the hold-out/non-missing data with the training data with missingness
  df <- rbind(df, df_test)
  
  return(df)
}

MNAR_sim <- function(df, p, direction = NULL, gamma = 1.5, by_facility = T, max_missing_date = '2019-12-01'){
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)
  
  if(by_facility){
    # get the number of points to impute
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    
    # if removing low and high points
    if(is.null(direction)){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = median(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    # if removing high points
    }else if(direction == 'upper'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = min(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    # if removing low points
    }else if(direction == 'lower'){
      df = do.call('rbind', lapply(unique(df$facility), function(xx){
        tmp = df %>% filter(facility == xx)
        m = max(tmp$y)
        q = (abs(tmp$y - m))^gamma
        tmp$y[sample(nrow(tmp), num_impute, prob = q)] <- NA
        return(tmp)
      }))
    }
    
  }else{
    print('havent coded this part - is it necessary?')
    browser()
    #num_impute = round(p*nrow(df))
    #df$y[sample(nrow(df), num_impute)] <- NA
  }
  
  df <- rbind(df[,colnames(df_test)], df_test)
  
  return(df)
}

MAR_spatiotemporal_sim <- function(df, p, rho = 0.3, alpha = 0.3, tau2 = 1, by_facility = T, max_missing_date = '2019-12-01'){
  # make phi for df
  # for all, or for each facility,
  # sample according to expit(phi)
  df$y_true = df$y
  
  # split the data into a test/hold-out set and the training set
  df_test <- df %>%
    filter(date > max_missing_date)
  df <- df %>%
    filter(date <= max_missing_date)
  
  # get all facility names
  facilities = as.character(unique(df$facility))
  
  # get all dates
  dates = unique(df$date) %>% sort()
  
  # make the spatio-temporal precision matrix
  Q = make_precision_mat(df, rho = rho)
  
  tryCatch({
    V = tau2*solve(Q)
  }, error = function(e){
    print(e)
    print('havent dealt with non-invertible precision matrices yet')
    browser()
  })
  
  # checking ordering of facilities matches
  if(!identical(colnames(V), facilities)){
    browser()
    stop('names of covariances and facilities dont match')
  }
  
  phi = matrix(0, nrow = length(dates), ncol = length(facilities))
  colnames(phi) = facilities
  
  # first time step
  phi[1,] = MASS::mvrnorm(n = 1, mu = rep(0, nrow(V)), Sigma = V)
  
  # cycle through other time steps, using auto-correlated priors
  for(i in 2:length(dates)){
    phi[i,] = MASS::mvrnorm(n = 1, mu = alpha*phi[i-1,], Sigma = V)
  }
  
  # convert to matching format
  phi = as.data.frame(phi)
  phi$date = dates
  phi_df = tidyr::gather(phi, facility, phiM, setdiff(colnames(phi), c('facility','date')))
  
  # convert to expit for probability
  phi_df$prob.sample = exp(phi_df$phiM)/(1 + exp(phi_df$phiM))
  
  # merge the phi values into the original data frame
  df = merge(df, phi_df, by = c('date','facility'))
  
  if(by_facility){
    num_impute = round(p*nrow(df)/length(unique(df$facility)))
    df = do.call('rbind', lapply(unique(df$facility), function(xx){
      tmp = df %>% filter(facility == xx)
      tmp$y[sample(nrow(tmp), num_impute, prob = tmp$prob.sample)] <- NA
      return(tmp)
    }))
  }else{
    num_impute = round(p*nrow(df))
    df$y[sample(nrow(df), num_impute, prob = df$prob.sample)] <- NA
  }
  
  # combine the hold-out/non-missing data with the training data with missingness
  df <- rbind(df[,colnames(df_test)], df_test)
  
  return(df)
}
