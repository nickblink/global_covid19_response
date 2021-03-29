## FUNCIONS FOR MODEL FITTING ##
library(dplyr)
library(tidyverse) 
library(readr) 
library(reshape2) 
library(lubridate)
library(ciTools)
library(MASS)


###################################################
################# MAIN FUNCTIONS ##################
###################################################

# MAIN FUNCTION 1. FACILITY-LEVEL COUNTS AND PROPORTION
fit.site.specific.denom.pi <- function(data, # data frame with site_name, indicator_var, denom_var, date_var, site_var
                                      site_var, # name of health facility/site (if only one site, create site_var variable with one name)
                                      site_name, # name of site to apply modeling to
                                      extrapolation_date, # date that the evaluation period begins
                                      indicator_var, # name of indicator variable
                                      date_var, # name of date variable ("Date" format)
                                      denom_var=NULL, # name of denominator variable (if applicable)
                                      counts_only=FALSE, # TRUE if no denominator variable
                                      R=250){
  
  if(counts_only==FALSE & is.null(denom_var)){
    warning("You must supply a denominator variable")
  }
  
  # inputs 
  period = 12 #can be updated for daily, weekly, quarterly data
  
  # create data frame
  if (counts_only==TRUE){
    data %>% 
      dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
      filter(site==site_name) %>%
      dplyr::select(date,indicator) %>%
      arrange(date) %>% 
      dplyr::mutate(year = year(date) - min(year(date)) + 1,
                    total=1:n(),
                    cos1 = cos(2*1*pi*total/period),
                    sin1 = sin(2*1*pi*total/period),
                    cos2 = cos(2*2*pi*total/period),
                    sin2 = sin(2*2*pi*total/period),
                    cos3 = cos(2*3*pi*total/period),
                    sin3 = sin(2*3*pi*total/period)) -> data.new
    
    # remove NAs and dates after extrapolation date 
    data.new %>% filter(!is.na(indicator)) %>% 
      filter(date < extrapolation_date) -> data.base
  } else {
    data %>% 
      dplyr::rename(site=site_var,indicator=indicator_var,date=date_var,denom=denom_var) %>%
      filter(site==site_name) %>%
      dplyr::select(date,indicator,denom) %>%
      arrange(date) %>% 
      dplyr::mutate(year = year(date) - min(year(date)) + 1,
                    total=1:n(),
                    cos1 = cos(2*1*pi*total/period),
                    sin1 = sin(2*1*pi*total/period),
                    cos2 = cos(2*2*pi*total/period),
                    sin2 = sin(2*2*pi*total/period),
                    cos3 = cos(2*3*pi*total/period),
                    sin3 = sin(2*3*pi*total/period)) -> data.new
    
    #raw counts model 
    data.new %>% filter(!is.na(indicator)) %>% # UPDATE: remove no NA in denom
      filter(date < extrapolation_date) -> data.base 
    
  }
  
  #browser()
  # model specifications
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  formula_prop = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + offset(log(denom))"
  
  # STEP 1. Fit Counts Model
  tryCatch({
    
    mod_counts <- suppressWarnings(glm.nb(as.formula(formula_counts), data=data.base)) #warning if large theta value 
    
    #check overdispersion
    mod_counts_pois <- suppressWarnings(glm(as.formula(formula_counts), data=data.base, family=poisson)) # warning for non-integer values (possible with imputation)
    overdisp <- AER::dispersiontest(mod_counts_pois)$p.value
    
    if(overdisp < 0.05){
      
      count.predict <- exp(predict(mod_counts,newdata=data.new))
      pi.counts <- glm.nb.pi(mod_counts,R=R,prop=FALSE,data_all=data.new)
      
      print(paste0("Overdispersion detected for ",site_name,". Negative binomial will be used."))
      
    } else {
      
      mod_counts <- mod_counts_pois
      count.predict <- exp(predict(mod_counts,newdata=data.new))
      
      data.new.filled <- data.new %>% mutate(indicator = ifelse(is.na(indicator),count.predict,indicator))
      pi <- suppressWarnings(add_pi(data.new.filled,mod_counts)) #warning is just about approximation
      pi.low <- pi$LPB0.025 
      pi.up <- pi$UPB0.975
      
      pi.counts <- data.frame(prediction=pi,pi.low=pi.low,pi.up=pi.up)
    }
    
    # STEP 2. Fit Proportion Model (if selected / denominator available)
    if (counts_only==FALSE){
      
      tryCatch({
        
        # Impute missing denominator values
        if(sum(is.na(data.new$denom))>0){
          
          formula_denom = "denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
          mod_denom <- suppressWarnings(glm.nb(as.formula(formula_denom), data=data.base))
          alpha_hat <- mod_denom$coefficients
          alpha_vcov <- vcov(mod_denom)
          
          sapply(1:R, function(r) {
            set.seed(10*r)
            alpha_boot <- MASS::mvrnorm(1,alpha_hat,alpha_vcov)
            
            pred_boot_denom <- (data.new %>% 
                                  mutate(intercept=1) %>%
                                  dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                                  as.matrix())%*%as.matrix(alpha_boot)
            pred_boot_denom_exp <- exp(pred_boot_denom)
            
            denom_updated <- data.new$denom
            denom_updated[is.na(data.new$denom)] <- pred_boot_denom_exp[is.na(data.new$denom)]
            
            denom_updated
            
          }) -> denom_no_miss
          
        } else { 
          
          denom_no_miss <- matrix(data.new$denom, length(data.new$denom), R) 
          
        } 
        
        # Fit model and calculate prediction intervals
        mod_prop <- suppressWarnings(glm.nb(as.formula(formula_prop),data=data.base))
        
        # check for overdispersion
        mod_counts_pois_prop <- glm(as.formula(formula_prop), data=data.base, family=poisson)
        overdisp_adj <- AER::dispersiontest(mod_counts_pois)$p.value
        
        if(overdisp_adj < 0.05){
          
          pi.prop <- glm.nb.pi(mod_prop,R=R,prop=TRUE,data_all=data.new,denom_values=denom_no_miss)
          
        } else {
          
          pi.prop <- glm.pois.prop.pi(mod_prop,R=R,data_all=data.new,denom_values=denom_no_miss)
          
        }
        
        # STEP 3. Compile Results
        # create data frame with predictions, 95% CIs, and observed values 
        results.df <- data.frame(site=site_name, #ADD < back in
                                 date=as.Date(data.new$date),
                                 observed=data.new$indicator,
                                 est_count=count.predict,
                                 ci_count_low=pi.counts$pi.low,
                                 ci_count_up=pi.counts$pi.up,
                                 observed_prop=data.new$indicator/data.new$denom,
                                 est_prop=pi.prop$prediction,
                                 ci_prop_low=pi.prop$pi.low,
                                 ci_prop_up=pi.prop$pi.up) 
        
      }, warning = function(w){
        # this code is used if site does not have sufficient data to fit
        # this is expected as we are including all sites in initial analysis
        # we filter out "bad" sites in the figure file
        
        results.df <<- data.frame(site = site_name,
                                  date=as.Date(unique(data.new$date)),
                                  observed=data.new$indicator,
                                  est_count=NA,
                                  ci_count_low=NA,
                                  ci_count_up=NA,
                                  observed_prop=data.new$indicator/data.new$denom,
                                  est_prop=NA,
                                  ci_prop_low=NA,
                                  ci_prop_up=NA)
      }, error = function(e){
        results.df <<- data.frame(site = site_name,
                                  date=as.Date(unique(data.new$date)),
                                  observed=data.new$indicator,
                                  est_count=NA,
                                  ci_count_low=NA,
                                  ci_count_up=NA,
                                  observed_prop=data.new$indicator/data.new$denom,
                                  est_prop=NA,
                                  ci_prop_low=NA,
                                  ci_prop_up=NA)
      })
    }
    else{
      
      results.df <<- data.frame(site=site_name,
                                date=as.Date(data.new$date),
                                est_count=count.predict,
                                ci_count_low=pi.counts$pi.low,
                                ci_count_up=pi.counts$pi.up,
                                observed=data.new$indicator,
                                observed_prop=NA,
                                est_prop=NA,
                                ci_prop_low=NA,
                                ci_prop_up=NA)
    }},warning = function(w){
      # this code is used if site does not have sufficient data to fit
      # this is expected as we are including all sites in initial analysis
      # we filter out "bad" sites in the figure file
      print(sprintf('warning for site %s: %s', site_name, w))
      
      results.df <<- data.frame(site = site_name,
                                date=as.Date(unique(data.new$date)),
                                observed=data.new$indicator,
                                est_count=NA,
                                ci_count_low=NA,
                                ci_count_up=NA,
                                observed_prop=NA,
                                est_prop=NA,
                                ci_prop_low=NA,
                                ci_prop_up=NA)
    }, error = function(e){
      
      print(sprintf('error for site %s: %s', site_name, e))
      
      results.df <<- data.frame(site = site_name,
                                date=as.Date(unique(data.new$date)),
                                observed=data.new$indicator,
                                est_count=NA,
                                ci_count_low=NA,
                                ci_count_up=NA,
                                observed_prop=NA,
                                est_prop=NA,
                                ci_prop_low=NA,
                                ci_prop_up=NA)
    })
  
  return(results.df)
  
}


# MAIN FUNCTION 2. AGGREGATED COUNTS AND PROPORTION 
fit.aggregate.pi.boot <- function(data, # data frame with indicator_var, denom_var, date_var, site_var
                                  indicator_var, # name of indicator variable
                                  denom_var, # name of denominator variable (if not applicable, can use any name)
                                  date_var, # name of date variable
                                  site_var, # name of smaller geographic unit variable (e.g. health facility)
                                  counts_only=FALSE, # TRUE if no denominator variable
                                  R=250){
  
  # INPUT
  period = 12 #can be updated for daily, weekly, quarterly data
  
  # FORMAT DATA FRAME
  data %>% 
    dplyr::rename(indicator = indicator_var,
                  denom = denom_var,
                  date = date_var,
                  site = site_var) %>% 
    dplyr::select(date,site,indicator,denom) -> data.new
  
  facility_complete_list <- data.new %>% distinct(site) %>% pull()
  
  #browser()
  
  # MODEL FITS (COEFFICIENTS & VARIANCE FOR BOOTSTRAP)
  # denominator (if proportions)
  if(counts_only==FALSE){
    lapply(facility_complete_list,function(x){
      print(x)
      site.model.coefficients(data=data.new %>% dplyr::select(-indicator),
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="denom",
                              site_var="site",
                              date_var="date",
                              R=R)
      
    }) -> facility_denom
    
    # indicator
    lapply(facility_complete_list,function(x){
      
      site.model.coefficients(data=data.new,
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="indicator",
                              denom_var="denom",
                              site_var="site",
                              date_var="date",
                              R=R)
      
    }) -> facility_indicator
    
  } else {
    
    # indicator only
    lapply(facility_complete_list,function(x){
      
      site.model.coefficients(data=data.new,
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="indicator",
                              site_var="site",
                              date_var="date",
                              R=R)
      
    }) -> facility_indicator
    
  }
  
  # BOOTSTRAP  
  lapply(1:R, function(r){
    set.seed(10*r)
    
    # counts
    sapply(1:length(facility_complete_list), function(z){
      
      model_matrix <- facility_indicator[[z]]$model_matrix
      beta_hat <- facility_indicator[[z]]$beta_estimates
      beta_vcov <- facility_indicator[[z]]$beta_var
      theta <- facility_indicator[[z]]$theta
      
      beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
      pred_boot <- (model_matrix %>% 
                      mutate(intercept=1) %>%
                      dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                      as.matrix())%*%as.matrix(beta_boot)
      pred_boot_exp <- exp(pred_boot)
      pred_boot_exp[which(is.na(pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
      
      x = MASS::rnegbin(n = nrow(model_matrix),mu = pred_boot_exp,theta = theta)
      
      x
      
    }) -> counts.boot
    
    sum.counts.boot <- apply(counts.boot,1,sum)
    
    # proportions (if applicable)
    if(counts_only==FALSE){
      
      lapply(1:length(facility_complete_list), function(z){
        
        # are there missing denominator values for this facility?
        data.new %>% 
          filter(site == facility_complete_list[z]) %>%
          dplyr::summarize(na_denom = sum(is.na(denom)))%>% 
          pull(na_denom) -> na_denom
        
        if(na_denom>0){
          
          model_matrix <- facility_denom[[z]]$model_matrix
          alpha_hat <- facility_denom[[z]]$beta_estimates
          alpha_vcov <- facility_denom[[z]]$beta_var
          
          alpha_boot <- MASS::mvrnorm(1,alpha_hat,alpha_vcov)
          
          pred_boot_denom <- (model_matrix %>% 
                                mutate(intercept=1) %>%
                                dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                                as.matrix())%*%as.matrix(alpha_boot)
          pred_boot_denom_exp <- exp(pred_boot_denom)
          
          denom_updated <- model_matrix$indicator
          denom_updated[is.na(model_matrix$indicator)] <- pred_boot_denom_exp[is.na(model_matrix$indicator)]
          
          denom_no_miss <- as.matrix(denom_updated)
          
        } else { 
          
          denom_no_miss <- data.new %>% filter(site == facility_complete_list[z]) %>% pull(denom) %>% as.matrix()
          
        } 
        
        # bootstrap step
        model_matrix <- facility_indicator[[z]]$model_matrix
        beta_hat_adj <- facility_indicator[[z]]$beta_estimates_adj
        beta_vcov_adj <- facility_indicator[[z]]$beta_var_adj
        theta_adj <- facility_indicator[[z]]$theta_adj
        
        beta_boot <- MASS::mvrnorm(1,beta_hat_adj,beta_vcov_adj)
        pred_boot <- (model_matrix %>% 
                        mutate(intercept=1) %>%
                        dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                        as.matrix())%*%as.matrix(beta_boot)
        pred_boot_exp <- exp(pred_boot + log(denom_no_miss)) # UPDATE
        pred_boot_exp[which(is.na(pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
        
        y = MASS::rnegbin(n = nrow(model_matrix),mu = pred_boot_exp,theta = theta_adj)
        
        cbind(as.matrix(y)/denom_no_miss,denom_no_miss)
        
      }) -> prop.boot
      
      # total number of outpatient visits each month
      N_r <- apply(do.call(cbind,lapply(1:length(facility_complete_list), function(z) prop.boot[[z]][,2])),1,sum)
      
      # weighted average of proportions (based on facility size)
      mean.prop.boot.list <- lapply(1:length(facility_complete_list), function(z) prop.boot[[z]][,1]*(prop.boot[[z]][,2]/N_r)) 
      mean.prop.boot <- apply(do.call(cbind,mean.prop.boot.list),1,sum)
      
      list(sum.counts.boot,mean.prop.boot)
      
    } else {
      
      sum.counts.boot
      
    }
  }) -> results.boot.list
  
  # collate results
  if(counts_only==TRUE){
    
    results.boot <- do.call(cbind,results.boot.list)
    
    pred <- apply(results.boot,1,median)
    pi.low <- apply(results.boot,1,quantile,.025)
    pi.up <- apply(results.boot,1,quantile,.975)
    
    pred.prop <- NA
    pi.low.prop <- NA
    pi.up.prop <- NA
    
  } else {
    
    results.boot <- do.call(cbind, lapply(results.boot.list, `[[`, 1))
    pred <- apply(results.boot,1,median)
    pi.low <- apply(results.boot,1,quantile,.025)
    pi.up <- apply(results.boot,1,quantile,.975)
    
    results.boot.prop <- do.call(cbind, lapply(results.boot.list, `[[`, 2))
    pred.prop <- apply(results.boot.prop,1,median)
    pi.low.prop <- apply(results.boot.prop,1,quantile,.025)
    pi.up.prop <- apply(results.boot.prop,1,quantile,.975)
    
  }
  
  data_predicted <- data.frame(date=unique(data.new$date),
                               est_count=pred,
                               ci_count_low=pi.low,
                               ci_count_up=pi.up,
                               est_prop=pred.prop,
                               ci_prop_low=pi.low.prop,
                               ci_prop_up=pi.up.prop)
  
  # Calculate observed (fill in missing with predicted values)
  if(counts_only==FALSE){
    
    do.call(rbind,lapply(1:length(facility_complete_list), function(z) {
      
      data.frame(date=unique(data.new$date),
                 site=facility_complete_list[z],
                 observed_indicator=facility_indicator[[z]]$model_matrix$indicator,
                 observed_denominator=facility_denom[[z]]$model_matrix$indicator,
                 pred=facility_indicator[[z]]$predictions,
                 pred_denominator=facility_denom[[z]]$predictions) %>%
        mutate(observed_indicator_new = case_when(is.na(observed_indicator) ~ pred,
                                                  TRUE ~ observed_indicator)) %>%
        mutate(observed_denominator_new = case_when(is.na(observed_denominator) ~ pred_denominator,
                                                    TRUE ~ observed_denominator)) %>% 
        dplyr::select(date,site,observed_indicator_new,observed_denominator_new) 
      
    })) -> observed
    
    
  } else {
    
    do.call(rbind,lapply(1:length(facility_complete_list), function(z) {
      
      data.frame(date=unique(data.new$date),
                 site=facility_complete_list[z],
                 observed_indicator=facility_indicator[[z]]$model_matrix$indicator,
                 pred=facility_indicator[[z]]$predictions) %>%
        mutate(observed_indicator_new = case_when(is.na(observed_indicator) ~ pred,
                                                  TRUE ~ observed_indicator)) %>%
        dplyr::select(date,site,observed_indicator_new) %>% 
        mutate(observed_denominator_new=NA)
      
    })) -> observed
    
  }
  
  observed %>% 
    group_by(date) %>% 
    dplyr::summarize(observed_indicator_sum = sum(observed_indicator_new),
                     observed_denominator_sum = sum(observed_denominator_new)) %>% 
    mutate(observed = observed_indicator_sum,
           observed_prop = observed_indicator_sum/observed_denominator_sum) %>% 
    dplyr::select(date,observed,observed_prop) -> data_observed
  
  
  # final data frame
  left_join(data_predicted,data_observed,by="date") %>%
    dplyr::select(date,observed,est_count,ci_count_low,ci_count_up,
                  observed_prop,est_prop,ci_prop_low,ci_prop_up) -> results
  
  return(results)
  
}


###################################################
##### HELPER FUNCTIONS FOR THE MAIN FUNCTIONS #####
###################################################

# Calculate prediction interval from a Negative Binomial model (counts or proportions) 
glm.nb.pi <- function(fit,R,prop,data_all,denom_values=1){
  
  if (prop == FALSE){
    denom_values <- matrix(1, nrow(data_all), R) 
  } 
  
  # extract information from model fit
  overdisp <- summary(fit)$theta
  beta_hat <- fit$coefficients
  beta_vcov <- vcov(fit)
  
  # bootstrap 
  sapply(1:R, function(r){
    beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
    pred_boot <- (data_all %>% 
                    mutate(intercept=1) %>%
                    dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                    as.matrix())%*%as.matrix(beta_boot)
    pred_boot_exp <- exp(pred_boot + log(denom_values[,r])) 
    pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
    x = MASS::rnegbin(n = nrow(data_all),mu = pred_boot_exp,theta = overdisp)
    
    x/denom_values[,r]
    
  }) -> sim.boot
  
  # create data frame
  pi.low <- apply(sim.boot,1,quantile,.025,na.rm=TRUE)
  pi.up <- apply(sim.boot,1,quantile,.975,na.rm=TRUE)
  predict <- apply(sim.boot,1,quantile,.50,na.rm=TRUE)
  
  return(data.frame(prediction=predict,
                    pi.low=pi.low,
                    pi.up=pi.up))
  
}

# Calculate prediction interval from a Poisson model with imputed denom values
glm.pois.prop.pi <- function(fit,R,data_all,denom_values){
  
  # extract information from model fit
  beta_hat <- fit$coefficients
  beta_vcov <- vcov(fit)
  
  # bootstrap 
  sapply(1:R, function(r){
    beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
    pred_boot <- (data_all %>% 
                    mutate(intercept=1) %>%
                    dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                    as.matrix())%*%as.matrix(beta_boot)
    pred_boot_exp <- exp(pred_boot + log(denom_values[,r])) 
    pred_boot_exp[which(is.na(pred_boot_exp))] <- 1
    x = MASS::rpois(n = nrow(data_all),pred_boot_exp)
    
    x/denom_values[,r]
    
  }) -> sim.boot
  
  # create data frame
  pi.low <- apply(sim.boot,1,quantile,.025,na.rm=TRUE)
  pi.up <- apply(sim.boot,1,quantile,.975,na.rm=TRUE)
  predict <- apply(sim.boot,1,quantile,.50,na.rm=TRUE)
  
  return(data.frame(prediction=predict,
                    pi.low=pi.low,
                    pi.up=pi.up))
  
}

# EXTRACT BETA COEFFICIENTS FROM SINGLE FACILITY COUNT & PROPORTION WITH PREDICTION INTERVAL 
site.model.coefficients <-function(data,
                                   site_name,
                                   extrapolation_date,
                                   indicator_var,
                                   denom_var=NA, # if proportions are desired
                                   site_var,
                                   date_var,
                                   R=250){
  
  # Inputs
  period = 12 #can be updated for daily, weekly, quarterly data

  # create data frame
  if (is.na(denom_var)){
    data %>% 
      dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
      filter(site==site_name) %>%
      dplyr::select(date,indicator) %>%
      arrange(date) %>% 
      dplyr::mutate(year = year(date) - min(year(date)) + 1,
                    total=1:n(),
                    cos1 = cos(2*1*pi*total/period),
                    sin1 = sin(2*1*pi*total/period),
                    cos2 = cos(2*2*pi*total/period),
                    sin2 = sin(2*2*pi*total/period),
                    cos3 = cos(2*3*pi*total/period),
                    sin3 = sin(2*3*pi*total/period)) -> data.new
    
    # remove NAs and dates after extrapolation date 
    data.new %>% filter(!is.na(indicator)) %>% 
      filter(date < extrapolation_date) -> data.base
  } else {
    data %>% 
      dplyr::rename(site=site_var,indicator=indicator_var,date=date_var,denom=denom_var) %>%
      filter(site==site_name) %>%
      dplyr::select(date,indicator,denom) %>%
      arrange(date) %>% 
      dplyr::mutate(year = year(date) - min(year(date)) + 1,
                    total=1:n(),
                    cos1 = cos(2*1*pi*total/period),
                    sin1 = sin(2*1*pi*total/period),
                    cos2 = cos(2*2*pi*total/period),
                    sin2 = sin(2*2*pi*total/period),
                    cos3 = cos(2*3*pi*total/period),
                    sin3 = sin(2*3*pi*total/period)) -> data.new
    
    #raw counts model 
    data.new %>% filter(!is.na(indicator)) %>% # UPDATE: remove no NA in denom
      filter(date < extrapolation_date) -> data.base 
    
  }
  
  # COUNT MODEL
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  mod_counts <- glm.nb(as.formula(formula_counts), data=data.base)
  
  theta <- mod_counts$theta
  beta_hat <- mod_counts$coefficients
  beta_vcov <- vcov(mod_counts)
  predictions <- exp(predict(mod_counts,data.new))
  
  
  # PROPORTION MODEL (if applicable)
  if(!is.na(denom_var)){
    
    formula_prop = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + offset(log(denom))"
    mod_prop <- glm.nb(as.formula(formula_prop), data=data.base)
    theta_adj <- mod_prop$theta
    
    beta_adj_hat <- mod_prop$coefficients
    beta_adj_vcov <- vcov(mod_prop)
    
  } else {
    
    beta_adj_hat <- NA
    beta_adj_vcov <- NA
    theta_adj <- NA
    mod_prop <- NA
  }
  
  results.list <- list(model_matrix = data.new,
                       predictions = predictions,
                       beta_estimates=beta_hat,
                       beta_var=beta_vcov,
                       theta=theta,
                       beta_estimates_adj=beta_adj_hat,
                       beta_var_adj=beta_adj_vcov,
                       theta_adj=theta_adj,
                       fitted.count=mod_counts,
                       fitted.prop=mod_prop)
  
  return(results.list)
}

