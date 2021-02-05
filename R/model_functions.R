## FUNCIONS FOR MODEL FITTING ##
library(dplyr)
library(tidyverse) 
library(readr) 
library(reshape2) 
library(lubridate)
library(lme4)
library(GLMMadaptive)
library(ciTools)



# SINGLE FACILITY COUNT & PROPORTION WITH PREDICTION INTERVAL
fit.site.specific.denom.pi <-function(data,site_name,extrapolation_date,
                                      indicator_var,site_var,date_var,denom_var=NULL,
                                      level="month",counts_only=FALSE,R=250){
  if(counts_only==FALSE & is.null(denom_var)){
    warning("You must supply a denominator variable")
  }
  
  # period 
  if (level == "month") {
    period = 12 
  } else if (level == "biweek"){
    period = 26 
  } else if (level == "week"){
    period = 52
  } else if (level == "day"){
    period= 365
  }
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
    data.new %>% filter(!is.na(indicator) & !is.na(denom)) %>% 
      filter(date < extrapolation_date) -> data.base 
    
  }
  # fit model (add 0.5 in case denom is zero)
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  formula_prop = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + offset(log(denom))"
  
  
  #fit unadjusted counts model
  tryCatch({
    
    ### First, filling in the indicator missing values and generating prediction intervals around the indicator.
    mod_counts <- glm(as.formula(formula_counts), data=data.base, family=quasipoisson)
    
    #check dispersion parameter
    overdisp <- summary(mod_counts)$dispersion
    
    # if there is overdispersion, keep the quasipoisson results. Otherwise, keep the poisson.
    if(overdisp > 1){
      
      # Nick comment: this is the same as doing predict(..., type = 'response')
      count.predict <- exp(predict(mod_counts,newdata=data.new))
      
    } else {
 
      mod_counts <- glm(as.formula(formula_counts), data=data.base, family=poisson)
      count.predict <- exp(predict(mod_counts,newdata=data.new))
      print(paste0("underdispersion detected for ",site_name," with dispersion ",overdisp))
    }
    #add_pi() function cannot have missing data - fill in with predicted value 
    data.new.filled <- data.new %>% mutate(indicator = ifelse(is.na(indicator),count.predict,indicator))
    
    #prediction interval - what to do with 
    pi <- suppressWarnings(add_pi(data.new.filled,mod_counts)) #warning is just about approximation
    pi.low <- pi$LPB0.025 
    pi.up <- pi$UPB0.975
    
    ### Running both count and proportion on fixed indicator variable
    if (counts_only==FALSE){
      #if denominator data given, calculate adjusted counts and proportion estimates
      tryCatch({
        # print(sum(is.na(data.new$denom))) <- not many. Max of 4. It would be interesting to mess around with this number to test missing data methods.
        
        # filling in missing denominator values
        if(sum(is.na(data.new$denom))>0){
          
          formula_denom = "denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
          mod_denom <- glm(as.formula(formula_denom), data=data.base, family=quasipoisson)
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
          
          # denom_no_miss is a matrix of the all the bootstrapped iterations filling in the missing data. So all the non-missing values (rows) will be identical but the previously missing values will have different predicted values now based off of the originally trained model and multivariate normal sampling from that model (aka parametric bootstrap).
          
        } else { 
          
          denom_no_miss <- matrix(data.new$denom, length(data.new$denom), R) 
          
        } # END UPDATE: add NA in denom
        
        # training without the denominator NAs filled in.
        mod_prop <- glm(as.formula(formula_prop), data=data.base, family=quasipoisson)
        
        # model doesn't return the offset - there's not beta for that. Hence the addition of the offset in the later part.
        beta_hat <- mod_prop$coefficients
        beta_vcov <- vcov(mod_prop)
        overdisp <- summary(mod_prop)$dispersion

        sapply(1:R, function(r){
          
          #indicator 
          beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
          pred_boot <- (data.new %>% 
                          mutate(intercept=1) %>%
                          dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                          as.matrix())%*%as.matrix(beta_boot)
          pred_boot_exp <- exp(pred_boot + log(denom_no_miss[,r])) # UPDATE
          pred_boot_exp[which(is.na(pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
          se_boot <- pred_boot_exp/(overdisp - 1)
          x = MASS::rnegbin(n = nrow(data.new),mu = pred_boot_exp, theta = se_boot)
          
          x/denom_no_miss[,r] # UPDATE
          
        }) -> sim.boot
        
        prop.pi.low <- apply(sim.boot,1,quantile,.025,na.rm=TRUE)
        prop.pi.up <- apply(sim.boot,1,quantile,.975,na.rm=TRUE)
        prop.predict <- apply(sim.boot,1,quantile,.50,na.rm=TRUE)
        
        # create data frame with predictions, 95% CIs, and observed values 
        results.df <<- data.frame(site=site_name,
                                  date=as.Date(data.new$date),
                                  est_raw_counts=count.predict,
                                  ci_raw_counts_low=pi.low,
                                  ci_raw_counts_up=pi.up,
                                  observed=data.new$indicator,
                                  est_prop=prop.predict,
                                  ci_low_prop=prop.pi.low,
                                  ci_up_prop=prop.pi.up,
                                  observed_prop=data.new$indicator/data.new$denom) 
        
      }, warning = function(w){
        # this code is used if site does not have sufficient data to fit 
        # this is expected as we are including all sites in initial analysis 
        # we filter out "bad" sites in the figure file
        
        results.df <<- data.frame(site = site_name,
                                  date=as.Date(unique(data.new$date)),
                                  est_raw_counts=NA,
                                  ci_raw_counts_low=NA,
                                  ci_raw_counts_up=NA,
                                  observed=data.new$indicator,
                                  est_prop=NA,
                                  ci_low_prop=NA,
                                  ci_up_prop=NA,
                                  observed_prop=data.new$indicator/data.new$denom)
      }, error = function(e){
        results.df <<- data.frame(site = site_name,
                                  date=as.Date(unique(data.new$date)),
                                  est_raw_counts=NA,
                                  ci_raw_counts_low=NA,
                                  ci_raw_counts_up=NA,
                                  observed=data.new$indicator,
                                  est_prop=NA,
                                  ci_low_prop=NA,
                                  ci_up_prop=NA,
                                  observed_prop=data.new$indicator/data.new$denom)
      }) 
    }
    else{
      
      results.df <<- data.frame(site=site_name,
                                date=as.Date(data.new$date),
                                est_raw_counts=count.predict,
                                ci_raw_counts_low=pi.low,
                                ci_raw_counts_up=pi.up,
                                observed=data.new$indicator,
                                est_prop=NA,
                                ci_low_prop=NA,
                                ci_up_prop=NA,
                                observed_prop=NA)
    }},warning = function(w){
      # this code is used if site does not have sufficient data to fit 
      # this is expected as we are including all sites in initial analysis 
      # we filter out "bad" sites in the figure file
      
      results.df <<- data.frame(site = site_name,
                                date=as.Date(unique(data.new$date)),
                                est_raw_counts=NA,
                                ci_raw_counts_low=NA,
                                ci_raw_counts_up=NA,
                                observed=data.new$indicator,
                                est_prop=NA,
                                ci_low_prop=NA,
                                ci_up_prop=NA,
                                observed_prop=NA)
    }, error = function(e){
      results.df <<- data.frame(site = site_name,
                                date=as.Date(unique(data.new$date)),
                                est_raw_counts=NA,
                                ci_raw_counts_low=NA,
                                ci_raw_counts_up=NA,
                                observed=data.new$indicator,
                                est_prop=NA,
                                ci_low_prop=NA,
                                ci_up_prop=NA,
                                observed_prop=NA)
    })
  
  return(results.df)
  
}



# CLUSTER PROPORTION AND COUNTS WITH PREDICTION INTERVALS (Bootstrap!)
# EXTRACT BETA COEFFICIENTS (and overdispersion)
# SINGLE FACILITY COUNT & PROPORTION WITH PREDICTION INTERVAL 
#
# why is R here? There is no bootstrapping in this function
site.model.coefficients <-function(data,
                                   site_name,
                                   extrapolation_date,
                                   indicator_var,
                                   denom_var=NA, # if proportions are desired
                                   site_var,
                                   date_var,
                                   level="month"){
  
  # period 
  if (level == "month") {
    period = 12 
  } else if (level == "biweek"){
    period = 26 
  } else if (level == "week"){
    period = 52
  } else if (level == "day"){
    period= 365
  }
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
  mod_counts <- glm(as.formula(formula_counts), data=data.base, family=quasipoisson)
  overdisp <- summary(mod_counts)$dispersion
  if(overdisp < 1){
    
    mod_counts <- glm(as.formula(formula_counts), data=data.base, family=poisson)
    print(paste0("underdispersion detected for ",site_name," with dispersion ",overdisp))
    
  }
  beta_hat <- mod_counts$coefficients
  beta_vcov <- vcov(mod_counts)
  predictions <- exp(predict(mod_counts,data.new))
  
  
  # PROPORTION MODEL (if applicable)
  if(!is.na(denom_var)){
    
    formula_prop = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + offset(log(denom))"
    mod_prop <- glm(as.formula(formula_prop), data=data.base, family=quasipoisson)
    overdisp_adj <- summary(mod_prop)$dispersion
    if(overdisp_adj < 1){
      mod_prop <- glm(as.formula(formula_prop), data=data.base, family=poisson)
      print(paste0("underdispersion detected for ",site_name," with dispersion ",overdisp))
    }
    
    beta_adj_hat <- mod_prop$coefficients
    beta_adj_vcov <- vcov(mod_prop)
    
  } else {
    
    beta_adj_hat <- NA
    beta_adj_vcov <- NA
    overdisp_adj <- NA
  }
  
  results.list <- list(model_matrix = data.new,
                       predictions = predictions,
                       beta_estimates=beta_hat,
                       beta_var=beta_vcov,
                       overdispersion=overdisp,
                       beta_estimates_adj=beta_adj_hat,
                       beta_var_adj=beta_adj_vcov,
                       overdispersion_adj=overdisp_adj)
  
  return(results.list)
}

# selects sites to keep based on missing criteria and minimum count
sites_included_cluster <- function(data,indicator,p_miss_eval,p_miss_base,n_count_base){
  
  data %>% 
    rename(observed = indicator) %>%
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(facility) %>% 
    dplyr::summarize(mean = mean(is.na(observed))) %>% 
    filter(mean <= p_miss_base) %>% dplyr::select(facility) %>% pull(facility) -> sites.base
  
  data %>% 
    rename(observed = indicator) %>%
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(facility) %>% 
    dplyr::summarize(median = median(observed,na.rm=TRUE)) %>% 
    filter(median > n_count_base) %>% dplyr::select(facility) %>% pull(facility) -> sites.sparse.base
  
  data %>% 
    rename(observed = indicator) %>%
    filter(date >= as.Date(extrapolation_date)) %>% 
    group_by(facility) %>% 
    dplyr::summarize(mean = mean(is.na(observed))) %>% 
    filter(mean <= p_miss_eval) %>% dplyr::select(facility) %>% pull(facility) -> sites.eval 
  
  sites.baseline <- intersect(sites.base,sites.sparse.base)
  sites.included <- intersect(sites.eval,sites.baseline)
  
  return(sites.included)
  
}

# For Debugging

# data = data
# indicator_var = "indicator_count_ari_total"
# denom_var = "indicator_denom"
# site_var = "facility"
# date_var = "date"
# counts_only=FALSE
# n_count_base = 0
# p_miss_base = 0.2
# p_miss_eval = 0.5
# R=250

# CLUSTER PROPORTION AND COUNTS WITH PREDICTION INTERVALS 
fit.cluster.pi <- function(data,
                           indicator_var,
                           denom_var,
                           date_var,
                           site_var,
                           counts_only=FALSE,
                           n_count_base,p_miss_base,p_miss_eval,
                           R=250){
  
  
  # Filtering out sites with not enough data
  site_list <- sites_included_cluster(data,indicator_var,p_miss_eval,p_miss_base,n_count_base)
  data %>% filter(facility %in% site_list) -> data_filtered
  
  # FORMAT DATA FRAME
  data_filtered %>% 
    dplyr::rename(indicator = indicator_var,
                  denom = denom_var,
                  date = date_var,
                  site = site_var) %>% 
    dplyr::select(date,site,indicator,denom) -> data.new

  
  facility_complete_list <- data.new %>% distinct(site) %>% pull()
  
  # MODEL FITS (COEFFICIENTS & VARIANCE FOR BOOTSTRAP)
  # denominator (if proportions)
  if(counts_only==FALSE){
    lapply(facility_complete_list,function(x){
      
      site.model.coefficients(data=data.new %>% dplyr::select(-indicator),
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="denom",
                              site_var="site",
                              date_var="date")
      
    }) -> facility_denom
    
    # indicator
    lapply(facility_complete_list,function(x){
      
      site.model.coefficients(data=data.new,
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="indicator",
                              denom_var="denom",
                              site_var="site",
                              date_var="date")
      
    }) -> facility_indicator
    
  } else {
    
    # indicator only
    lapply(facility_complete_list,function(x){
      
      site.model.coefficients(data=data.new,
                              site_name=x,
                              extrapolation_date=extrapolation_date,
                              indicator_var="indicator",
                              site_var="site",
                              date_var="date")
      
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
      overdisp <- facility_indicator[[z]]$overdispersion
      
      beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
      pred_boot <- (model_matrix %>% 
                      mutate(intercept=1) %>%
                      dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                      as.matrix())%*%as.matrix(beta_boot)
      pred_boot_exp <- exp(pred_boot)
      pred_boot_exp[which(is.na(pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
      
      # negative binomial or poisson?
      if (overdisp >= 1){
        se_boot <- pred_boot_exp/(overdisp - 1)
        x = MASS::rnegbin(n = nrow(model_matrix),mu = pred_boot_exp,theta = se_boot)
      } else {
        x = rpois(n = nrow(model_matrix),pred_boot_exp)
      }
      
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
        
        # fill in missing data
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
        overdisp_adj <- facility_indicator[[z]]$overdispersion_adj
        
        beta_boot <- MASS::mvrnorm(1,beta_hat_adj,beta_vcov_adj)
        pred_boot <- (model_matrix %>% 
                        mutate(intercept=1) %>%
                        dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                        as.matrix())%*%as.matrix(beta_boot)
        pred_boot_exp <- exp(pred_boot + log(denom_no_miss)) # UPDATE (with offset)
        pred_boot_exp[which(is.na(pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
        
        # negative binomial or poisson?
        if (overdisp_adj >= 1){
          se_boot <- pred_boot_exp/(overdisp_adj - 1)
          y = MASS::rnegbin(n = nrow(model_matrix),mu = pred_boot_exp,theta = se_boot)
        } else {
          y = rpois(n = nrow(model_matrix),pred_boot_exp)
        }
        
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
    pi.low <- apply(results.boot,1,quantile,.025,na.action=na.omit,na.rm=TRUE)
    pi.up <- apply(results.boot,1,quantile,.975,na.action=na.omit,na.rm=TRUE)
    
    pred.prop <- NA
    pi.low.prop <- NA
    pi.up.prop <- NA
    
  } else {
    
    results.boot <- do.call(cbind, lapply(results.boot.list, `[[`, 1))
    pred <- apply(results.boot,1,median)
    pi.low <- apply(results.boot,1,quantile,.025,na.action=na.omit,na.rm=TRUE)
    pi.up <- apply(results.boot,1,quantile,.975,na.action=na.omit,na.rm=TRUE)
    
    results.boot.prop <- do.call(cbind, lapply(results.boot.list, `[[`, 2))
    pred.prop <- apply(results.boot.prop,1,median)
    pi.low.prop <- apply(results.boot.prop,1,quantile,.025,na.action=na.omit,na.rm=TRUE)
    pi.up.prop <- apply(results.boot.prop,1,quantile,.975,na.action=na.omit,na.rm=TRUE)
    
  }
  
  data_predicted <- data.frame(date=unique(data.new$date),
                               est_raw_counts=pred,
                               ci_raw_counts_low=pi.low,
                               ci_raw_counts_up=pi.up,
                               est_prop=pred.prop,
                               ci_low_prop=pi.low.prop,
                               ci_up_prop=pi.up.prop)
  
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
    mutate(ci_low_prop = ifelse(ci_low_prop < 0, 0, ci_low_prop),
           ci_up_prop = ifelse(ci_up_prop > 1, 0, ci_up_prop),
           ci_raw_counts_low=ifelse(ci_raw_counts_low < 0,0,ci_raw_counts_low),
           ci_raw_counts_up=ifelse(ci_raw_counts_up < 0,0,ci_raw_counts_up)) -> results
  
  results <- results %>% mutate(site = data %>% distinct(county) %>% pull())
  
  return(results)
  
}





is.zero <- function(x) { as.numeric(x == 0)}
is.complete <- function(x,prop) { as.numeric(x <= prop)}


# return sites that have a missigness less than or equal to some threshold during eval period
observed <- function(data,indicator_var,date_var,site_var,prop_miss=0,site_list=TRUE){
  
  data %>% 
    dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
    filter(date >= as.Date(extrapolation_date)) %>% 
    group_by(site) %>% dplyr::summarize(mean = mean(is.na(indicator))) -> ds 
  
  
  if(site_list == TRUE){
    
    ds %>%  filter(mean <= prop_miss) %>% dplyr::select(site) %>% pull(site) -> site_list
    
    return(site_list)
    
  } else {
    
    ds %>% dplyr::summarize(sum=sum(is.complete(mean,prop=prop_miss)),
                            prop=mean(is.complete(mean,prop=prop_miss))) -> ds
    
    return(ds)
    
  }
}

# return sites that have a missingness less than or equal to some threshold
complete <- function(data,indicator_var,date_var,site_var,prop_miss=0.2,site_list=TRUE){
  
  data %>% 
    dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(site) %>% dplyr::summarize(mean = mean(is.na(indicator))) -> ds
  
  
  if(site_list == TRUE){
    
    ds %>%  filter(mean <= prop_miss) %>% dplyr::select(site) %>% pull(site) -> site_list
    
    return(site_list)
    
  } else {
    
    ds %>% dplyr::summarize(sum=sum(is.complete(mean,prop=prop_miss)),
                            prop=mean(is.complete(mean,prop=prop_miss))) -> ds
    
    return(ds)
    
  }
}

# return sites that have a median count above a threshold
sparse <- function(data,indicator_var,date_var,site_var,count=0,site_list=TRUE){
  
  data %>% 
    dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(site) %>% dplyr::summarize(median = median(indicator,na.rm=TRUE)) -> ds
  
  if(site_list == TRUE){
    
    ds %>%  filter(median > count) %>% dplyr::select(site) %>% pull(site) -> site_list
    
    return(site_list)
    
  } else {
    
    ds %>% dplyr::summarize(sum=sum(1-is.complete(median,prop=count)),
                            prop=mean(1-is.complete(median,prop=count))) -> ds
    
    return(ds)
    
  }
}



# return sites that have a missigness less than or equal to some threshold during eval period
site.inclusion <- function(data,
                           p_miss_eval=0.5,
                           p_miss_base=0.2,
                           n_count_base=0,
                           prop=FALSE,
                           site_var="site"){
  
  if(prop == TRUE){ 
    data %>% 
      dplyr::select(-observed) %>% dplyr::rename(observed = `observed_prop`) -> data
  }
  
  if(site_var != "site"){
    data %>% dplyr::rename(site = site_var) -> data
  }
  
  data %>% 
    filter(date >= as.Date(extrapolation_date)) %>% 
    group_by(site) %>% 
    dplyr::summarize(mean = mean(is.na(observed))) %>% 
    filter(mean <= p_miss_eval) %>% dplyr::select(site) %>% pull(site) -> sites.eval 
  
  data %>% 
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(site) %>% 
    dplyr::summarize(mean = mean(is.na(observed))) %>% 
    filter(mean <= p_miss_base) %>% dplyr::select(site) %>% pull(site) -> sites.base
  
  data %>% 
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(site) %>% 
    dplyr::summarize(median = median(observed,na.rm=TRUE)) %>% 
    filter(median > n_count_base) %>% dplyr::select(site) %>% pull(site) -> sites.sparse.base
  
  intersect(intersect(sites.base,sites.eval),sites.sparse.base) -> sites_predict
  
  return(sites_predict)
  
}

# return sites that have a missingness less than or equal to some threshol

incomplete_sites <- function(data,complete_sites_list,indicator_var,date_var,site_var){
  data[[paste0(`indicator_var`)]] %>% 
    dplyr::rename(site=site_var,date=date_var) %>%
    filter(!(site %in% complete_sites_list[[`indicator_var`]])) -> ds
  return(ds)
}

complete_sites <- function(data,complete_sites_list,indicator_var,date_var,site_var){
  data[[paste0(`indicator_var`)]] %>% 
    dplyr::rename(site=site_var,date=date_var) %>%
    filter((site %in% complete_sites_list[[`indicator_var`]])) -> ds
  return(ds)
}

median_counts <- function(data,indicator_var,date_var,site_var){
  
  data %>% 
    dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) %>%
    filter(date < as.Date(extrapolation_date)) %>% 
    group_by(site) %>% dplyr::summarize(med = median(indicator,na.rm=TRUE)) %>% 
    dplyr::summarize(med = median(med,na.rm=TRUE)) -> ds
  
  return(ds)
  
}
