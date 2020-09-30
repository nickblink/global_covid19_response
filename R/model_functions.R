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
    data.new %>% filter(!is.na(indicator)) %>% 
      filter(date < extrapolation_date) -> data.base 
    
    # remove NAs and dates after extrapolation date 
    data.new %>% filter(!is.na(indicator)) %>% 
      filter(!is.na(denom)) %>% 
      filter(date < extrapolation_date) %>% 
      dplyr::select(-indicator) -> data.base.denom
    
    data.new %>% dplyr::select(-indicator) -> data.new.denom
  }
  # fit model (add 0.5 in case denom is zero)
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  formula_denom = "denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  
  
  #fit unadjusted counts model
  tryCatch({
    
    mod_counts <- glm(as.formula(formula_counts), data=data.base, family=quasipoisson)
    count.predict <- exp(predict(mod_counts,newdata=data.new))
    
    #add_pi() function cannot have missing data - fill in with predicted value 
    data.new.filled <- data.new %>% mutate(indicator = ifelse(is.na(indicator),count.predict,indicator))
    
    #prediction interval 
    pi <- suppressWarnings(add_pi(data.new.filled,mod_counts)) #warning is just about approximation
    pi.low <- pi$LPB0.025 
    pi.up <- pi$UPB0.975
    
    if (counts_only==FALSE){
      #if denominator data given, calculate adjusted counts and proportion estimates
      tryCatch({
        
        #model for denom
        mod_counts_denom <- suppressWarnings(glm(as.formula(formula_denom), data=data.base.denom, family=quasipoisson))
        count.predict.denom <- exp(predict(mod_counts_denom,newdata=data.new.denom))
        
        #add_pi() function cannot have missing data - fill in with predicted value 
        data.new.denom.filled <- data.new.denom %>% mutate(denom = ifelse(is.na(denom),count.predict.denom,denom))
        
        #prediction interval 
        pi.denom <- suppressWarnings(add_pi(data.new.denom.filled,mod_counts_denom)) #warning is just about approximation
        pi.low.denom <- pi.denom$LPB0.025 
        pi.up.denom <- pi.denom$UPB0.975
        
        #get empirical prediction interval for proportion using parametric bootstrap 
        #TO DO: double check pred interval is approx normal (sample from betas instead?)
        
        se.indicator <- (pi.up-pi.low)/(2*qnorm(.975))
        # se.denom <- (pi.up.denom-pi.low.denom)/(2*qnorm(.975))
        prop.boot <- sapply(1:R, function(x){
          boot.indicator <- rnorm(nrow(data.new),count.predict,se.indicator)
          # boot.denom <- rnorm(nrow(data.new),count.predict.denom,se.denom)
          boot.denom <- data.new.denom %>% 
            mutate(se_denom = (as.numeric(pi.up.denom)-as.numeric(pi.low.denom))/(2*qnorm(.975))) %>%
            mutate(boot_value = case_when(is.na(denom) ~ rnorm(1,count.predict.denom,se_denom),
                                          TRUE ~ as.numeric(denom))) %>% pull(boot_value)
          boot.indicator/boot.denom
        })
        
        prop.predict <- apply(prop.boot,1,quantile,.50)
        
        prop.pi.low <- apply(prop.boot,1,quantile,.025)
        prop.pi.low[which(prop.pi.low<0)] <- 0
        
        prop.pi.up <- apply(prop.boot,1,quantile,.975)
        prop.pi.up[which(prop.pi.up>1)] <- 1
        
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
fit.cluster.pi <- function(data,
                           indicator_var,
                           denom_var,
                           date_var,
                           site_var,
                           denom_results_all, # list
                           indicator_results_all, # list
                           counts_only=FALSE,
                           n_count_base,p_miss_base,p_miss_eval,
                           R=250){
  
  # filter & arrange date in both data sources
  data %>% arrange(date) -> data.new
  
  indicator_results_all[[`indicator_var`]] %>% 
    filter(site %in% unique(data.new$facility)) %>% arrange(date) -> indicator_results
  
  # remove sites with sparse or missing indicator and/or denom data
  sparse(data,indicator_var=indicator_var,date_var="date",site_var="facility",count=n_count_base) -> site_notsparse
  observed(data,indicator_var=indicator_var,date_var="date",site_var="facility",prop_miss=p_miss_eval) -> site_observed
  complete(data,indicator_var=indicator_var,date_var="date",site_var="facility",prop_miss=p_miss_base) -> site_complete
  indicator_results %>% filter(!is.na(est_raw_counts)) %>% distinct(site) %>% pull() -> results_notNA
  
  
  if (counts_only==FALSE){  
    denom_results_all[[`denom_var`]] %>%
      filter(site %in% unique(data.new$facility)) %>% arrange(date) -> denom_results
    
    sparse(data,indicator_var=denom_var,date_var="date",site_var="facility",count=n_count_base) -> site_notsparse_denom
    observed(data,indicator_var=denom_var,date_var="date",site_var="facility",prop_miss=p_miss_eval) -> site_observed_denom
    complete(data,indicator_var=denom_var,date_var="date",site_var="facility",prop_miss=p_miss_base) -> site_complete_denom
    denom_results %>% filter(!is.na(est_raw_counts)) %>% distinct(site) %>% pull() -> results_notNA_denom
    
    facility_complete_list <- Reduce(intersect, list(site_observed,site_complete,site_notsparse,results_notNA,
                                                     site_observed_denom,site_complete_denom,site_notsparse_denom,results_notNA_denom))
    
  } else { 
    
    facility_complete_list <- Reduce(intersect, list(site_observed,site_complete,site_notsparse,results_notNA))
    
    
  }
  
  p_facility_excluded <- 1-length(facility_complete_list)/length(unique(data.new$facility))
  facility_excluded <- unique(data.new$facility)[which(!unique(data.new$facility) %in% facility_complete_list)]
  
  if(length(facility_complete_list)==0){ #exit function if no facilities are included 
    
    data.frame(date=as.Date(unique(data.new$date)),
               est_raw_counts=NA,
               ci_raw_counts_low=NA,
               ci_raw_counts_up=NA,
               observed=NA,
               est_prop=NA,
               ci_low_prop=NA,
               ci_up_prop=NA,
               observed_prop=NA) -> df
    
    return(df)
    
  }
  
  data %>% 
    dplyr::rename(indicator = indicator_var,
                  site = site_var,
                  date = date_var) %>% 
    filter(facility %in% facility_complete_list) -> data.new
  
  
  # start bootstrap 
  lapply(1:R, function(r){
    
    # counts
    sapply(facility_complete_list, function(z){
      
      # subset data to site
      indicator_results2 <- indicator_results %>% filter(site == z)
      
      indicator_results2 %>% 
        mutate(se_indicator = (as.numeric(ci_raw_counts_up)-as.numeric(ci_raw_counts_low))/(2*qnorm(.975))) %>%
        mutate(boot_value = rnorm(nrow(indicator_results2),est_raw_counts,se_indicator)) %>% pull(boot_value)
      
      
    }) -> counts.boot
    
    sum.counts.boot <- apply(counts.boot,1,sum)
    
    # proportions (if applicable)
    if(counts_only==FALSE){
      
      sapply(facility_complete_list, function(z){
        
        # subset data to site
        denom_results2 <- denom_results %>% filter(site == z)
        
        # draw from normal approx (but only for those with missing values - this is "fixed")
        denom_results2 %>% 
          mutate(se_denom = (as.numeric(ci_raw_counts_up)-as.numeric(ci_raw_counts_low))/(2*qnorm(.975))) %>%
          mutate(boot_value = case_when(is.na(observed) ~ rnorm(1,est_raw_counts,se_denom),
                                        TRUE ~ as.numeric(observed))) %>% pull(boot_value)
        
      }) -> denom.boot
      
      sum.denom.boot<- apply(denom.boot,1,sum) #TO DO: does this need to be embedded in the above bootstrap or no?
      prop.boot <- sum.counts.boot/sum.denom.boot
      
      list(sum.counts.boot,prop.boot,sum.denom.boot)
      
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
    
    results.boot.denom <- do.call(cbind, lapply(results.boot.list, `[[`, 3))
    pred.denom <- apply(results.boot.denom,1,median)
    
  }
  
  data_predicted <- data.frame(date=as.Date(unique(data.new$date)),
                               est_raw_counts=pred,
                               ci_raw_counts_low=pi.low,
                               ci_raw_counts_up=pi.up,
                               est_prop=pred.prop,
                               ci_low_prop=pi.low.prop,
                               ci_up_prop=pi.up.prop)
  
  # Calculate observed (fill in missing with predicted values) 
  indicator_results %>% 
    filter(site %in% facility_complete_list) %>%
    mutate(observed = case_when(is.na(observed) ~ est_raw_counts,
                                TRUE ~ as.numeric(observed)))  %>% 
    dplyr::select(site,date,observed) -> indicator_results_nomiss
  
  if (counts_only==TRUE){
    
    indicator_results_nomiss %>% 
      group_by(date) %>%
      dplyr::summarize(observed_count=sum(observed),
                       observed_prop=NA) -> data_observed
    
  } else {
    
    denom_results %>% 
      filter(site %in% facility_complete_list) %>%
      mutate(observed = case_when(is.na(observed) ~ est_raw_counts,
                                  TRUE ~ as.numeric(observed))) %>% 
      dplyr::select(site,date,observed) -> denom_results_nomiss
    
    left_join(indicator_results_nomiss,denom_results_nomiss,by=c("site","date")) %>%
      group_by(date) %>%
      dplyr::summarize(observed_count=sum(observed.x),
                       observed_denom=sum(observed.y),
                       observed_prop=observed_count/observed_denom) -> data_observed
    
  }
  
  
  # final data frame
  left_join(data_predicted,data_observed,by="date") %>%
    mutate(ci_low_prop = ifelse(ci_low_prop < 0, 0, ci_low_prop),
           ci_up_prop = ifelse(ci_up_prop > 1, 0, ci_up_prop),
           ci_raw_counts_low=ifelse(ci_raw_counts_low < 0,0,ci_raw_counts_low),
           ci_raw_counts_up=ifelse(ci_raw_counts_up < 0,0,ci_raw_counts_up)) -> results
  
  results <- results %>% mutate(site = data %>% rename(site = site_var) %>% distinct(site) %>% pull())
  
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
