# Preparing Data 

data=df
site_name="JJ Dossen Hospital"
extrapolation_date=extrapolation_date
indicator_var="indicator_count_ari_total"
denom_var="indicator_denom"
site_var="facility"
date_var="date"
R=R
counts_only = TRUE
period=12


df %>% 
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


data.new %>% filter(!is.na(indicator)) %>% 
  filter(date < extrapolation_date) -> data.base



# Function to pass into tsboot
model <- function(d){
  
  data.base <- d
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  mod_counts <- glm(as.formula(formula_counts), data=data.base, family=quasipoisson)
  return(mod_counts$coefficients)
  
}

#model(data.base)


#Prediction Intervals using glm output
result <- fit.site.specific.denom.pi(data=df,
                                     site_name="JJ Dossen Hospital",
                                     extrapolation_date=extrapolation_date,
                                     indicator_var="indicator_count_ari_total",
                                     denom_var="indicator_denom", 
                                     site_var="facility",
                                     date_var="date",
                                     R=R,
                                     counts_only = TRUE)

naive_pi <- function(data.new){
  formula_counts = "indicator ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3"
  mod_counts <- glm(as.formula(formula_counts), data=data.new, family=quasipoisson)
  
  result_alpha_hat <- mod_counts$coefficients
  result_alpha_vcov <- vcov(mod_counts)
  overdisp <- summary(mod_counts)$dispersion
  result_pred <- matrix(nrow = 48,ncol=500)
  for (s in 1:500){
    result_alpha_boot <- MASS::mvrnorm(1,result_alpha_hat,result_alpha_vcov)
    result_pred_boot <- (data.new %>% 
                           mutate(intercept=1) %>%
                           dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                           as.matrix())%*%as.matrix(result_alpha_boot)
    
    result_pred_boot_exp <- exp(result_pred_boot)
    result_pred_boot_exp[which(is.na(result_pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
    result_se_boot <- result_pred_boot_exp/(overdisp - 1)
    
    result_pred[,s] <- MASS::rnegbin(n = nrow(data.new),mu = result_pred_boot_exp,theta = result_se_boot)
  }
  
  result.pi.low <- apply(result_pred,1,quantile,.025,na.rm=TRUE)
  result.pi.up <- apply(result_pred,1,quantile,.975,na.rm=TRUE)
  result.predict <- apply(result_pred,1,quantile,.50,na.rm=TRUE)
  
  # create data frame with predictions, 95% CIs, and observed values 
  result.results.df <- data.frame(site=site_name,
                                  date=as.Date(data.new$date),
                                  est_raw_counts=result.predict,
                                  ci_raw_counts_low=result.pi.low,
                                  ci_raw_counts_up=result.pi.up,
                                  observed=data.new$indicator,
                                  est_prop=NA,
                                  ci_low_prop=NA,
                                  ci_up_prop=NA,
                                  observed_prop=NA)
  return(result.results.df)
  
}



# Prediction Intervals using tsboot mbb output
mbb_pi <- function(blocksize,data.new){
  mbb <- tsboot(data.new,model,500,sim="fixed",blocksize)
  mbb_alpha_hat <- mbb$t0
  mbb_alpha_vcov <- cov(mbb$t,mbb$t)
  mbb_pred <- matrix(nrow = 48,ncol=500)
  for (s in 1:500){
    mbb_alpha_boot <- MASS::mvrnorm(1,mbb_alpha_hat,mbb_alpha_vcov)
    mbb_pred_boot <- (data.new %>% 
                        mutate(intercept=1) %>%
                        dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                        as.matrix())%*%as.matrix(mbb_alpha_boot)
    
    mbb_pred_boot_exp <- exp(mbb_pred_boot)
    mbb_pred_boot_exp[which(is.na(mbb_pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
    mbb_se_boot <- mbb_pred_boot_exp/(overdisp - 1)
    
    mbb_pred[,s] <- MASS::rnegbin(n = nrow(data.new),mu = mbb_pred_boot_exp,theta = mbb_se_boot)
  }
  
  mbb.pi.low <- apply(mbb_pred,1,quantile,.025,na.rm=TRUE)
  mbb.pi.up <- apply(mbb_pred,1,quantile,.975,na.rm=TRUE)
  mbb.predict <- apply(mbb_pred,1,quantile,.50,na.rm=TRUE)
  
  # create data frame with predictions, 95% CIs, and observed values 
  mbb.results.df <- data.frame(site=site_name,
                               date=as.Date(data.new$date),
                               est_raw_counts=mbb.predict,
                               ci_raw_counts_low=mbb.pi.low,
                               ci_raw_counts_up=mbb.pi.up,
                               observed=data.new$indicator,
                               est_prop=NA,
                               ci_low_prop=NA,
                               ci_up_prop=NA,
                               observed_prop=NA)
  return(mbb.results.df)
}




# Prediction Intervals using tsboot stationary output
sb_pi <- function(blocksize,data.new){
  sb <- tsboot(data.new,model,500,sim="geom",blocksize)
  sb_alpha_hat <- sb$t0
  sb_alpha_vcov <- cov(sb$t,sb$t)
  sb_pred <- matrix(nrow = 48,ncol=500)
  for (s in 1:500){
    sb_alpha_boot <- MASS::mvrnorm(1,sb_alpha_hat,sb_alpha_vcov)
    sb_pred_boot <- (data.new %>% 
                       mutate(intercept=1) %>%
                       dplyr::select(intercept,year,ends_with("1"),ends_with("2"),ends_with("3")) %>%
                       as.matrix())%*%as.matrix(sb_alpha_boot)
    
    sb_pred_boot_exp <- exp(sb_pred_boot)
    sb_pred_boot_exp[which(is.na(sb_pred_boot_exp))] <- 1 # have to do this so the rnegbin runs - will be ignored in final step
    sb_se_boot <- sb_pred_boot_exp/(overdisp - 1)
    
    sb_pred[,s] <- MASS::rnegbin(n = nrow(data.new),mu = sb_pred_boot_exp,theta = sb_se_boot)
  }
  
  sb.pi.low <- apply(sb_pred,1,quantile,.025,na.rm=TRUE)
  sb.pi.up <- apply(sb_pred,1,quantile,.975,na.rm=TRUE)
  sb.predict <- apply(sb_pred,1,quantile,.50,na.rm=TRUE)
  
  # create data frame with predictions, 95% CIs, and observed values 
  sb.results.df <- data.frame(site=site_name,
                              date=as.Date(data.new$date),
                              est_raw_counts=sb.predict,
                              ci_raw_counts_low=sb.pi.low,
                              ci_raw_counts_up=sb.pi.up,
                              observed=data.new$indicator,
                              est_prop=NA,
                              ci_low_prop=NA,
                              ci_up_prop=NA,
                              observed_prop=NA)
  return(sb.results.df)
}




plot_site <- function(input,site_name=0,ylab="Number of New Cases",xlab="Date",text_size=14){
  
  if(site_name!=0){ input %>% filter(site==site_name) -> input }
  
  ggplot(input) + 
    geom_line(aes(date,est_raw_counts),color="red") + 
    geom_ribbon(aes(x=date,ymin=ci_raw_counts_low,ymax=ci_raw_counts_up),fill="red",alpha=.2) + 
    geom_line(aes(date,observed),color="black") + 
    ggtitle(paste0(`site_name`))+
    geom_vline(data=input, mapping=aes(xintercept=as.Date(extrapolation_date)),linetype="dashed") + 
    pretty_plot() + theme(text = element_text(size=text_size)) + 
    ylab(ylab) + xlab(xlab)
  
}
plot_site(result,"JJ Dossen Hospital")
plot_site(result.results.df,"JJ Dossen Hospital")
plot_site(mbb.results.df,"JJ Dossen Hospital")
plot_site(sb.results.df,"JJ Dossen Hospital")



#true facility counts model
start_date <- "2016-01-01"
end_date <- "2019-12-01"
period <- 12
X <- data.frame(date = seq.Date(from=as.Date(start_date),to=as.Date(end_date),by="month"))
X <- X %>% dplyr::mutate(total = 1:n(),Intercept = 1,year = year(date) - min(year(date)) + 1,
                         cos1 = cos(2*1*pi*total/period),
                         sin1 = sin(2*1*pi*total/period),
                         cos2 = cos(2*2*pi*total/period),
                         sin2 = sin(2*2*pi*total/period),
                         cos3 = cos(2*3*pi*total/period),
                         sin3 = sin(2*3*pi*total/period))



coef_1_counts <- c(5.772575205, -0.076948755,  0.043700424, -0.005004144,  0.059776802, -0.071432970, -0.063551902, -0.054741817)
true_X_1_counts <- data.frame(exp(as.matrix(X[,3:10]) %*% as.matrix(coef_1_counts))) 
names(true_X_1_counts) <- "observed"
true_counts <- data.frame(true_X_1_counts)
names(true_counts) <- "observed"
true_counts_facility <- true_X_1_counts

bootstrap_success_counts_naive <- matrix(nrow=48,ncol=1000)
bootstrap_success_counts_mbb <- matrix(nrow=48,ncol=1000)
bootstrap_success_counts_sb <- matrix(nrow=48,ncol=1000)


for (s in 1:100){
  # X_1_counts <- rpois(48,true_X_1_counts$observed)
  X_1_counts <- rnegbin(48,true_X_1_counts$observed,true_X_1_counts$observed/(14.39184 - 1))
  
  test_df <- data.frame(date=rep(X$date),county = c(rep("A",48)),district = c(rep("AA",48)),
                        facility = c(rep("AAA",48)))
  test_df <- cbind(test_df,data.frame(indicator = X_1_counts))

  test_df %>% 
    dplyr::select(date,indicator) %>%
    arrange(date) %>% 
    dplyr::mutate(year = year(date) - min(year(date)) + 1,
                  total=1:n(),
                  cos1 = cos(2*1*pi*total/period),
                  sin1 = sin(2*1*pi*total/period),
                  cos2 = cos(2*2*pi*total/period),
                  sin2 = sin(2*2*pi*total/period),
                  cos3 = cos(2*3*pi*total/period),
                  sin3 = sin(2*3*pi*total/period)) -> test_df.new
  all_sites_test <- test_df %>% distinct(facility) %>% pull()
  
  # sim_naive_pi <- naive_pi(test_df.new)
  sim_mbb_pi <- mbb_pi(12,test_df.new)
  sim_sb_pi <- sb_pi(12,test_df.new)
  
  # coverage_counts_naive <- (sim_naive_pi %>% filter(ci_raw_counts_low <= observed & observed <= ci_raw_counts_up) %>% nrow())/48
  coverage_counts_mbb <- (sim_mbb_pi %>% filter(ci_raw_counts_low <= observed & observed <= ci_raw_counts_up) %>% nrow())/48
  coverage_counts_sb <- (sim_sb_pi %>% filter(ci_raw_counts_low <= observed & observed <= ci_raw_counts_up) %>% nrow())/48
  
  
  # bootstrap_success_counts_naive[s] <- coverage_counts_naive
  bootstrap_success_counts_mbb[s] <- coverage_counts_mbb
  bootstrap_success_counts_sb[s] <- coverage_counts_sb
  
  
  print(paste0("Trial ",s))
  # print(mean(bootstrap_success_counts_naive,na.rm=TRUE))
  print(mean(bootstrap_success_counts_mbb,na.rm=TRUE))
  print(mean(bootstrap_success_counts_sb,na.rm=TRUE))
  
}

