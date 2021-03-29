library(ggfortify)
library(forecast)
library(ggplot2)
library(lmtest)

source("R/model_functions.R")


# RESIDUALS PLOT
plot_residuals <- function(input,type,extrapolation_date,fontsize=14,title=""){
  
  xlab="Date"
  ylab="Residuals"
  input %>% filter(date < as.Date(extrapolation_date)) -> input2
  
  if(type == "count"){
    
    input2 %>%
      mutate(residuals_count = observed - est_count) %>%
      ggplot(.) + 
      geom_point(aes(date,residuals_count)) +
      geom_hline(data=input2, mapping=aes(yintercept=0)) + 
      ggtitle(title)+
      theme_bw() + 
      theme(text = element_text(size=fontsize)) + 
      ylab(ylab) + 
      xlab("Date")
  }
  
  else if(type == "proportion"){
    
    input2 %>%
      mutate(residuals_prop = observed_prop - est_prop) %>%
      ggplot(.) + 
      geom_point(aes(date,residuals_prop)) +
      geom_hline(data=input2, mapping=aes(yintercept=0)) + 
      ggtitle(title)+
      theme_bw() + 
      theme(text = element_text(size=fontsize)) + 
      ylab(ylab) + 
      xlab("Date")
    }
  
  else{
    print("Please specify the plot type as either count or proportion")
  }
}

# PACF PLOTS
plot_pacf <- function(input,type,extrapolation_date,fontsize=14,title=""){
  
  input %>% filter(date < as.Date(extrapolation_date)) -> input2
  input2 %>% distinct(date) %>% nrow() -> lag.max
  
  if(type == "count"){
    input2 %>%
      mutate(residuals_count = observed - est_count) %>% 
      pull(residuals_count) -> residuals
  }
  else if(type == "proportion"){
    input2 %>%
      mutate(residuals_prop = observed_prop - est_prop) %>%
      pull(residuals_prop) -> residuals
  }  else{
    print("Please specify the plot type as either count or proportion")
  }
  
  autoplot(Pacf(residuals, lag.max=lag.max, plot=FALSE)) + 
    ylim(c(-1,1)) + theme_bw() + theme(text = element_text(size=fontsize)) + ggtitle(title)
  
}

# ACF PLOTS
plot_acf <- function(input,type,extrapolation_date,fontsize=14,title=""){
  
  input %>% filter(date < as.Date(extrapolation_date)) -> input2
  input2 %>% distinct(date) %>% nrow() -> lag.max
  
  if(type == "count"){
    input2 %>%
      mutate(residuals_count = observed - est_count) %>% 
      pull(residuals_count) -> residuals
  }
  else if(type == "proportion"){
    input2 %>%
      mutate(residuals_prop = observed_prop - est_prop) %>%
      pull(residuals_prop) -> residuals
  }  else{
    print("Please specify the plot type as either count or proportion")
  }
  
  autoplot(Acf(residuals, lag.max=lag.max, plot=FALSE)) + 
    ylim(c(-1,1)) + theme_bw() + theme(text = element_text(size=fontsize)) + ggtitle(title)
  
}


# BREUSCH GODFREY TEST
bgtest_multiple_pval <- function(data,
                                 site_var,
                                 extrapolation_date,
                                 indicator_var,
                                 denom_var=NA, # if proportions are desired
                                 date_var,
                                 type, #count or proportion
                                 order=12){
  
  data %>% dplyr::rename(site=site_var,indicator=indicator_var,date=date_var) -> data.new
  
  if(type == "count"){
    sapply(all_sites,function(x) {
      lmtest::bgtest(site.model.coefficients(data=data.new,
                                             site_name=x,
                                             extrapolation_date,
                                             indicator_var="indicator",
                                             denom_var=denom_var,
                                             site_var="site",
                                             date_var="date",
                                             R=1)$fitted.count,
                     order=order)[4]
    }) -> bgtest.result
  } else if(type == "prop"){
    sapply(all_sites,function(x) {
      lmtest::bgtest(site.model.coefficients(data=data.new,
                                             site_name=x,
                                             extrapolation_date,
                                             indicator_var="indicator",
                                             denom_var=denom_var,
                                             site_var="site",
                                             date_var="date",
                                             R=1)$fitted.prop,
                     order=order)[4]
    }) -> bgtest.result
    
  }
  do.call(rbind,bgtest.result)
}
