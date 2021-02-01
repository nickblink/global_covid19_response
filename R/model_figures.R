

library(tidyverse) 
library(BuenColors)
library(ggplot2)
library(patchwork)
library(scales)
library(lubridate)
library(prettyGraphs)



plot_site <- function(input,type){
  
  xlab="Date"
  text_size=14
  site_name <- input %>% distinct(site) %>% pull()
  if(site_name!=0){ input %>% filter(site==site_name) -> input }
  
    
    if(type == "count"){
      ylab="Number of New Cases"
      ggplot(input) + 
        geom_line(aes(date,est_raw_counts),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_raw_counts_low,ymax=ci_raw_counts_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed),color="black") + 
        ggtitle(paste0(`site_name`))+
        geom_vline(data=input, mapping=aes(xintercept=as.Date(extrapolation_date)),linetype="dashed") + 
        pretty_plot() + theme(text = element_text(size=text_size)) + 
        ylab(ylab) + xlab(xlab)
    }
    
    else if(type == "proportion"){
      ylab="Proportion of New Cases"
      ggplot(input) +  
        geom_line(aes(date,est_prop),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_low_prop,ymax=ci_up_prop),fill="red",alpha=.2) + 
        geom_line(aes(date,observed_prop),color="black") + 
        ggtitle(paste0(`site_name`))+
        geom_vline(data=input, mapping=aes(xintercept=as.Date(extrapolation_date)),linetype="dashed") + 
        pretty_plot() + theme(text = element_text(size=text_size)) + 
        ylab(ylab) + xlab(xlab) }
    
    else{
      print("Please specify the plot type as either count or proportion")
    }

  
}
 
  
  


plot_facet <- function(input,type){
  
  counts_var = "est_raw_counts"
  ylab="Number of New Cases"
  xlab="Date"
  text_size=10
  input <- input %>% dplyr::rename(facet_var=`site`)
  
  if(type=="count"){
    if (grepl("raw",counts_var)){
      ggplot(input) + geom_line(aes(date,est_raw_counts),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_raw_counts_low,ymax=ci_raw_counts_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed),color="black") + 
        facet_wrap(~as.factor(facet_var),scale="free_y") + 
        geom_vline(data=input, mapping=aes(xintercept=as.numeric(as.Date(extrapolation_date))),linetype="dashed") + 
        pretty_plot() + theme(text = element_text(size=text_size)) + 
        ylab(ylab) + xlab(xlab)
    }
    else {
      ggplot(input) + geom_line(aes(date,est_adj_counts),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_adj_counts_low,ymax=ci_adj_counts_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed_prop),color="black") + 
        facet_wrap(~as.factor(facet_var),scale="free_y") + 
        geom_vline(data=input, mapping=aes(xintercept=as.numeric(as.Date(extrapolation_date))),linetype="dashed") + 
        pretty_plot() + theme(text = element_text(size=text_size)) + 
        ylab(ylab) + xlab(xlab)
    }
  }
  
  else if(type=="proportion"){
    ggplot(input) + geom_line(aes(date,est_prop),color="red") + 
      geom_ribbon(aes(x=date,ymin=ci_low_prop,ymax=ci_up_prop),fill="red",alpha=.2) + 
      geom_line(aes(date,observed_prop),color="black") + 
      facet_wrap(~as.factor(facet_var),scale="free_y") + 
      geom_vline(data=input, mapping=aes(xintercept=as.numeric(as.Date(extrapolation_date))),linetype="dashed") + 
      pretty_plot() + theme(text = element_text(size=text_size)) + 
      ylab(ylab) + xlab(xlab)
  }
  
  else{
    print("Please specify the plot type as either count or proportion")
  }
  
}



# CREATE HEATMAP (similar to NYT graphs)

plot_heatmap <- function(input,type){
  
  if(type=="count"){
    input %>% 
      filter(date >= extrapolation_date)  %>% 
      mutate(date_new = paste0("0",month(date),"-2020")) %>% 
      mutate(deviation_final = (observed - est_raw_counts)/est_raw_counts,
             sig = case_when(
               observed >= ci_raw_counts_up ~ 1,
               observed <= ci_raw_counts_low ~ 1,
               TRUE ~ 0)) %>% 
      ggplot(., aes(date_new,site)) + 
      geom_tile(aes(height = 0.3, width=.95, fill = deviation_final, color=as.factor(sig)),size=.3) +
      ylab("") +
      scale_color_manual(values = c("white","black"),guide=FALSE) +
      scale_x_discrete(position = "top") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            text = element_text(size=8)) + 
      scale_fill_gradientn(colors = jdb_palette("ocean_brick"),
                           labels = c("lower","as expected","higher"),
                           breaks=c(-2.5,0,2.5),
                           limits=c(-2.5,2.5)) + 
      labs(fill="Deviation") 
  }
  
  else if(type == "proportion"){
    input %>% 
      filter(date >= extrapolation_date)  %>% 
      mutate(date_new = paste0("0",month(date),"-2020")) %>% 
      mutate(deviation_final = (observed_prop - est_prop)/est_prop,
             sig = case_when(
               observed_prop >= ci_up_prop ~ 1,
               observed_prop <= ci_low_prop ~ 1,
               TRUE ~ 0)) %>% 
      ggplot(., aes(date_new,site)) + 
      geom_tile(aes(height = 0.3, width=.95, fill = deviation_final, color=as.factor(sig)),size=.3) +
      ylab("") +
      scale_color_manual(values = c("white","black"),guide=FALSE) +
      scale_x_discrete(position = "top") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            text = element_text(size=8)) + 
      scale_fill_gradientn(colors = jdb_palette("ocean_brick"),
                           labels = c("lower","as expected","higher"),
                           breaks=c(-2.5,0,2.5),
                           limits=c(-2.5,2.5)) + 
      labs(fill="Deviation") 
  }
  
  else{
    print("Please specify the plot type as either count or proportion")
  }
  
}
