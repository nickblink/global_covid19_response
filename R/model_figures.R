library(tidyverse) 
library(ggplot2)
library(patchwork)
library(scales)
library(lubridate)
library(prettyGraphs)


plot_site <- function(input,type,fontsize=14,title=""){
  
  xlab="Date"
    
    if(type == "count"){
      ylab="Number of New Cases"
      ggplot(input) + 
        geom_line(aes(date,est_count),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_count_low,ymax=ci_count_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed),color="black") + 
        ggtitle(title)+
        geom_vline(data=input, mapping=aes(xintercept=as.Date(extrapolation_date)),linetype="dashed") + 
        theme_bw() + 
        theme(text = element_text(size=fontsize)) + 
        ylab(ylab) + 
        xlab("Date")
    }
    
    else if(type == "proportion"){
      ylab="Proportion of New Cases"
      ggplot(input) +  
        geom_line(aes(date,est_prop),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_prop_low,ymax=ci_prop_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed_prop),color="black") + 
        ggtitle(title)+
        geom_vline(data=input, mapping=aes(xintercept=as.Date(extrapolation_date)),linetype="dashed") + 
        theme_bw() + 
        theme(text = element_text(size=fontsize)) + 
        ylab(ylab) +
        xlab("Date") }
    
    else{
      print("Please specify the plot type as either count or proportion")
    }

  
}
 
  
  
plot_facet <- function(input,type,fontsize=14,title=""){
  
  input <- input %>% dplyr::rename(facet_var=`site`)
  
  if(type=="count"){
    
    ylab="Number of New Cases"
    
      ggplot(input) + geom_line(aes(date,est_count),color="red") + 
        geom_ribbon(aes(x=date,ymin=ci_count_low,ymax=ci_count_up),fill="red",alpha=.2) + 
        geom_line(aes(date,observed),color="black") + 
        facet_wrap(~as.factor(facet_var),scale="free_y") + 
        geom_vline(data=input, mapping=aes(xintercept=as.numeric(as.Date(extrapolation_date))),linetype="dashed") + 
        ggtitle(title)+
        theme_bw() + 
        theme(text = element_text(size=fontsize)) + 
        ylab(ylab) +
        xlab("Date")
      
    } else if(type=="proportion"){
      
      ylab="Proportion of New Cases"
      
    ggplot(input) + geom_line(aes(date,est_prop),color="red") + 
      geom_ribbon(aes(x=date,ymin=ci_prop_low,ymax=ci_prop_up),fill="red",alpha=.2) + 
      geom_line(aes(date,observed_prop),color="black") + 
      facet_wrap(~as.factor(facet_var),scale="free_y") + 
      geom_vline(data=input, mapping=aes(xintercept=as.numeric(as.Date(extrapolation_date))),linetype="dashed") + 
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



# CREATE HEATMAP 

plot_heatmap <- function(input,type,fontsize=14,title=""){
  
  if(type=="count"){
    input %>% 
      filter(date >= extrapolation_date)  %>% 
      mutate(date_new = paste0("0",month(date),"-2020")) %>% 
      mutate(deviation_final = (observed - est_count)/est_count,
             sig = case_when(
               observed >= ci_count_up ~ 1,
               observed <= ci_count_low ~ 1,
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
            text = element_text(size=fontsize)) + 
    scale_fill_gradient2(low = muted("blue"),
                         mid = "grey90",
                         high = muted("red"),
                         breaks = c(-2.5,0,2.5),
                         limits = c(-2.5,2.5),
                         labels = c("lower","as expected","higher")) + 
    labs(fill="Deviation") +
    ggtitle(title)
      
  }
  
  else if(type == "proportion"){
    input %>% 
      filter(date >= extrapolation_date)  %>% 
      mutate(date_new = paste0("0",month(date),"-2020")) %>% 
      mutate(deviation_final = (observed_prop - est_prop)/est_prop,
             sig = case_when(
               observed_prop >= ci_prop_up ~ 1,
               observed_prop <= ci_prop_low ~ 1,
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
            text = element_text(size=fontsize)) + 
      scale_fill_gradient2(low = muted("blue"),
                           mid = "white",
                           high = muted("red"),
                           labels = c("lower","as expected","higher")) + 
      labs(fill="Deviation") +
      ggtitle(title)
  }
  
  else{
    print("Please specify the plot type as either count or proportion")
  }
  
}


