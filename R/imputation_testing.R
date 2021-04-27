library(cutoffR)
library(rstanarm)
library(CARBayesST)
library(cowplot)
library(dplyr)

# load in functions and data
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source("R/model_functions.R")
source("R/model_figures.R")
source('R/imputation_functions.R')

# county-level data
data2 <- readRDS('data/ari_total_county_revisions.rds') %>%
  filter(date <  as.Date("2020-01-01"))
mean(is.na(data)) # ok no missingness

# data for all facilities
D = readRDS('data/liberia_cleaned_NL.rds')

##### Modeling all counties #####
# Declare this for all functions
extrapolation_date <- "2020-01-01"
D = D %>% filter(date < extrapolation_date)

# check the distribution of district size
district_size = D %>% group_by(district) %>% summarize(n = length(unique(facility))) %>% arrange(n)
table(district_size$n)

# checking facility missingness
facility_miss = D %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            ari_count = sum(indicator_count_ari_total, na.rm = T),
            ari_nonzero = sum(indicator_count_ari_total > 0, na.rm = T))

# only keep facilities with < 80% missing
# warning('removing Kingdom Care Medical Clinic because it has some highly unusual values')
# facility_list = facility_miss %>% filter(ari_miss < 0.8, ari_count > 0, facility != 'Kingdom Care Medical Clinic') %>% arrange(ari_miss) %>% distinct(facility) %>% pull()
facility_list = facility_miss %>% filter(ari_nonzero >= 10) %>% arrange(ari_miss) %>% distinct(facility) %>% pull()

D = D %>%
  filter(facility %in% facility_list)

### My imputations
# (0) Trying CARBayes on errything
# system.time({
#   tmp <- CARBayes_imputation(D, col = "indicator_count_ari_total", return_type = 'all', AR = 1)
# })
# 15 minutes (if not doing separate betas)
# D = tmp[[1]]
# 
# tt = tmp[[2]]$samples$beta

# save(D, file = 'results/allcounties_samebeta_full_CAR_imputation_04212021.RData')

# don't save it all! So big!
# save(tmp, file = 'results/allcounties_full_CAR_imputation_04212021.RData')

# trying it by each county. How it probably should be done anyway.
if(F){
system.time({
  tmp2 <- lapply(unique(D$county), function(xx){
    print(xx)
    D2 = D %>% filter(county == xx)
    print(dim(D2))
    tmp_lst <- CARBayes_imputation(D2, col = "indicator_count_ari_total", return_type = 'all', AR = 1)
    return(tmp_lst)
  })
}) # 3.5 hours, almost all by Monteserado 
# save(tmp2, file = 'C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/allcounties_CAR_imputation_04222021.RData')
}

### (1) Impute all the values. Now Bayesian imputation! <- version 1
# D <- bayes_periodic_imputation(D, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility')

system.time({
  tmp <- lapply(unique(D$county), function(xx){
    print(xx)
    D2 = D %>% filter(county == xx)
    tmp_df <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility', harmonic_priors = F)
    return(tmp_df)
  })
}) # ~ 6 minutes with 500 iterations (need to do more eventually)

D = do.call('rbind',tmp)

### (2) Impute all the values. Now Bayesian imputation! <- version 2
system.time({
  tmp2 <- lapply(unique(D$county), function(xx){
    print(xx)
    D2 = D %>% filter(county == xx)
    tmp_df <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility', harmonic_priors = T)
    return(tmp_df)
  })
}) # ~ 6 minutes with 500 iterations

D = do.call('rbind',tmp2)

### (3) Impute all values with periodic imputation
D <- periodic_imputation(D, col = "indicator_count_ari_total", family = 'NB', group = 'facility')

# merge in the CAR data
load('C:/Users/nickl/Dropbox/Nick_Cranston/HSPH/Research/Hedt_Synd_Surveillance_Project/allcounties_CAR_imputation_04222021.RData')

tmp <- lapply(tmp2, function(xx){
  xx[[1]]
})
D2 = do.call('rbind',tmp) %>%
  select(date, facility, indicator_count_ari_total_CARBayes_ST_0.5)

D = merge(D, D2, by = c('date','facility')) %>%
  select(date, facility, district, county, indicator_count_ari_total, 
         CARBayes_ST = indicator_count_ari_total_CARBayes_ST_0.5,
         bayes_harmonic_prior_miss = indicator_count_ari_total_pred_bayes_harmonic_miss, 
         bayes_harmonic = indicator_count_ari_total_pred_bayes_harmonic,
         pred_harmonic = indicator_count_ari_total_pred_harmonic)

# rename the columns

D = as.data.frame(D)
# save(D, file = 'results/periodic_bayes_CAR_allcounties_04252021.RData')


##### Plotting all counties #####
load('results/periodic_bayes_CAR_allcounties_04252021.RData')

facility_fits = D %>% group_by(facility, county) %>%
  summarize(CAR = sum(CARBayes_ST),
            harm = sum(pred_harmonic),
            bayes = sum(bayes_harmonic_prior_miss)) %>%
  arrange(desc(CAR))

max(facility_fits$CAR)
max(facility_fits$harm)
max(facility_fits$bayes)
# all good. No crazy values

# par(mfrow = c(3,3))
# for(i in 1:9){
#   f = facility_fits$facility[i]
#   tt = D %>% filter(facility == f)
#   plot(tt$date, tt$indicator_count_ari_total)
# }


#warning('excluding 9 facilities because of the Bayes method, but can revisit later')
facility_list = D %>% group_by(facility) %>%
  summarize(ari_miss = mean(is.na(indicator_count_ari_total)),
            CAR_miss = mean(is.na(CARBayes_ST)),
            Bayes_miss = mean(is.na(bayes_harmonic_prior_miss)),
            harmonic_miss = mean(is.na(pred_harmonic))) %>%
  filter(CAR_miss + Bayes_miss + harmonic_miss == 0) %>%
  ungroup() %>%
  distinct(facility) %>% pull

# df of counties ordered by missingness for plotting
county_missing = D %>% 
  group_by(county) %>%
  summarize(ari_miss = mean(is.na(indicator_count_ari_total))) %>%
  arrange(ari_miss)

# # rename the columns of D
# select(date, facility, district, county, indicator_count_ari_total, 
#        CARBayes_ST = indicator_count_ari_total_CARBayes_ST_0.5,
#        bayes_harmonic_prior_miss = indicator_count_ari_total_pred_bayes_harmonic_miss, 
#        bayes_harmonic = indicator_count_ari_total_pred_bayes_harmonic,
#        pred_harmonic = indicator_count_ari_total_pred_harmonic)


plot_list = lapply(county_missing$county, function(xx){
  D2 = D %>% filter(county == xx,
                    facility %in% facility_list)
  
  missingness = county_missing$ari_miss[county_missing$county == xx]
  
  # call plot_county
  p1 = plot_county_imputations(D2, 
                               imp_vec = c('CARBayes_ST','bayes_harmonic_prior_miss', 'pred_harmonic'), 
                               color_vec= c('blue','red','forestgreen'), title = sprintf('%s (%0.3f)', xx, missingness)) +
    theme(legend.position = 'none')
  
  p2 = p1 + theme(legend.position = "bottom", legend.text=element_text(size=8))
  ggsave(p2, file = sprintf('figures/individual_counties/%s_imputations_04252021.png',xx))
  
  return(p1)
  
})

legend <- get_legend(
  # create some space to the left of the legend
  plot_county_imputations(D %>% filter(county == 'Bomi'), 
                          imp_vec = c('CARBayes_ST','bayes_harmonic_prior_miss', 'pred_harmonic'), color_vec= c('blue','red','forestgreen'), title = 'tmp', labels = c('CARBayes', 'bayes harmonic', 'harmonic')) + theme(legend.position = "bottom") 
)

plot_grid(plot_grid(plotlist = plot_list), legend, ncol = 1, rel_heights = c(1, 0.1))

# ggsave('figures/all_counties_imputation_04252021.png', scale = 3)


### Doing the fitted models rather than just imputed points
plot_list = lapply(county_missing$county, function(xx){
  D2 = D %>% filter(county == xx,
                    facility %in% facility_list)
  
  missingness = county_missing$ari_miss[county_missing$county == xx]
  
  # call plot_county
  p1 = plot_county_fits(D2, imp_vec = c('CARBayes_ST','bayes_harmonic_prior_miss', 'pred_harmonic'),  color_vec= c('blue','red','forestgreen'), title = sprintf('%s (%0.3f)', xx, missingness)) +
    theme(legend.position = 'none')
  
  p2 = p1 + theme(legend.position = "bottom", legend.text=element_text(size=8))
  ggsave(p2, file = sprintf('figures/individual_counties/%s_fits_04252021.png',xx))
  
  return(p1)
  
})

legend <- get_legend(
  # create some space to the left of the legend
  plot_county_fits(D %>% filter(county == 'Bomi'), 
                          imp_vec = c('CARBayes_ST','bayes_harmonic_prior_miss', 'pred_harmonic'), color_vec= c('blue','red','forestgreen'), title = 'tmp', labels = c('CARBayes', 'bayes harmonic', 'harmonic')) + theme(legend.position = "bottom")
)

plot_grid(plot_grid(plotlist = plot_list), legend, ncol = 1, rel_heights = c(1, 0.1))

#ggsave('figures/all_counties_fits_04252021.png', scale = 3)


# plotting county level results (from paper)
if(FALSE){
  # updating the missing values
  D2$ari_bayes_imp2 = D2$indicator_count_ari_total
  D2$ari_bayes_imp2[is.na(D2$ari_bayes_imp2)] = D2$indicator_count_ari_total_pred_bayes_harmonic_miss[is.na(D2$ari_bayes_imp2)]
  
  # (2) Run the same models as before 
  facility.results.bayes <- do.call(rbind, lapply(facility_list,function(x){
    fit.site.specific.denom.pi(data=D2,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var="ari_bayes_imp",
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500,
                               counts_only = T)
  })
  ) 
  plot_facet(facility.results.bayes, type="count")
  
  
  facility.results.bayes2 <- do.call(rbind, lapply(facility_list,function(x){
    fit.site.specific.denom.pi(data=D2,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var="ari_bayes_imp2",
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500)
  })
  ) 
  
  plot_facet(facility.results.bayes2, type="count")
  
  # ok some issues here - to return later
  tmp = D2 %>% filter(facility == 'Gonjeh Clinic')
}


##### Modeling Bomi #####
# Declare this for all functions
extrapolation_date <- "2020-01-01"
D = D %>% filter(date < extrapolation_date)

# filter to only look at one county
# Bomi has 0.255 ARI missingness. Not the highest but reasonably high
D2 = D %>% filter(county == 'Bomi')

# check the distribution of facility missingness
#   any over 0.2? 
facility_miss = D2 %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))
# 13/27 have > 20% missing ARI. And 4 (?) have > 20% missing denom (not sure if the code does anything about that, though)

# only keep facilities with < 80% missing
facility_list = facility_miss %>% filter(ari_miss < 0.8) %>% arrange(ari_miss) %>% distinct(facility) %>% pull()

D2 = D2 %>%
  filter(facility %in% facility_list)

# run models from the paper
if(FALSE){
  
  
  # full model fit (remember this is excluding > 0.2 missing)
  county.fit = fit.aggregate.pi.boot(D2, 
                                     indicator_var = "indicator_count_ari_total",
                                     denom_var = "indicator_denom",
                                     site_var = "facility",
                                     date_var = "date",
                                     counts_only=FALSE,
                                     R=250)
  
  # district fits
  uni_districts = unique(D2$district)
  district.fits <- lapply(uni_districts, function(xx){
    fit.cluster.pi(data = D2 %>% filter(district == xx), 
                   indicator_var = "indicator_count_ari_total",
                   denom_var = "indicator_denom",
                   site_var = "facility",
                   date_var = "date",
                   counts_only=FALSE,
                   R=250) %>%
      mutate(site = xx)
  })
  
  # clinic fits
  facility.results <- do.call(rbind, lapply(facility_list,function(x){
    fit.site.specific.denom.pi(data=D2,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var="indicator_count_ari_total",
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500)
  })
  ) 
  
  # individual facility results
  plot_facet(facility.results,type="count")
  
  # full county results
  plot_site(county.fit, "count", title="Bomi Aggregated Results")
}


### My imputations
# (0) Trying CARBayes!
tmp <- CARBayes_imputation(D2, col = "indicator_count_ari_total", return_type = 'all')
D2 <- tmp[[1]]

# (1) Impute all the values. Now Bayesian imputation! <- version 1
D2 <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility')

# updating the missing values
D2$ari_bayes_imp = D2$indicator_count_ari_total
D2$ari_bayes_imp[is.na(D2$ari_bayes_imp)] = D2$indicator_count_ari_total_pred_bayes_harmonic[is.na(D2$ari_bayes_imp)]

# (2) Impute all the values. Now Bayesian imputation! <- version 2
D2 <- bayes_periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', iterations = 500, group = 'facility', harmonic_priors = T)

# (3) Impute all values with periodic imputation
D2 <- periodic_imputation(D2, col = "indicator_count_ari_total", family = 'NB', group = 'facility')

# plotting county level results
if(FALSE){
  # updating the missing values
  D2$ari_bayes_imp2 = D2$indicator_count_ari_total
  D2$ari_bayes_imp2[is.na(D2$ari_bayes_imp2)] = D2$indicator_count_ari_total_pred_bayes_harmonic_miss[is.na(D2$ari_bayes_imp2)]
  
  # (2) Run the same models as before 
  facility.results.bayes <- do.call(rbind, lapply(facility_list,function(x){
    fit.site.specific.denom.pi(data=D2,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var="ari_bayes_imp",
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500,
                               counts_only = T)
  })
  ) 
  plot_facet(facility.results.bayes, type="count")
  
  
  facility.results.bayes2 <- do.call(rbind, lapply(facility_list,function(x){
    fit.site.specific.denom.pi(data=D2,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var="ari_bayes_imp2",
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500)
  })
  ) 
  
  plot_facet(facility.results.bayes2, type="count")
  
  # ok some issues here - to return later
  tmp = D2 %>% filter(facility == 'Gonjeh Clinic')
}

D2 = as.data.frame(D2)

plot_imputations <- function(df, imp_vec, color_vec, fac_list = NULL){
  # get facility list if not supplied
  if(is.null(fac_list)){
    fac_list = unique(df$facility)
  }
  
  # rename columns of df to make it easier
  for(col in imp_vec){
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
    for(j in 1:length(imp_vec)){
      col = imp_vec[j]
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
      
      tmp2 = tmp %>% dplyr::select(date, y = imp_vec[j]) %>% mutate(method = col)
      df_f = rbind(df_f, tmp2)
      
      #scale_colour_manual(name="Error Bars",values=cols)
      #p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = imp_vec[j]), color = color_vec[j]))
     # p1 = suppressWarnings(p1 + geom_line(data = tmp, aes_string(x = 'date', y = imp_vec[j], color = color_vec[j])))
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

bomi_plot <- plot_imputations(D2, 
                 imp_vec = c('CARBayes_ST','pred_bayes_harmonic_miss', 'pred_harmonic'), 
                 color_vec= c('blue','red','orange'), 
                 facility_list)
bomi_plot

#ggsave('figures/Bomi_imputation_04072021.png',scale = 3)
#cowplot::save_plot(filename = 'figures/Bomi_imputation_04072021_v2.png',plot = bomi_plot)

plot_county_imputations <- function(df, imp_vec, color_vec){
  
  # rename columns of df to make it easier
  for(col in imp_vec){
    ind = grep(col, colnames(df))
    if(length(ind) != 1){browser()}
    colnames(df)[ind] = col
  }
  
  # initialize the data frame to store final results
  df_f = NULL
  
  for(j in 1:length(imp_vec)){
    col = imp_vec[j]
    
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
  
  # plot them all!
  p1 <- ggplot(data = df_f, aes(x = date, y = y, group = method, color = method)) + 
    geom_line() +
    ggtitle('Aggregated Bomi County Predictions') + 
    scale_color_manual(values = c(color_vec)) + theme(legend.position = 'bottom', legend.text=element_text(size=20))
  
  return(p1)
}

plot_county_imputations(D2, 
                        imp_vec = c('CARBayes_ST','pred_bayes_harmonic_miss', 'pred_harmonic'), 
                        color_vec= c('blue','red','orange'))

# ggsave('figures/Bomi_County_plot_04112021.png',scale = 3)

# note the original bayesian method will probably give really similar betas but with smaller confidence intervals, which is sort of a lie.




##### Comparing my OG imputation to github code #####
tmp <- data %>% filter(facility == 'Facility D')
tt <- periodic_imputation(tmp, 'indicator_count_ari_total') %>%
  dplyr::select(date, indicator_count_ari_total, indicator_count_ari_total_pred_harmonic)
ss <- fit.site.specific.denom.pi(tmp, 'Facility D',extrapolation_date = as.Date("2020-01-01"), 'indicator_count_ari_total', site_var = 'facility', date_var = 'date', counts_only = T) %>%
  dplyr::select(date, est_raw_counts)

test = merge(tt, ss, by = 'date')

# ok good it works. This doesn't account for denominator, but it works.

##### Running imputation (older) #####

# res = impute_wrapper(data2, col = 'observed_count', p_vec = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), cutoff_method = NA, bayes_iterations = 2000, group = 'county')
# plot_results(res)
# 

#save(res, file = 'results/allcounties_bayes_imputation.RData')

data_miss = simulate_imputation(data2, 'observed_count', p = 0.9, group = 'county')

tt = bayes_periodic_imputation(data_miss, df_OG = data2, col = 'observed_count', group = 'county', family = 'NB', period = 12, iterations = 500, harmonic_priors = T)

# 
# res = impute_wrapper(data2, col = 'observed_count', p_vec = c(0.8), cutoff_method = NA, bayes_iterations = 200, group = 'county')
