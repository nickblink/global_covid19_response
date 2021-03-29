# testing the code for the facility-level modeling of syndromatic surveillance

# set the working directory
setwd('C:/Users/nickl/Documents/global_covid19_response/')

# load in functions and data
source("R/model_functions.R")
source('R/model_diagnostics.R')
source("R/model_figures.R")
# data <- readRDS("data/data_example_singlecounty.rds")
data <- readRDS("data/data_example.rds")

# exploring data
str(data)
table(data$county)
table(data$district)
table(data$facility)
# ok so 6 districts and 25 facilities. Now, do we aggregate at district or county-level? That is the question. Also, remember that district might be related facility, as in facilities in the same district migh be related

for(col in 5:10){
  print(sum(is.na(data[,col])))
  print(mean(is.na(data[,col])))
}
# ok so not all the same.

##### Correlation plotting #####



data$facility2 = sapply(1:nrow(data), function(xx) sprintf('%s (%s)', data$facility[xx], data$district[xx]))

# missingness
tmp = data %>% dplyr::select(date, facility, indicator_denom)
data_spread = tidyr::spread(tmp, facility, indicator_denom)
tmp2 = data_spread[,-1]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = as.matrix(tmp2)
heatmap(tmp2, Rowv = NA)
# ok other than M not much of continuous missingness in columns


# looking at correlation of data
corrplot::corrplot(cor(data_spread[,-1], use = 'complete.obs'), order = 'hclust', type = 'upper')
# so there is definitely enough correlation to do some analysis. At least for some/most of the sites. Some of the other sites are kinda screwed.

# make unique district plots
par(mfrow = c(2,3))
for(dist in unique(data$district)){
  tmp = data %>% filter(district == dist) %>% dplyr::select(date, facility, indicator_denom)
  tmp = tidyr::spread(tmp, facility, indicator_denom)
  corrplot::corrplot(cor(tmp[,-1], use = 'complete.obs'), order = 'hclust', type = 'upper', main = dist)
}

# doing this the long way
fac_dist = data %>% dplyr::select(facility, district) %>% unique()
within = c()
without = c()
for(i in 1:(nrow(fac_dist)-1)){
  fac = fac_dist$facility[i]
  for(j in (i+1):nrow(fac_dist)){
    fac2 = fac_dist$facility[j]
    if(fac_dist$district[i] == fac_dist$district[j]){
      within = c(within, cor(data_spread[,fac], data_spread[,fac2], use = 'complete.obs'))
    }else{
      without = c(without, cor(data_spread[,fac], data_spread[,fac2], use = 'complete.obs'))
    }
  }
}
par(mfrow = c(1,2))
hist(within, main = 'within-district correlations')
hist(without, main = 'between-district correlations')


#
##### github examples #####

# is under5 the age? That seems pretty important doesn't it? That is actually quite nice for looking at COVID. I would actually expect that to be more important than ari total. Or rather, that over5 would be more important than ari total.

# Declare this for all functions
extrapolation_date <- "2020-01-01"

# Run Facility-level Model
example_1_results <- fit.site.specific.denom.pi(data=data,
                                                site_name="Facility B",
                                                extrapolation_date=extrapolation_date,
                                                indicator_var="indicator_count_ari_total",
                                                denom_var="indicator_denom", 
                                                site_var="facility",
                                                date_var="date",
                                                R=500)   # Number of Boostrap resamples

# Run Facility-level Model
fit.site.specific.denom.pi(data=data,
                           site_var="facility",
                           site_name="Facility K",
                           extrapolation_date=extrapolation_date,
                           indicator_var="ari_cases",
                           denom_var="total_visits", 
                           date_var="date",
                           R=500) -> single.facility.results

# plot 'em!
plot_site(single.facility.results,type="count", title="Acute Respiratory Infections at Facility K")
plot_site(single.facility.results,type="proportion", title="Acute Respiratory Infections at Facility K")


# get all sites
all_sites <- data %>% distinct(facility) %>% pull()

# loop over all syndromic surveillance indicators and facilities
# groups all facilities into one data frame. The list is over the surveillance indicators
# loop over all  facilities
do.call(rbind, lapply(all_sites,function(x){
  fit.site.specific.denom.pi(data=data,
                             site_name=x,
                             extrapolation_date=extrapolation_date,
                             indicator_var="ari_cases",
                             denom_var="total_visits",   # corresponding denominator indicator needed for proportions
                             site_var="facility",
                             date_var="date",
                             R=500)
})
) -> facility.results

head(facility.results)
plot_facet(facility.results,type="count")
plot_heatmap(facility.results,type="count")

data %>% filter(is.na(ari_cases) | is.na(total_visits))

## # A tibble: 3 x 5
##   date       facility   county       ari_cases total_visits
##   <date>     <chr>      <chr>            <dbl>        <dbl>
## 1 2018-08-01 Facility K County Alpha       106           NA
## 2 2019-07-01 Facility Q County Alpha        NA           NA
## 3 2019-10-01 Facility K County Alpha       125           NA

fit.aggregate.pi.boot(data,
                      indicator_var = "ari_cases",
                      denom_var = "total_visits",
                      date_var = "date",
                      site_var = "facility",
                      R=500) -> aggregate.results

head(aggregate.results)

plot_site(aggregate.results, "count", title="Facility K and Q Aggregated Results")

plot_residuals(aggregate.results,
               type="count",
               extrapolation_date="2020-01-01",
               title="Residuals from Aggregated model")

# lazy way - Nick's code
tt = facility.list[[1]]
for(ss in all_sites){
  df = tt %>% filter(site == ss)
  print(ss)
  print(sum(is.na(df)))
}

df = tt %>% filter(site == 'Facility M')

# label each dataframe in the output with the respective indicator names; note this example has only 1
names(facility.list) <- c("indicator_count_ari_total")  

# data frame of all the results
head(facility.list[["indicator_count_ari_total"]])
dim(facility.list[["indicator_count_ari_total"]])

# loops over all denominator indicators, like healthcare utilization, and runs the regression model on these. Creates estimate of raw counts as well as CI on these. Does not estimate proportion because that doesn't make sense here. Though then the question is, how are these CIs used, since before they must have just used the observed counts as is, no? Ah, I think they only fill in missing denominators for the county-level models

# can have a list of more utilization indicators
lapply(c("indicator_denom"), function(y){     
  
  do.call(rbind, lapply(all_sites,function(x)
    fit.site.specific.denom.pi(data=data,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var=y,
                               site_var="facility",
                               date_var="date",
                               counts_only=TRUE)))
  
}) -> facility.list.denom


names(facility.list.denom) <- c("indicator_denom")

head(facility.list.denom[["indicator_denom"]])

# plot the denominator (utilization) and the counts (could also do proportion)
plot_facet(facility.list.denom[[1]], "count")
plot_facet(facility.list[["indicator_count_ari_total"]], "count")

# plot the heat map
plot_heatmap(facility.list[["indicator_count_ari_total"]],"count")

### County level stuff
# getting county results, excluding >20% missing baseline and >50% missing evaluation period.
county_results <- fit.cluster.pi(data = data,
                                 indicator_var = "indicator_count_ari_total",
                                 denom_var = "indicator_denom",
                                 site_var = "facility",
                                 date_var = "date",
                                 counts_only=FALSE,
                                 n_count_base = 0,
                                 p_miss_base = 0.2,
                                 p_miss_eval = 0.5,
                                 R=250)   # Number of Bootstrap resamples)

# look at results - only one county here
head(county_results)

# plot the results
plot_site(county_results, "count")
plot_site(county_results,"proportion")
# there's a spike! Heyo

# alright that's it. Now I should just look into the functions to see what they do.



##### New github examples #####

