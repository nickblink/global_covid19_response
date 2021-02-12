# testing the code for the facility-level modeling of syndromatic surveillance

# set the working directory
setwd('C:/Users/nickl/Documents/global_covid19_response/')

# load in functions and data
source("R/model_functions.R")
source("R/model_figures.R")
data <- readRDS("data/data_example_singlecounty.rds")

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

tmp = data %>% dplyr::select(date, facility, indicator_denom)
data_missing = tidyr::spread(tmp, facility, indicator_denom)
tmp2 = data_missing[,-1]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = as.matrix(tmp2)
heatmap(tmp2, Rowv = NA)
# ok other than M not much of continuous missingness in columns

# is under5 the age? That seems pretty important doesn't it? That is actually quite nice for looking at COVID. I would actually expect that to be more important than ari total. Or rather, that over5 would be more important than ari total.

# Declare this for all functions
extrapolation_date <- "2020-01-01"

# Run Facility-level Model
example_1_results <- fit.site.specific.denom.pi(data=data,
                                                site_name="Facility K",
                                                extrapolation_date=extrapolation_date,
                                                indicator_var="indicator_count_ari_total",
                                                denom_var="indicator_denom", 
                                                site_var="facility",
                                                date_var="date",
                                                R=500)   # Number of Boostrap resamples


# plot 'em!
plot_site(example_1_results, 'count')
plot_site(example_1_results, 'proportion')

# get all sites
all_sites <- data %>% distinct(facility) %>% pull()

# loop over all syndromic surveillance indicators and facilities
# groups all facilities into one data frame. The list is over the surveillance indicators
lapply(c("indicator_count_ari_total"), function(y){    # can have a list of more indicators than just ARI
  do.call(rbind, lapply(all_sites, function(x)
    fit.site.specific.denom.pi(data=data,
                               site_name=x,
                               extrapolation_date=extrapolation_date,
                               indicator_var=y,
                               denom_var="indicator_denom",   # corresponding denominator indicator needed for proportions
                               site_var="facility",
                               date_var="date",
                               R=500)))
}
) -> facility.list

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


