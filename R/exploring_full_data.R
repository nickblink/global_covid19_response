library(dplyr)
library(lubridate)
library(ggplot2)
library(gtools)

setwd('C:/Users/nickl/Documents/global_covid19_response/')

D = readRDS('data/liberia_cleaned_NL.rds')

##### Counts of basic things #####

length(unique(D$county)) # 15
length(unique(D$district)) # 92
length(unique(D$facility)) # 912

mean(D$indicator_count_ari_total, na.rm = T) # 42
mean(D$indicator_denom, na.rm = T) # 523

mean(is.na(D$indicator_count_ari_total)) # 0.36
mean(is.na(D$indicator_denom)) # 0.27

##### Missingness at different spatial hierarchical groupings ####

facility_miss = D %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

district_miss = D %>% 
  group_by(district) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

county_miss = D %>% 
  group_by(county) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

par(mfrow=c(2,3))

# plot(density(county_miss$ari_miss), main = 'county ARI missing')
# plot(density(district_miss$ari_miss), main = 'district ARI')
# plot(density(facility_miss$ari_miss), main = 'facility ARI')
# 
# plot(density(county_miss$denom_miss), main = 'county denom missing')
# plot(density(district_miss$denom_miss), main = 'district denom')
# plot(density(facility_miss$denom_miss), main = 'facility denom')

hist(county_miss$ari_miss, main = 'county ARI missing', xlab = 'proportion missing')
hist(district_miss$ari_miss, main = 'district ARI')
hist(facility_miss$ari_miss, main = 'facility ARI')

hist(county_miss$denom_miss, main = 'county denom missing', xlab = 'proportion missing')
hist(district_miss$denom_miss, main = 'district denom')
hist(facility_miss$denom_miss, main = 'facility denom')

##### Missingness by year and season #####
D = D %>%
  mutate(month = month(date),
         year = year(date))

date_miss = D %>% 
  group_by(date) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

month_miss = D %>% 
  group_by(month) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

year_miss = D %>% 
  group_by(year) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)))

# par(mfrow = c(2,3))
# plot(date_miss$date, date_miss$ari_miss, xlab = 'date', ylab = 'proportion missing', main = 'ARI missingness by date')
plot_list = list()
plot_list[[1]] <- ggplot(date_miss, aes(x = date, y = ari_miss)) + 
  geom_point() + 
  ggtitle('ARI missingness by date') + 
  xlab('date') +
  ylab('proportion missing')

plot_list[[2]] <- ggplot(month_miss, aes(x = month, y = ari_miss)) + 
  geom_bar(stat = 'identity') + 
  ggtitle('ARI missingness by month')

plot_list[[3]] <- ggplot(year_miss, aes(x = year, y = ari_miss)) + 
  geom_bar(stat = 'identity') + 
  ggtitle('ARI missingness by year')

plot_list[[4]] <- ggplot(date_miss, aes(x = date, y = denom_miss)) + 
  geom_point() + 
  ggtitle('denom missingness by date') + 
  xlab('date') +
  ylab('proportion missing')

plot_list[[5]] <- ggplot(month_miss, aes(x = month, y = denom_miss)) + 
  geom_bar(stat = 'identity') + 
  ggtitle('denom missingness by month')

plot_list[[6]] <- ggplot(year_miss, aes(x = year, y = denom_miss)) + 
  geom_bar(stat = 'identity') + 
  ggtitle('denom missingness by year')


cowplot::plot_grid(plotlist = plot_list, ncol = 3)


##### (not done) Missingness vs. CCA values #####


##### Visualizing the missingness pattern #####
### indicator denom

districts = district_miss %>%
  filter(denom_miss > 0.3) %>%
  # filter(denom_miss < 0.2, denom_miss > .05) %>%
  dplyr::select(district) %>% pull

par(mfrow = c(3,3))
for(d in sample(districts, 3)){
  df = D %>% filter(district == d) %>%
    dplyr::select(date, district, facility, indicator_denom)
  
  df_spread = tidyr::spread(df, facility, indicator_denom)
  
  tmp2 = df_spread[,-c(1,2)]
  for(col in colnames(tmp2)){
    tmp2[,col] = as.integer(is.na(tmp2[,col]))
  }
  tmp2 = as.matrix(tmp2)
  
  rownames(tmp2) = as.character(df_spread$date)
  heatmap(tmp2, keep.dendro = F, Rowv = NA)
  # 
  # tt = as.matrix(is.na(df_spread[,-c(1,2)]))
  # tt = apply(tt, 2, as.integer)
  # heatmap(tt)
}

### ARI
districts = district_miss %>%
  filter(ari_miss > 0.3) %>%
  dplyr::select(district) %>% pull

par(mfrow = c(3,3))
for(d in sample(districts, 3)){
  df = D %>% filter(district == d) %>%
    dplyr::select(date, district, facility, indicator_count_ari_total)
  
  df_spread = tidyr::spread(df, facility, indicator_count_ari_total)
  
  tmp2 = df_spread[,-c(1,2)]
  for(col in colnames(tmp2)){
    tmp2[,col] = as.integer(is.na(tmp2[,col]))
  }
  tmp2 = as.matrix(tmp2)
  
  rownames(tmp2) = as.character(df_spread$date)
  heatmap(tmp2, keep.dendro = F, Rowv = NA)
  # 
  # tt = as.matrix(is.na(df_spread[,-c(1,2)]))
  # tt = apply(tt, 2, as.integer)
  # heatmap(tt)
}

##### spatial correlation of missingness #####

# look at correlation of missing values within each district? And within county?

### Indicator denom
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D$district)){
  # spread the data for this district
  df = D %>% filter(district == d) %>%
    select(date, district, facility, indicator_denom)
  df_spread = tidyr::spread(df, facility, indicator_denom)
  
  if(sum(is.na(df_spread)) == 0){
    print(sprintf('skipping district %s because nothing is missing', d))
    next
  }
  
  # get the correlation matrix across facilities
  tmp <- cor(is.na(df_spread[,-c(1,2)]))
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_list = c(cor_list, na.omit(as.vector(tmp)))
  
  # if the average correlation is NaN (say if there is only one missing value)
  if(is.nan(mean(tmp, na.rm = T))){
    print(sprintf('skipping district %s because of NaN correlations', d))
    next
  }
  
  # print(mean(tmp, na.rm=  T))
  res_row = data.frame(district = d, p_miss = mean(is.na(df_spread[,-c(1,2)])), cor_miss = mean(tmp, na.rm=  T))
  res = rbind(res, res_row)
}

res_denom = res
cor_list_denom = cor_list

### ARI
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D$district)){
  # spread the data for this district
  df = D %>% filter(district == d) %>%
    select(date, district, facility, indicator_count_ari_total)
  df_spread = tidyr::spread(df, facility, indicator_count_ari_total)
  
  if(sum(is.na(df_spread)) == 0){
    print(sprintf('skipping district %s because nothing is missing', d))
    next
  }
  
  # get the correlation matrix across facilities
  tmp <- cor(is.na(df_spread[,-c(1,2)]))
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_list = c(cor_list, na.omit(as.vector(tmp)))
  
  # if the average correlation is NaN (say if there is only one missing value)
  if(is.nan(mean(tmp, na.rm = T))){
    print(sprintf('skipping district %s because of NaN correlations', d))
    next
  }
  
  # print(mean(tmp, na.rm=  T))
  res_row = data.frame(district = d, p_miss = mean(is.na(df_spread[,-c(1,2)])), cor_miss = mean(tmp, na.rm=  T))
  res = rbind(res, res_row)
}

res_ARI = res
cor_list_ARI = cor_list

# plotting district means of correlations
par(mfrow = c(1,2))
hist(res_ARI$cor_miss, main = 'within-district ARI', xlab = 'mean correlation')

hist(res_denom$cor_miss, main = 'within-district denom', xlab = 'mean correlation')

# plotting facility correlations
hist(cor_list_ARI, main = 'within-district ARI', xlab = 'correlation')
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.081

hist(cor_list_denom, main = 'within-district denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) # 0.093

#
##### Within county (not district) spatial correlation #####
res = NULL
cor_list_denom = c()
cor_list_ARI = c()

# sorry, but I'm shamefully doing for loops on for loops on for loops. It's slow because of that
for(cc in unique(D$county)){
  # spread the data for this district
  df = D %>% filter(county == cc) %>%
    select(date, county, district, facility, indicator_denom, indicator_count_ari_total)
  #df_spread = tidyr::spread(df, facility, indicator_denom)
  
  # the long way - a for loop!
  for(f in unique(df$facility)){
    target = df %>% filter(facility == f)
    
    # skip if the facility has no missingness
    if(sum(is.na(target)) == 0){
      next
    }
    
    d = target$district[1]
    compare_facs = df %>%
      filter(district != d) %>%
      pull(facility) %>% 
      unique()
    for(f2 in compare_facs){
      compare = df %>% filter(facility == f2)
      cor_list_denom = c(cor_list_denom, cor(is.na(target$indicator_denom), is.na(compare$indicator_denom)))
      cor_list_ARI = c(cor_list_ARI, cor(is.na(target$indicator_count_ari_total), is.na(compare$indicator_count_ari_total)))
    }
  }
}

cor_list_denom = na.omit(cor_list_denom)
cor_list_ARI = na.omit(cor_list_ARI)

# plotting facility correlations
hist(cor_list_ARI, main = 'within-county ARI', xlab = 'correlation')
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.02

hist(cor_list_denom, main = 'within-county denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) 
# 0.048

#
##### temporal correlation of missingness #####

### denom
missing_facs = facility_miss %>% 
  filter(denom_miss > 0) %>%
  select(facility) %>% pull()

# cycle through each facility and pull out blocks of missingness
count = c()
for(m in missing_facs){
  tmp <- D %>% filter(facility == m) %>% select(indicator_denom) %>% pull()
  counter = 0
  for(i in 1:length(tmp)){
    if(is.na(tmp[i])){
      counter = counter + 1
    }else{
      count = c(count, counter)
      counter = 0
    }
  }
}

# remove 0s and look at results
count = count[count != 0]
table(count)


### ARI
missing_facs = facility_miss %>% 
  filter(ari_miss > 0) %>%
  select(facility) %>% pull()

# cycle through each facility and pull out blocks of missingness
count = c()
for(m in missing_facs){
  tmp <- D %>% filter(facility == m) %>% select(indicator_count_ari_total) %>% pull()
  counter = 0
  for(i in 1:length(tmp)){
    if(is.na(tmp[i])){
      counter = counter + 1
    }else{
      count = c(count, counter)
      counter = 0
    }
  }
}

# remove 0s and look at results
count = count[count != 0]
table(count)

##### (Not done) Facility Size and Missingness #####

df = D %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            fac_size = mean(indicator_denom, na.rm = T)) %>%
  filter(!is.nan(fac_size))

plot(density(df$fac_size))
plot(density(log(df$fac_size)))
quantile(df$fac_size)

plot(log(df$fac_size), df$ari_miss)

# maybe barplots by quartile?
df$group = quantcut(df$fac_size, q = 10)

# group by quantile
tt = df %>%
  group_by(group) %>%
  summarise(denom_miss_avg = mean(denom_miss),
            ari_miss_avg = mean(ari_miss))

# plot it!
p1 <- ggplot(tt, aes(x = group, y = ari_miss_avg)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x=element_text(angle=-30)) +
  xlab('facility size quantile') + 
  ylab('average proportion missing') +
  ggtitle('ARI missingness')

p2 <- ggplot(tt, aes(x = group, y = denom_miss_avg)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x=element_text(angle=-30)) +
  xlab('facility size quantile') + 
  ylab('average proportion missing') +
  ggtitle('denominator missingness')

cowplot::plot_grid(p1, p2)
