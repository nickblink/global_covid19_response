library(dplyr)
library(lubridate)
library(ggplot2)
library(gtools)

#setwd('github_projects/global_covid19_response/')
source('R/imputation_functions.R')
D = readRDS('data/liberia_cleaned_NL.rds')
D2 = readRDS('data/liberia_cleaned_01-06-2021.rds')

remove_outliers <- function(x, k = 5){
  x <- na.omit(x)
  ind <- which(x < median(x) - k*sd(x) | x > median(x) + k*sd(x))
  while(length(ind) > 0){
    #print(x[ind])
    x <- x[-ind]
    ind <- which(x < median(x) - k*sd(x) | x > median(x) + k*sd(x))
  }
  return(x)
}

##### Exploring Maryland County Data for paper (05/15/2025) #####
D <- readRDS('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/liberia_cleaned_01-06-2021.rds')

D2 <- D %>%
  filter(county == 'Maryland') %>%
  select(date, facility, district, ari = indicator_count_ari_total) %>%
  mutate(na = as.integer(is.na(ari))) %>%
  group_by(facility) %>%
  mutate(
    ari_scaled = (ari - min(ari, na.rm = TRUE)) / 
      (max(ari, na.rm = TRUE) - min(ari, na.rm = TRUE)),
  )
D2$na <- as.integer(is.na(D2$ari))

D2$facility_clean <- D2$facility |>
  gsub("Clinic", "", x = _) |>
  gsub("Medical Center", "", x = _) |>
  gsub("Health Center", "", x = _) |>
  gsub("\\s*\\(.*?\\)", "", x = _) |>
  trimws()

# ari un-scaled
library(dplyr)
library(ggplot2)
library(lubridate)

# Define date range
start_date <- as.Date("2016-01-01")
end_date <- as.Date("2019-12-31")
expected <- 48  # expected number of months

# Step 1–3: Count observed and assign generic names
facility_labels <- D2 %>%
  filter(date >= start_date & date <= end_date) %>%
  group_by(facility_clean) %>%
  summarise(observed = sum(!is.na(ari)), .groups = "drop") %>%
  arrange(facility_clean) %>%  # or sort by observed if you prefer
  mutate(facility_number = paste0("Facility ", row_number()),
         label = paste0(facility_number, " (", observed, "/", expected, ")"))

# Ensure the label is an ordered factor
facility_labels$label <- factor(facility_labels$label, levels = facility_labels$label)

# Join labels back to the main data
D2 <- D2 %>%
  left_join(facility_labels %>% select(facility_clean, label), by = "facility_clean")

# Plot
ggplot(data = D2) + 
  geom_line(aes(x = date, y = ari)) + 
  geom_vline(xintercept = as.Date("2020-01-01"), linetype = "dashed", color = "red") +
  facet_wrap(~label, scales = "free_y") +
  labs(y = "ARI") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# ari scaled
# ggplot(data = D2) + 
#   geom_line(aes(x = date, y = ari_scaled)) + 
#   geom_vline(xintercept = as.Date("2020-01-01"), linetype = "dashed", color = "red") +
#   facet_wrap(~facility) + 
#   theme_minimal()

ggsave(filename = 'figures/Maryland_facilities_raw_05152025.pdf', width = 10, height = 7)

### Making the maps of locations
library(sf)
library(rnaturalearth)
library(patchwork)

# loading facility locs
D_fac <- read_csv('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/maryland_county_coordinates.csv')

# loading district shapefile
# # Read in GADM admin level 2 shapefile (this only showed two districts.)
# adm2 <- st_read("C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/Liberia_district_shapefiles/gadm41_LBR_2.shp")
# 
# # Filter to Maryland County
# maryland_districts <- adm2 %>% filter(NAME_1 == "Maryland")
adm2 <- st_read('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/Data/geoBoundaries-LBR-ADM2.geojson')

# Filter to Maryland County
maryland_districts <- adm2 %>% 
  filter(shapeName %in% c('Karluway #1',
                          'Karluway #2',
                          'Harper',
                          'Pleebo/Sodoken',
                          #'Whojah',
                          'Gwelekpoken',
                          'Nyorken'))


# Convert D_fac to sf using Longitude and Latitude
D_fac_sf <- D_fac %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 

# Get Liberia country outline
liberia <- ne_countries(scale = "medium", country = "Liberia", returnclass = "sf")

liberia_admin1 <- ne_states(country = "Liberia", returnclass = "sf")

# Filter Maryland County.
maryland <- liberia_admin1 %>% filter(name == "Maryland")

# Get the bounding box for Maryland.
#maryland_bbox <- st_bbox(maryland_districts)
maryland_bbox <- st_bbox(maryland)
bbox_poly <- st_as_sfc(maryland_bbox)


# Create p1: Liberia with Maryland highlighted and bounding box
p1 <- ggplot() +
  geom_sf(data = liberia, fill = "white", color = "black") +
  #geom_sf(data = maryland_districts, fill = "lightblue", color = "blue") +
  geom_sf(data = maryland, fill = 'lightblue') +
  geom_sf(data = bbox_poly, fill = NA, color = "red", size = 3) +
  labs(title = "Liberia") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Create p2: Maryland with facility points
p2 <- ggplot() +
  geom_sf(data = maryland, fill = "white", color = "black") +
  #geom_sf(data = maryland_districts, fill = 'white', color = 'black') +
  geom_sf(data = D_fac_sf, color = 'red', size = 2) +
  geom_sf(data = bbox_poly, fill = NA, color = "red", size = 3) + 
  labs(title = "Maryland County") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# # option with districts colored
# p2 <- ggplot() +
#   geom_sf(data = maryland, fill = "white", color = "black") +
#   geom_sf(data = D_fac_sf, aes(color = `Health  District`), size = 2) +
#   geom_sf(data = bbox_poly, fill = NA, color = "red", size = 3) + 
#   labs(title = "Maryland County") +
#   theme_minimal() +
#   theme(
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank(),
#     legend.position = "none"
#   )


# Combine plots
p1 + p2

# ggsave(filename = 'figures/Liberia_Maryland_map_nodist_05152025.pdf', width = 6, height = 4)

#
##### Theta values from negative binomial #####
D <- read.csv('data/coef_nb_1_3.csv')
thetas <- D[,9]
thetas <- thetas[!is.na(thetas)]
plot(density(thetas, na.rm = T))

tt <- fitdistrplus::fitdist(thetas, distr = 'gamma', method = 'mle')

plot(tt)

rgamma(10, shape = 2.5, rate = 1/3)


#
##### Looking at district sizes #####
district_count = D %>% 
  filter(date < '2020-01-01') %>%
  group_by(district) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            n = length(unique(facility)))

district_count$n[district_count$n > 20] <- 20
ggplot(district_count, aes(x = n)) + 
  geom_histogram(breaks = c(0,5, 10, 15, 20)) +
  scale_x_continuous(labels = c(0, 5, 10, 15, '20+')) + 
  ylab('# of districts') + 
  xlab('district size (number of facilities)') + 
  ggtitle('District sizes in Liberia')

#
##### Looking at variance across seasonal terms #####
load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/all_indicators_all_facilities_12012023.RData")
# now looking at the seasonal terms
apply(res_df, 2, median, na.rm =T)
apply(res_df, 2, mean, na.rm =T)
# ok so the mean is non-zero. Interesting. But the median is basically 0. So that's promising that this is right
apply(res_df, 2, sd, na.rm =T)
# woah. High standard deviation. 
plot(density(res_df[,7], na.rm = T))
tmp = res_df[,7]; tmp[tmp < -10 | tmp > 10] = 0
plot(density(tmp, na.rm = T))

# par(mfrow = c(1,1))
for(i in 1:8){
  plot(density(res_df[,i], na.rm = T))
}
# ah extreme outliers. This isn't really what I want though. I want the variance across facilities given a specific indicator rather than the various across indicators.
# "indicator_count_ari_total"
# "indicator_count_allvisits"
# "indicator_count_pneumo_1"
tt = res[["indicator_count_pneumo_1"]]
cor(tt[,3:10], use = 'complete.obs')
# ok so the seasons are highly correlated. Not surprising. That's what I would expect. Not sure how to simulate that properly, but that's ok.

for(i in 3:10){
  plot(density(tt[,i], na.rm = T))
}
# ok so a little more balanced

for(i in 3:10){
  print(names(tt)[i])
  col = na.omit(tt[,i])
  print(sd(col))
  x <- remove_outliers(col)
  print(sd(x))
  print(mean(x))
  print(median(x))
}
# standard deviation is around .11 - .48 (.07 - 0.13 excluding outliers) for indicator_count_allvisits

# standard deviation is around .11 - .48 (.07 - 0.13 excluding outliers) for indicator_count_ari_total 5-20 (1.2 - 5.5 excluding outliers)

# HERE! REVISIT the outliers - NEED TO CHECK DENSITY PLOTS WITHOUT THESE - AM I REALLY EXCLUDING OUTLIERS?

#
##### Fitting all the indicators #####

indicators <- grep('indicator', colnames(D2), value = T)

system.time({
res <- lapply(indicators, function(ind){
  tmp <- fit_WF_model(D2, outcome = ind, family = 'nb')
  tmp
})
})
names(res) <- indicators
# k it'll be about 30 minutes to run them all. Exactamundo

# res_df <- t(sapply(res, function(xx){colMeans(xx[3:10], na.rm = T)}))

res_df <- t(sapply(res, function(xx){apply(xx[,3:10], 2, function(tt){median(tt, na.rm = T)})}))

# save(res, res_df, file = "C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Syndromic Surveillance/results/all_indicators_all_facilities_12012023.RData")

col = res_df[,1]
col[col<= -50] <- -50
col[col>50] <- 50
plot(density(col, na.rm = T))

max(col[col<50], na.rm = T) # 37

res_df[98,]
# ok so super steep slope

col = res_df[,1]
col[col<= -10] <- -10
col[col>10] <- 10
plot(density(col, na.rm = T))

# all visits - centered around 6 but the year is smaller than I have done
res_df["indicator_count_allvisits",]
res_df["indicator_denom",]
res_df["indicator_count_ari_total",]

# wow ARI total is small! Damn really? I could implement these numbers...that's tiny.
tt <- res[["indicator_count_ari_total"]]
cor(na.omit(tt)[,3:10])
# so intercept and year are negatively correlated. Interesting. So the higher the year, the more the negative correlation.

E4 = tt[,3] + 4*tt[,4]
plot(density(E4[E4>-5 & E4 < 10], na.rm = T))
# centered around 3. Huh

#

##### Comparing data sources #####
dim(D)
dim(D2)

D2 = D2[,colnames(D)]
D2 = D2 %>%
  filter(facility %in% unique(D$facility),
         date <= max(D$date))

res <- fit_WF_model(D)
res2 <- fit_WF_model(D2)

#
##### Fitting all facilities WF/freqGLM #####

# fitting the mean model
res <- fit_WF_model(D)

colMeans(abs(res[,-1]), na.rm = T)
colMeans(res[,-1], na.rm = T)

res2 <- res %>%
  na.omit()

colMeans(abs(res2[,-1]), na.rm = T)
colMeans(res2[,-1], na.rm = T)

res3 <- res2 %>%
  filter(num_missing < 20)

colMeans(abs(res3[,-1]), na.rm = T)
colMeans(res3[,-1], na.rm = T)

res4 <- res3 %>% 
  filter(`(Intercept)` > 0)

colMeans(abs(res4[,-1]), na.rm = T)
colMeans(res4[,-1], na.rm = T)

apply(res4[,-1], 2, sd)

write.csv(res4, file = 'results/all_facility_betas_filtered_09052022.csv', row.names = F)

tt <- read.csv('results/all_facility_betas_filtered_09052022.csv')
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
            ari_miss = mean(is.na(indicator_count_ari_total)),
            ari_miss_count = sum(is.na(indicator_count_ari_total)))

district_miss = D %>% 
  group_by(district) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            n = length(unique(facility)))

county_miss = D %>% 
  group_by(county) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            n = length(unique(facility)))


hist(facility_miss$ari_miss, main = 'facility ARI missingness', xlab = 'proportion missing')
abline(v = 0.2, col = 'red')


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

##### Analyzing non-missing vs. missing datasets #####

facility_miss = D %>%
  filter(date < '2020-01-01') %>%
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            ari_miss_count = sum(is.na(indicator_count_ari_total)))

district_miss = D %>% 
  filter(date < '2020-01-01') %>%
  group_by(district) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            n = length(unique(facility)))


dist_miss <- district_miss %>% 
  filter(ari_miss > 0)

dist_comp <- district_miss %>% 
  filter(ari_miss == 0)

mean(dist_comp$n)

mean(dist_miss$n)

D_miss <- D %>%
  filter(date < '2020-01-01') %>%
  filter(district %in% dist_miss$district)

D_comp <- D %>%
  filter(date < '2020-01-01') %>%
  filter(district %in% dist_comp$district)

mean(D_miss$indicator_count_ari_total, na.rm = T)
mean(D_comp$indicator_count_ari_total)

fac_miss = facility_miss %>%
  filter(ari_miss > 0)

fac_comp = facility_miss %>%
  filter(ari_miss == 0)

fac_miss <- D %>%
  filter(date < '2020-01-01') %>%
  filter(facility %in% fac_miss$facility)

fac_comp <- D %>%
  filter(date < '2020-01-01') %>%
  filter(facility %in% fac_comp$facility)

mean(fac_miss$indicator_count_ari_total, na.rm = T)
mean(fac_comp$indicator_count_ari_total)

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
  # filter(ari_miss > 0.2) %>%
  # filter(denom_miss < 0.2, denom_miss > .05) %>%
  dplyr::select(district) %>% pull

par(mfrow = c(3,3))
set.seed(1)
for(d in sample(districts, 3)){
  print(d)
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

df = D %>% filter(district == "Somalia Drive District") %>%
  dplyr::select(date, district, facility, indicator_count_ari_total)

df_spread = tidyr::spread(df, facility, indicator_count_ari_total)

tmp2 = df_spread[,-c(1,2)]
for(col in colnames(tmp2)){
  tmp2[,col] = as.integer(is.na(tmp2[,col]))
}
tmp2 = as.matrix(tmp2)

rownames(tmp2) = as.character(df_spread$date)
heatmap(tmp2, keep.dendro = F, Rowv = NA, ylab = 'date', xlab = 'facilities', scale = 'none', labCol = '')

gplots::heatmap.2(tmp2, dendrogram = 'none', Rowv = NA, Colv = T, xlab = 'facilities', trace = 'none', labCol = '', key = F)


### ARI
districts = district_miss %>%
  filter(ari_miss > 0.3) %>%
  dplyr::select(district) %>% pull

par(mfrow = c(3,3))
for(d in sample(districts, 3)){
  print(d)
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

warning('I am storing the values incorrectly. Look at the spatial correlations in the next section')

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
##### Within county (not district) spatial correlation of missingness #####
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

##### Making R01 ARI spatial correlation figs #####

facility_miss = D %>% 
  group_by(facility) %>%
  summarize(denom_miss = mean(is.na(indicator_denom)),
            ari_miss = mean(is.na(indicator_count_ari_total)),
            ari_count = sum(!is.na(indicator_count_ari_total)))

exclude = facility_miss %>% filter(ari_count < 20) %>% select(facility) %>% pull()

D = add_periodic_cov(D)
uni_fac = unique(D$facility)

### raw values
res = NULL
cor_list = c()

for(d in unique(D$district)){
  # get the county
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D %>% filter(district == d, !(facility %in% exclude)) %>%
    select(date, district, facility, indicator_count_ari_total)
  df_spread = tidyr::spread(df, facility, indicator_count_ari_total)
  
  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

res_ARI_raw = res
cor_list_ARI_raw = cor_list

### get the deviance residuals
uni_fac = unique(D$facility)
D$ARI_resid = NA

# for each facility, run the harmonic model
test = lapply(uni_fac, function(xx){
  D3 = D %>% filter(facility == xx)
  
  if(sum(!is.na(D3$indicator_count_ari_total)) < 20){
    print(sprintf('skipping %s because of %s missing', xx, sum(is.na(D3))))
    return(NULL)
  }
  #print(xx)
  
  # create the formulas
  count_ARI_form = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  count_denom_form = as.formula("indicator_denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  
  # run the models
  D3$ARI_resid[!is.na(D3$indicator_count_ari_total)] = tryCatch({
    mod_1 = MASS::glm.nb(count_ARI_form, data = D3)
    summary(mod_1)$deviance.resid
  }, error = function(e){
    rep(NA, sum(!is.na(D3$indicator_count_ari_total)))
  })
  
  
  # only keep relevant columns  
  D3 = D3 %>%
    dplyr::select(date, county, district, facility, ARI_resid)
  
  # return the results
  D3
})

# convert to one data frame
D2 = do.call('rbind',test)

### ARI analyze them like you did with all the other values
res = NULL
cor_list = c()

for(d in unique(D2$district)){
  # get the county
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D2 %>% filter(district == d) %>%
    select(date, district, facility, ARI_resid)
  df_spread = tidyr::spread(df, facility, ARI_resid)
  
  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_list = c(cor_list, na.omit(as.vector(tmp)))
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

within_district_cors = res %>% group_by(county) %>%
  summarize(cor_within = mean(cor))

res_ARI_resid = res
cor_list_ARI_resid = cor_list

### plot 'em!
par(mfrow = c(1,2))
hist(cor_list_ARI_raw, main = 'within-district ARI values', xlab = 'correlation', breaks = 20, freq = F)
abline(v = mean(cor_list_ARI_raw), col = 'red')
mean(cor_list_ARI_raw) # 0.1288



ggplot(data = cor_list_ARI_raw) + geom_histogram()

hist(cor_list_ARI_resid, main = 'within-district ARI residuals', xlab = 'correlation', breaks = 20, freq = F)
abline(v = mean(cor_list_ARI_resid), col = 'red')
mean(cor_list_ARI_resid) # 0.0733


##### Spatial correlation of raw values #####

D = add_periodic_cov(D)
uni_fac = unique(D$facility)

### Indicator denom
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D$district)){
  # get the county 
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D %>% filter(district == d) %>%
    dplyr::select(date, district, facility, indicator_denom)
  df_spread = tidyr::spread(df, facility, indicator_denom)

  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

res_denom = res
cor_list_denom = cor_list

### ARI
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D$district)){
  # get the county
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D %>% filter(district == d) %>%
    select(date, district, facility, indicator_count_ari_total)
  df_spread = tidyr::spread(df, facility, indicator_count_ari_total)
  
  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

res_ARI = res
cor_list_ARI = cor_list

within_district_cors = res %>% group_by(county) %>%
  summarize(cor_within = mean(cor))

# plotting district means of correlations
par(mfrow = c(1,2))
hist(res_ARI$cor, main = 'within-district ARI', xlab = 'mean correlation')
# hist(res_denom$cor_denom, main = 'within-district denom', xlab = 'mean correlation')
# pretty high. Interesting.

# plotting facility correlations
hist(cor_list_ARI, main = 'within-district ARI values', xlab = 'correlation', breaks = 20)
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.102

hist(cor_list_denom, main = 'within-district denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) # 0.151

# interesting. So not as high as first glance. 
# However, this does show that there is clear correlation of deviance. That is a good sign moving forward for spatial correlation.


### Within-county, not district
res2 = NULL
cor_list_denom = c()
cor_list_ARI = c()

# sorry, but I'm shamefully doing for loops on for loops on for loops. It's slow because of that
for(cc in unique(D$county)){
  # spread the data for this district
  df = D %>% filter(county == cc) %>%
    select(date, county, district, facility, indicator_denom, indicator_count_ari_total)

  # the long way - a for loop!
  for(f in unique(df$facility)){
    target = df %>% filter(facility == f)
    
    
    d = target$district[1]
    compare_facs = df %>%
      filter(district != d) %>%
      pull(facility) %>% 
      unique()
    for(f2 in compare_facs){
      compare = df %>% filter(facility == f2)
      cor_d = cor(target$indicator_denom, compare$indicator_denom, use = 'pairwise.complete.obs')
      cor_ARI = cor(target$indicator_count_ari_total, compare$indicator_count_ari_total, use = 'pairwise.complete.obs')
      
      res2 = rbind(res2, data.frame(county = cc, district = d, cor_ARI = cor_ARI, cor_denom = cor_d))
      cor_list_denom = c(cor_list_denom, cor_d)
      cor_list_ARI = c(cor_list_ARI, cor_ARI)
    }
  }
}

without_district_cors = res2 %>% group_by(county) %>%
  summarize(cor_without = mean(cor_ARI, na.rm = T))

cor_list_denom = na.omit(cor_list_denom)
cor_list_ARI = na.omit(cor_list_ARI)

# plotting facility correlations
hist(cor_list_ARI, main = 'within-county ARI values', xlab = 'correlation')
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.103

hist(cor_list_denom, main = 'within-county denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) 
# 0.13

# merging two for comparison
comparison = merge(within_district_cors, without_district_cors) %>%
  arrange(cor_within)

# write.csv(comparison, row.names = F, file = 'results/ari_value_correlation_04192021.csv')

# save(comparison, res, res2, file = 'results/ari_value_correlation_04192021.RData')

##### Spatial correlation of deviance residuals #####

D = add_periodic_cov(D)

# D2 = D %>% filter(district == 'Dowein')
# uni_fac = unique(D2$facility)

uni_fac = unique(D$facility)
D$ARI_resid = D2$denom_resid = NA

# for each facility, run the harmonic model
test = lapply(uni_fac, function(xx){
  D3 = D %>% filter(facility == xx)
  
  if(sum(!is.na(D3$indicator_count_ari_total)) < 20 | sum(!is.na(D3$indicator_denom)) < 20){
    print(sprintf('skipping %s because of %s missing', xx, sum(is.na(D3))))
    return(NULL)
  }
  #print(xx)
  
  # create the formulas
  count_ARI_form = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  count_denom_form = as.formula("indicator_denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
  
  # run the models
  D3$ARI_resid[!is.na(D3$indicator_count_ari_total)] = tryCatch({
    mod_1 = MASS::glm.nb(count_ARI_form, data = D3)
    summary(mod_1)$deviance.resid
  }, error = function(e){
    rep(NA, sum(!is.na(D3$indicator_count_ari_total)))
  })
  
  D3$denom_resid[!is.na(D3$indicator_denom)] = tryCatch({
    mod_2 = MASS::glm.nb(count_denom_form, data = D3)
    summary(mod_2)$deviance.resid
  }, error = function(e){
    rep(NA, sum(!is.na(D3$indicator_denom)))
  })
  
  # mod_1 = MASS::glm.nb(count_ARI_form, data = D3)
  # mod_2 = MASS::glm.nb(count_denom_form, data = D3)
  # 
  # # get the deviance residuals
  # D3$ARI_resid[!is.na(D3$indicator_count_ari_total)] = summary(mod_1)$deviance.resid
  # D3$denom_resid[!is.na(D3$indicator_denom)] = summary(mod_2)$deviance.resid

  # only keep relevant columns  
  D3 = D3 %>%
    dplyr::select(date, county, district, facility, denom_resid, ARI_resid)
  
  # return the results
  D3
})

# convert to one data frame
D2 = do.call('rbind',test)

# analyze them like you did with all the other values

### Indicator denom
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D2$district)){
  # get the county
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D2 %>% filter(district == d) %>%
    dplyr::select(date, district, facility, denom_resid)
  df_spread = tidyr::spread(df, facility, denom_resid)
  
  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

res_denom = res
cor_list_denom = cor_list

### ARI
res = NULL
cor_list = c()

# sorry, but I'm doing a for loop
for(d in unique(D2$district)){
  # get the county
  county = D %>% filter(district == d) %>% distinct(county) %>% pull()
  
  # spread the data for this district
  df = D2 %>% filter(district == d) %>%
    select(date, district, facility, ARI_resid)
  df_spread = tidyr::spread(df, facility, ARI_resid)
  
  # get the correlation matrix across facilities
  tmp <- cor(df_spread[,-c(1,2)], use = 'pairwise.complete.obs')
  tmp[lower.tri(tmp, diag = T)] <- NA
  
  # store the correlations
  cor_list = c(cor_list, na.omit(as.vector(tmp)))
  
  # store the correlations
  cor_district = na.omit(as.vector(tmp))
  cor_list = c(cor_list, cor_district)
  
  # store the results with county and district
  res_district = data.frame(county = rep(county, length(cor_district)), district = rep(d, length(cor_district)), cor = cor_district)
  res = rbind(res, res_district)
}

within_district_cors = res %>% group_by(county) %>%
  summarize(cor_within = mean(cor))

res_ARI = res
cor_list_ARI = cor_list

# plotting district means of correlations
par(mfrow = c(1,2))
hist(res_ARI$cor_ARI, main = 'within-district ARI residuals', xlab = 'mean correlation')
hist(res_denom$cor_denom, main = 'within-district denom', xlab = 'mean correlation')
# pretty high. Interesting.

# plotting facility correlations
hist(cor_list_ARI, main = 'within-district ARI residuals', xlab = 'correlation')
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.074

hist(cor_list_denom, main = 'within-district denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) # 0.133

# interesting. So not as high as first glance. 
# However, this does show that there is clear correlation of deviance. That is a good sign moving forward for spatial correlation.


### Within-county, not district
res2 = NULL
cor_list_denom = c()
cor_list_ARI = c()

# sorry, but I'm shamefully doing for loops on for loops on for loops. It's slow because of that
for(cc in unique(D2$county)){
  # spread the data for this district
  df = D2 %>% filter(county == cc) %>%
    select(date, county, district, facility, denom_resid, ARI_resid)
  #df_spread = tidyr::spread(df, facility, indicator_denom)
  
  # the long way - a for loop!
  for(f in unique(df$facility)){
    target = df %>% filter(facility == f)
    
    # # skip if the facility has no missingness
    # if(sum(is.na(target)) == 0){
    #   next
    # }
    
    d = target$district[1]
    compare_facs = df %>%
      filter(district != d) %>%
      pull(facility) %>% 
      unique()
    for(f2 in compare_facs){
      compare = df %>% filter(facility == f2)
      # cor_list_denom = c(cor_list_denom, cor(target$denom_resid, compare$denom_resid, use = 'pairwise.complete.obs'))
      # cor_list_ARI = c(cor_list_ARI, cor(target$ARI_resid, compare$ARI_resid, use = 'pairwise.complete.obs'))
      
      compare = df %>% filter(facility == f2)
      cor_d = cor(target$denom_resid, compare$denom_resid, use = 'pairwise.complete.obs')
      cor_ARI = cor(target$ARI_resid, compare$ARI_resid, use = 'pairwise.complete.obs')
      
      res2 = rbind(res2, data.frame(county = cc, district = d, cor_ARI = cor_ARI, cor_denom = cor_d))
      cor_list_denom = c(cor_list_denom, cor_d)
      cor_list_ARI = c(cor_list_ARI, cor_ARI)
      
    }
  }
}

without_district_cors = res2 %>% group_by(county) %>%
  summarize(cor_without = mean(cor_ARI, na.rm = T))

# merging two for comparison
comparison = merge(within_district_cors, without_district_cors) %>%
  arrange(cor_within)

# write.csv(comparison, row.names = F, file = 'results/ari_residual_correlation_04192021.csv')

# save(comparison, res, res2, file = 'results/ari_residual_correlation_04192021.RData')

cor_list_denom = na.omit(cor_list_denom)
cor_list_ARI = na.omit(cor_list_ARI)

# plotting facility correlations
hist(cor_list_ARI, main = 'within-county ARI residuals', xlab = 'correlation')
abline(v = mean(cor_list_ARI), col = 'red')
mean(cor_list_ARI) # 0.065

hist(cor_list_denom, main = 'within-county denom', xlab = 'correlation')
abline(v = mean(cor_list_denom), col = 'red')
mean(cor_list_denom) 
# 0.13

##### ARI correlation not within county #####

N = 5000
facility_data = D %>% select(county, district, facility) %>% distinct()

set.seed(1)
cor_ARI = c()
for(i in 1:N){
  f_sample = sample_n(facility_data, 2)
  # if within the same county, skip it
  if(f_sample$county[1] == f_sample$county[2]){
    next
  }else{
    # get the residual correlations
    target = D %>% filter(facility == f_sample$facility[1])
    compare = D %>% filter(facility == f_sample$facility[2])
    #cor_d_resid = cor(target_resid$denom_resid, compare_resid$denom_resid, use = 'pairwise.complete.obs')
    cor_ARI = c(cor_ARI, cor(target$indicator_count_ari_total, compare$indicator_count_ari_total, use = 'pairwise.complete.obs'))
  }
}

hist(cor_ARI, main = 'without-county ARI values', xlab = 'correlation', breaks = 20)
abline(v = mean(cor_ARI, na.rm = T), col = 'red')
mean(cor_ARI, na.rm = T) # 0.14

#### Residual ARI correlation not within county
# for sampling
N = 2000
facility_data = D2 %>% select(county, district, facility) %>% distinct()

cor_ARI_resid = c()
for(i in 1:N){
  f_sample = sample_n(facility_data, 2)
  # if within the same county, skip it
  if(f_sample$county[1] == f_sample$county[2]){
    next
  }else{
    # get the residual correlations
    target_resid = D2 %>% filter(facility == f_sample$facility[1])
    compare_resid = D2 %>% filter(facility == f_sample$facility[2])
    #cor_d_resid = cor(target_resid$denom_resid, compare_resid$denom_resid, use = 'pairwise.complete.obs')
    cor_ARI_resid = c(cor_ARI_resid, cor(target_resid$ARI_resid, compare_resid$ARI_resid, use = 'pairwise.complete.obs'))
  }
}

hist(cor_ARI_resid, main = 'without-county ARI values', xlab = 'correlation', breaks = 20)
abline(v = mean(cor_ARI_resid, na.rm = T), col = 'red')
mean(cor_ARI_resid, na.rm = T) # 0.14

# 
##### ICC #####

tmp = D %>% filter(!is.na(indicator_count_ari_total))
clus.rho(D$indicator_count_ari_total, D$district)

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

##### Plotting different facilites with ~ 20% missing
load('results/simulation_epi_MCARp2_R500_res_07142022.RData')

df <- imputed_list[[1]] %>%
  filter(facility == 'A1')
df$y2 <- df$y_pred_freqGLMepi
df$y2[!is.na(df$y)] <- NA

plot(df$y_true, type = 'l')
plot(df$y, type = 'l')
plot(df$y_true, type = 'l')

plot(df$y, type = 'l')
points(df$y2, col = 'red')





