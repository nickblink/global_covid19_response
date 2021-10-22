

library(surveillance)


##### Following the vignette #####
data('influMen')
print(fluMen <- disProg2sts(influMen))

plot(fluMen, type = observed ~ time | unit,
     same.scale = F,
     col = 'grey')


# reading in data for German flu counts.
flu.counts <- as.matrix(read.table(system.file("extdata/counts_flu_BYBW.txt",
                                               package = "surveillance"),
                                   check.names = F))

# neighborhood matrix. Obviously to look at later.
nhood <- as.matrix(read.table(system.file("extdata/neighbourhood_BYBW.txt",
                                          package = "surveillance"),
                              check.names = F))

library(Matrix)
print(image(Matrix(nhood)))

popfracs <- read.table(system.file("extdata/population_2001-12-31_BYBW.txt",
                                   package = "surveillance"),
                       header = T)$popFrac

# create the sts object that this package uses. I don't think the popfracs part is necessary. At least I hope, because I am not using that.
flu <- sts(flu.counts, start = c(2001, 1), frequency = 52, population = popfracs, neighbourhood = nhood)


### new data
data("measlesDE")
aggregate(measlesDE, nfreq = 26)
meningo <- fluMen[,'meningococcus']

f1 <- addSeason2formula(f = ~ 1, S = 1, period = 52)

# end is the formula for nu, the seasonal effects
result0 <- hhh4(meningo, control = list(end = list(f = f1), family = 'Poisson'))
summary(result0)

# 
a = result0$fitted.values[,1]
b = meningo@observed[,1]
cor(a, b[-1])
cor(a, b[-312])


result1 <- update(result0, family = 'NegBin1')

# the ar = list(f = ~-1) default excludes the autocorrelation parameter, but this line of code, includes it
result2 <- update(result1, ar = list(f = ~1))

result3 <- update(result2, ne = list(f = ~1, weights = neighbourhood(meningo)))

# I think currently this does not include the time trend (actually the written out model has a full time trend - does ours just have a yearly trend? Why not just a full time one?)

plot(result2)
# cool


### Flu and meningitis together
neighbourhood(fluMen)["meningococcus",'influenza'] <- 0
neighbourhood(fluMen)

f.end <- addSeason2formula(f = ~ -1 + fe(1, unitSpecific = T), S = c(3,1), period = 52)

m <- list(ar = list(f = ~-1 + fe(1, unitSpecific = T)),
          ne = list(f = ~ 1, weights = neighbourhood(fluMen)),
          end = list(f = f.end), 
          family = 'NegBinM')  # disease specific overdispersion
result <- hhh4(fluMen, control = m)

summary(result)

plot(result, units = 1:2)

# random effects of location-specific intercepts I believe, as well as linear time trend
f.end <- addSeason2formula(f = ~ -1 + ri(type = 'iid', corr = 'all') + I((t - 208)/100), S = 3, period = 52)

data('fluBYBW')
dim(fluBYBW) # 8 years, 140 districts

model.B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type = 'iid', corr = 'all'),
                           weights = neighbourhood(fluBYBW),
                           normalize = T),
                 end = list(f = f.end, offset = population(fluBYBW)),
                 family = 'NegBin1', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.B2 <- hhh4(fluBYBW, model.B2)


f.end <- addSeason2formula(f = ~ -1 + ri(type = 'iid', corr = 'all') + I(t), S = 3, period = 52)

data('fluBYBW')
dim(fluBYBW) # 8 years, 140 districts

model.B3 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type = 'iid', corr = 'all'),
                           weights = neighbourhood(fluBYBW),
                           normalize = T),
                 end = list(f = f.end, offset = population(fluBYBW)),
                 family = 'NegBin1', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.B3 <- hhh4(fluBYBW, model.B3)

### How does this deal with missingness?
f.end <- addSeason2formula(f = ~ -1 + ri(type = 'iid', corr = 'all') + I((t - 208)/100), S = 3, period = 52)

data('fluBYBW')
dim(fluBYBW) # 8 years, 140 districts

# make missing
p = 0.2
fluBYBW@observed[sample(length(fluBYBW@observed), round(length(fluBYBW@observed)*p))] <- NA
sum(is.na(fluBYBW@observed))

model.B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type = 'iid', corr = 'all'),
                           weights = neighbourhood(fluBYBW),
                           normalize = T),
                 end = list(f = f.end, offset = population(fluBYBW)),
                 family = 'NegBin1', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.B2 <- hhh4(fluBYBW, model.B2)

# interesting. So it still works. And these results are slightly different. So how is it doing it though? Especially with the neighborhood correlation

# looking at the documentation, I think this counts any NA as an NA for all neighbors. No bueno. And then I think it just excludes all of these. Ya that's not good.

##### So now how does it do predictions on missing data? #####
data('fluBYBW')
dim(fluBYBW) 

p = 0.3
fluBYBW.M = fluBYBW
fluBYBW.M@observed[sample(length(fluBYBW.M@observed), round(length(fluBYBW.M@observed)*p))] <- NA

f.end <- addSeason2formula(f = ~ -1 + ri(type = 'iid', corr = 'all') + I((t - 208)/100), S = 3, period = 52)

model.B1 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type = 'iid', corr = 'all'),
                           weights = neighbourhood(fluBYBW),
                           normalize = T),
                 end = list(f = f.end, offset = population(fluBYBW)),
                 family = 'NegBin1', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.B1 <- hhh4(fluBYBW, model.B1)

### No autoregression
model.B2 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ 1,
                           weights = neighbourhood(fluBYBW.M),
                           normalize = T),
                 end = list(f = f.end, offset = population(fluBYBW.M)),
                 family = 'NegBin1', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.B2 <- hhh4(fluBYBW.M, model.B2)


sum(is.na(result.B2$fitted.values))
sum(is.na(result.B1$fitted.values))
sum(is.na(fluBYBW.M@observed))
# interesting. Not sure what causing the missingness then. How are there more missing in the observed values than the fitted? But almost the same?

# ^But when I remove the autoregressive term in the model there is no missingness, so all missingness comes from that I guess. 

a = which(is.na(result.B2$fitted.values))
b = which(is.na(fluBYBW.M@observed))
length(intersect(a,b))
# so not even much of an intercept. Only 1/3. Weird yo. Then is it all based on neighbors?
# ah the fitted values are missing one row. Which? Probably the first one, right?

a = result.B2$fitted.values[,1]
b = fluBYBW.M@observed[,1]

cor(a, b[-1], use = 'complete.obs')
cor(a, b[-416], use = 'complete.obs')
# what????? I need to read the code. It should be -1 because the first data point can't be used.

test = fluBYBW.M@observed[-416,]
identical(is.na(test), is.na(result.B2$fitted.values)) #T
# So all in all they dont predict the last row and the only missingness is of the original data points. Got it. Still not sure why the last row (as opposed to the first) is not used and also what they do about missingness in spatially related areas

# it looks like it's probably not predicting the last value. Not sure why here, but oh well.
tp = oneStepAhead(result.B2, tp = nrow(fluBYBW))

# also, if the NA values only match the original NA values (not those of neighbors), how do they deal with neighboring NAs? Need to look this up.
##### Now using our data #####
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source('R/imputation_functions.R')

D <- simulate_data_spatiotemporal(district_sizes = c(4), R = 1, rho = 0.3, alpha = 0.5, tau = 0.5)[[1]][[1]]

D2 = D %>% dplyr::select(district, facility) %>% distinct()
W = full_join(D2, D2, by = 'district') %>%
  filter(facility.x != facility.y) %>%
  dplyr::select(-district) %>%
  igraph::graph_from_data_frame() %>%
  igraph::as_adjacency_matrix() %>%
  as.matrix()

ARI.counts <- D %>% 
  dplyr::select(date, facility, y) %>%
  tidyr::spread(facility,y) %>% 
  arrange(date) %>%
  dplyr::select(-date) %>%
  as.matrix()

## ARI.counts is a T x I matrix of observations. W is an I x I matrix of adjacency.
ARI <- sts(ARI.counts, start = c(2016, 1), frequency = 12, neighbourhood = W)

f.end <- addSeason2formula(f = ~ 1 + I(t/12), S = 3, period = 12)

model.1 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ 1,
                           weights = neighbourhood(ARI),
                           normalize = T),
                 end = list(f = f.end),
                 family = 'Poisson', verbose = T,
                 optimizer = list(variance = list(method = 'Nelder-Mead')))

result.1 <- hhh4(ARI, model.1)
  
coef(result.1, se = T,
     idx2Exp = T)
# crazy coefficients for the seasonal terms. That's a bit concerning but ok.

f.end <- addSeason2formula(f = ~ -1 + I(t/12) + fe(1, unitSpecific = T), S = 3, period = 12)

model.2 <- list(ar = list(f = ~ 1),
                ne = list(f = ~ 1,
                          weights = neighbourhood(ARI),
                          normalize = T),
                end = list(f = f.end),
                family = 'Poisson', verbose = T,
                optimizer = list(variance = list(method = 'Nelder-Mead')))

result.2 <- hhh4(ARI, model.2)

coef(result.2, se = T,
     idx2Exp = T)


f.end <- addSeason2formula(f = ~ -1 + fe(I(t/12), unitSpecific = T) + fe(1, unitSpecific = T), S = 3, period = 12)

model.3 <- list(ar = list(f = ~ 1),
                ne = list(f = ~ 1,
                          weights = neighbourhood(ARI),
                          normalize = T),
                end = list(f = f.end),
                family = 'Poisson', verbose = T,
                optimizer = list(variance = list(method = 'Nelder-Mead')))

result.3 <- hhh4(ARI, model.3)

coef(result.3, se = T,
     idx2Exp = T)

# gosh darnit I just need to do this by hand now.
# ^ Why? I can't do seasonal terms by the individual facility?
# How are they getting the standard errors of these estimates?

f.end <- addSeason2formula(f = ~ -1 + fe(I(t/12), unitSpecific = T) + fe(1, unitSpecific = T), S = c(3,3,3,3), period = 12)

model.4 <- list(ar = list(f = ~ 1, f = ~ 1, f = ~ 1, f = ~ 1, unitSpecific = T),
                ne = list(f = ~ 1,
                          weights = neighbourhood(ARI),
                          normalize = T),
                end = list(f = f.end),
                family = 'Poisson', verbose = T,
                optimizer = list(variance = list(method = 'Nelder-Mead')))

result.4 <- hhh4(ARI, model.4)
result.4

# can I do ar or ne terms by facility?
# How do they optimize?
# How do they find the variance? Some Fisher estimation


