library(rstanarm)
library(bayestestR)
library(easystats)
library(insight)



##### Testing 1 - using iris data #####
model <- lm(Sepal.Length ~ Petal.Length, data=iris)
summary(model)
insight::get_parameters(model)

model <- stan_glm(Sepal.Length ~ Petal.Length, data=iris)#, iter = 2000)
describe_posterior(model)
posteriors <- insight::get_parameters(model)

head(posteriors)

test <- posterior_predict(model, newdata = iris)
predictions_I_think <- colMeans(test)
plot(iris$Sepal.Length, predictions_I_think)
# ok good. This function works as expected. 

4.392924 + 1.4*0.3947724 -> aa
any(abs(test[,1] - aa) < 1e-5) # hmm. Ah it's predicted ya dingus. This should include variance!

### With the data I am now using
setwd('C:/Users/nickl/Documents/global_covid19_response/')
source("R/model_functions.R")
source("R/model_figures.R")
data <- readRDS("data/data_example_singlecounty.rds")

df = data %>% 
  filter(facility == 'Facility B',
         date <  as.Date("2020-01-01")) %>%
  na.omit()

period = 12
df = df %>%
  dplyr::mutate(year = year(date) - min(year(date)) + 1,
                month = month(date),
                cos1 = cos(2*1*pi*month/period),
                sin1 = sin(2*1*pi*month/period),
                cos2 = cos(2*2*pi*month/period),
                sin2 = sin(2*2*pi*month/period),
                cos3 = cos(2*3*pi*month/period),
                sin3 = sin(2*3*pi*month/period))

# OG quasipoisson
formula_col = as.formula("indicator_denom ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
mod_freq <- glm(formula_col, data = df, family=quasipoisson)
summary(mod_freq)

# bayes negative binomial now
mod_bayes <- stan_glm(formula_col, data = df, family = neg_binomial_2)
# can also do stan_glmer
describe_posterior(mod_bayes)
colMeans(insight::get_parameters(mod_bayes))
# ok good it's basically the same. How to predict though?

ss = predict(mod_freq, data = df, type = 'response')
tt = predict(mod_bayes, data = df, type = 'response')
cor(ss,tt)
cor(ss, df$indicator_denom)

plot(df$indicator_denom, type = 'l')
lines(ss, col = 'red')
lines(tt, col = 'blue')

# Ok so they match. The predictions are not that good though. Yikes.

# I didn't mess with the priors at all though.

# next up
  # mess with priors?
  # implement Bayesian imputation
  # implement imputation on all facilities and check results?
  # Make a darn report! Are ya darn tootin!



##### Testing 2 - using Liberia data #####

# function to add periodic terms
add_periodic_cov <- function(df, period = 12){
  df = df %>%
    dplyr::mutate(year = year(date) - min(year(date)) + 1,
                  month = month(date),
                  cos1 = cos(2*1*pi*month/period),
                  sin1 = sin(2*1*pi*month/period),
                  cos2 = cos(2*2*pi*month/period),
                  sin2 = sin(2*2*pi*month/period),
                  cos3 = cos(2*3*pi*month/period),
                  sin3 = sin(2*3*pi*month/period))
  return(df)
}

setwd('C:/Users/nickl/Documents/global_covid19_response/')

D = readRDS('data/liberia_cleaned_NL.rds')
D = D %>% filter(district == "Harper District" )

# 8 facilities. Good

D = add_periodic_cov(D, period = 12) %>% as.data.frame()

formula_col = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")
formula_glmm = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + (1|facility)")

tmp = D %>% filter(facility == 'Pullah Clinic')

# basic models
stan_glm(formula_col, data = tmp, family = neg_binomial_2, refresh = 0, iter = 500)

# glmm model (without the structured covariance)
glmer_fit <- stan_glmer(formula_glmm, data = D, family = neg_binomial_2, refresh = 0, iter = 500)

varCorr(glmer_fit)

# how to specify auto-correlation?

# How to do structured spatial random effects?
# Conditional auto-regressive something?
  # Maybe I should learn about this next.


# Messing with GlMMTMB (no Bayesian stuff RN)
head(D)
uni_fac = unique(D$facility)
locations = data.frame(facility = uni_fac, 
                       x = rnorm(length(uni_fac),10,5), 
                       y = rnorm(length(uni_fac), 10,5))
locations$pos = numFactor(scale(locations$x), scale(locations$y))

D2 = merge(D, locations)

library(glmmTMB)

model_1_fit <- glmmTMB(formula_glmm, data = D2, family = nbinom2)

formula_mat = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3 + mat(pos + 0|facility)")

model_2_fit <- glmmTMB(formula_mat, data = D2, family = nbinom2)
# so model convergeence. But maybe it still works?

# Ok so I got this to run. But not based on true distance. That's un problemo

# I want to be able to use conditional autoregression priors for the structured spatial random effects


### CARBayesST
library(CARBayesST)

# Make a W binary based off of district

# To try
#   ST.CARlinear: basic approach - this doesn't allow for covariates? But that's ok to start, I guess.
#   ST.CARar: allows autoregressive order 1
#   ST.CARadaptive: allows estimation of adjacency, W. That might be useful in our case. As in, allows for some nearby locations to be highly correlated and other to not be. So this adds an extra parameter to estimate that.
#   ST.CARclustrends: ?? Maybe. Not for now
#   Ah freak I might have to code my own function. Dasss annoying. 
# 
#   Huh I can also estimate the adjacency matrix using W.estimate(). Seems cool (this removes edges rather than adds - so i could try including full counties here)
# https://cran.r-project.org/web/packages/CARBayesST/vignettes/CARBayesST.pdf

library(CARBayesdata)
library(sp)
data("GGHB.IG")
data('pollutionhealthdata')
head(pollutionhealthdata)

# IG = Different locations

library(dplyr)
pollutionhealthdata <- pollutionhealthdata %>% mutate(
  SMR = pollutionhealthdata$observed / pollutionhealthdata$expected,
  logSMR = log(pollutionhealthdata$observed / pollutionhealthdata$expected))

library(GGally)
ggpairs(pollutionhealthdata, columns=c(9, 5:7))

group_IG <- group_by(pollutionhealthdata, IG)
SMR.av <- summarise(group_IG, SMR.mean = mean(SMR))
GGHB.IG@data$SMR <- SMR.av$SMR.mean

library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = SMR.av$IG)
W.list <- nb2listw(W.nb, style = "B")
W <- nb2mat(W.nb, style = "B")
head(W[,1:5])
diag(W)

formula <- observed ~ offset(log(expected)) + jsa + price + pm10

chain1 <- ST.CARar(formula = formula, family = "poisson",
                   data = pollutionhealthdata, W = W, burnin = 20000, n.sample = 220000,
                   thin = 100, AR=1)
# took 2.5 minutes

# using my data!
D = readRDS('data/liberia_cleaned_NL.rds')
D = add_periodic_cov(D, period = 12) %>% as.data.frame()
D = D %>% filter(county == "Bomi" )
# create adj by district

# create adjacency matrix!
D2 = D %>% dplyr::select(district, facility) %>% distinct()

W = full_join(D2,D2, by = 'district') %>%
  filter(facility.x != facility.y) %>%
  dplyr::select(-district) %>%
  igraph::graph_from_data_frame() %>%
  igraph::as_adjacency_matrix() %>%
  as.matrix()

# model formula
formula_1 = as.formula("indicator_count_ari_total ~ year + cos1 + sin1 + cos2 + sin2 + cos3 + sin3")

chain1 <- ST.CARar(formula = formula_1, family = "poisson",
                   data = D, W = W, burnin = 20000, n.sample = 40000,
                   thin = 10, AR=1)

# so it ran! Got that. Now what to do with it!
test = chain1$fitted.values
plot(D$indicator_count_ari_total, chain1$fitted.values)
# almost exact, huh? Overfitting much?



### To look into
# Different ST.CAR functions
# Scaling - is ST.CAR being dominated by large facilities?
# Does it make sense to use W.estimate() here?