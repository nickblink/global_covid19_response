library(rstanarm)
library(bayestestR)
library(easystats)
library(insight)

model <- lm(Sepal.Length ~ Petal.Length, data=iris)
summary(model)
insight::get_parameters(model)

model <- stan_glm(Sepal.Length ~ Petal.Length, data=iris)#, iter = 2000)
describe_posterior(model)
posteriors <- insight::get_parameters(model)

head(posteriors)

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






