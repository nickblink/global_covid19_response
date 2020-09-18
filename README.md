README
================

# Global COVID-19 Response

*Last updated: 15 Sep 2020* \#\# Table of Contents - [About](#About) -
[Goals](#Goals) - [Modelling technique](#Modelling-technique) -
[Facility-level models](#Facility-level-models) - [District and
county-level models](#District-and-county-level-models) - [Missing data
considerations](#Missing-data-considerations) - [Overview of folders and
files](#Overview-of-folders-and-files)

## About:

This repository contains code to follow the Data Processing Pipeline for
the Global Covid-19 Syndromic Surveillance Team - a partnership between
sites at Partners in Health, the Global Health Research Core at Harvard
Medical School, and Brigham and Women’s Hospital. The data has been
modified to respect the privacy of our sites, in hopes that other groups
can benefit from the functions we have written.

This repository contains data, code, and other items needed to reproduce
this work. Outputs include figures, tables, and Leaflet maps. Further
explanation of outputs and their construction is given in the “Overview
of folders and files” section, which includes detailed explanations of
the functions we have written.

## Goals:

The main goal of the Global COVID-19 Syndromic Survillance Team is to
monitor changes in indicators that may signal changes in COVID-19 case
numbers in health systems from our eight partnering countries: Haiti,
Lesotho, Liberia, Madagascar, Malawi, Mexico, Peru, and Rwanda. This is
accomplished through establishing a baseline using prior data, and
monitoring for deviations for relevant indicators. The data
visualization tools created using our functions allow identification of
local areas that are experiencing upticks in COVID-19-related symptoms.

## Modelling technique:

The process starting with the raw data and finishing with the various
outputs is referred to as the Data Processing Pipeline (see Figure 1
below):

After data has been cleaned, it is processed according to the level it
is available at (either on a facility of county/district basis) for each
indicator. This is done by taking data from a historic baseline period,
and then projecting it into the evaluation period. This then is compared
to the observed counts/proportions. A 95% confidence interval has been
chosen, and we have defined the baseline period to be data from January
2016-December 2019.

The functions included in this repository focus on the modelling and
processing stages.

### Facility-level models:

The following generalized linear model using the negative binomial
distribution with log link to account for overdispersion was used:

Where time, t, is defined as the cumulative month number (e.g. January
2016 is month 1). The year term captures trend, and the harmonic term,
k, captures seasonality. Note that year is a linear term and will only
capture monotonic trends. In theory, one could model each year with a
binary indicator to allow for flexible yearly deviations. We were unable
to do this because our baseline period does not contain any 2020 months,
which would render expected counts during this period unidentifiable.

### District and county-level models:

The following generalized linear mixed model was used to model the
expected counts at the district or county level. We included random
intercepts for each facility (denoted by j). The negative binomial
distribution with log link was used to account for overdispersion:

Importantly, we did not report facility-level estimates using this
model. Instead, we estimated the marginal (population-level) total count
by integrating over the random effect distribution. The facility-level
estimates from this model will NOT match results from the individual
facility-level models above. The facility-level models allowed for
facility-specific year and seasonality trends, whereas the GLMM assumed
common year and seasonality trends across all facilities in the district
or county, allowing only the baseline counts to vary by facility via the
random intercepts.

### Missing data considerations:

We excluded facilities from our analysis for two reasons: (1) missing
dates in the baseline period (creation of the expected counts model) (2)
missing observed counts in the evaluation period.

For the first reason, facilities with high levels of missing data (more
than 20% of baseline dates missing) were excluded. Although there are
statistical methods that can handle missing time series data, we decided
to only include sites that demonstrated ability to collect data over
time. A complete case (time) analysis was conducted with included
facilities, which assumes that the counts were missing completely at
random (MCAR). Specifically, we assumed the reason for missing counts
was independent of time and count value. If the MCAR assumption was
violated and we had information about the missing data mechanism, one
could impute values for the missing data and report unbiased expected
counts and correct inference.

For the second reason, facilities with ANY missing monthly data during
the evaluation period (January 2020 onward) were removed. As the
syndromic surveillance exercise hinges on comparing the observed counts
to the expected and flagging for deviations, we require complete
observed data during this period. In this context, it would be invalid
to impute observed counts based on information from the baseline period.
In theory, one could attempt to impute the observed count based on
information during the evaluation period.

## Overview of folders and files:

### data

This folder contains example data used to demonstrate functions.
\#\#\#\# data.example\_singlecounty.rds The facility-level dataset used
to demonstrate the functions throughout this repository. Note- specific
names and numbers have been altered to respect the privacy of our sites.

### R

This folder contains the functions used to create the key data
visualization figures and maps.

### maps

### figures

This folder contains figures that have been included in README.md.
