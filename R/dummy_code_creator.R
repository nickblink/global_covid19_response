### This script will create dummy data for the Github repo ###

#load libraries
library(dplyr)

#read in original data
dummy <- readRDS("../../global_covid19_ss/liberia/data/cleaned/liberia_cleaned_08-25-2020.rds")

#Change names/perturb numbers (note- only filtering for Maryland county)
dummy <- dummy %>% filter(county== "Maryland") %>%
  dplyr::select(date, county, district, facility, starts_with("indicator_count_ari"), indicator_denom, indicator_denom_under5, indicator_denom_over5) %>%
  mutate(facility = recode(facility,
                           "Barraken Clinic" = "Facility A",
                           "Boniken Clinic" = "Facility B",
                           "Cavalla Clinic" = "Facility C",
                           "Cavalla Rubber Plantation Medical Center" = "Facility D",
                           "Edith H. Wallace Health Center" = "Facility E",
                           "Feloken Clinic" = "Facility F",
                           "Fish Town Clinic" = "Facility G",
                           "Gbarwiliken Clinic" = "Facility H",
                           "Gbloken Clinic" = "Facility I",
                           "Glofarken Clinic" = "Facility J",
                           "JJ Dossen Hospital" = "Facility K",
                           "Juduken Clinic (Barrobo Whojah District)" = "Facility L",
                           "Litter Wlebo Clinic" = "Facility M",
                           "Manolu Clinic" = "Facility N",
                           "Newaken Clinic" = "Facility O",
                           "Pleebo Health Center" = "Facility P",
                           "Pougbaken Clinic" = "Facility Q",
                           "Pullah Clinic" = "Facility R",
                           "Rock Town Clinic" = "Facility S",
                           "Rock Town Kunokudi Clinic" = "Facility T",
                           "Sacred Heart Clinic" = "Facility U",
                           "Sodoken Clinic" = "Facility V",
                           "St Francis Clinic" = "Facility W",
                           "Yediaken Clinic" = "Facility X",
                           "Youkudi Clinic" = "Facility Y"),
         district = recode(district,
                           "Barrobo Farjah District" = "District 1",
                           "Barrobo Whojah District" = "District 2",
                           "Harper District" = "District 3",
                           "Karluway 1 District" = "District 4",
                           "Karluway 2 District" = "District 5",
                           "Pleebo District" = "District 6"),
         county = recode(county,
                         "Maryland" = "County Alpha"))

         dummy <-  dummy %>% mutate(indicator_count_ari_total = round(runif(n(),indicator_count_ari_total/1.25,indicator_count_ari_total*1.25),0),
         indicator_count_ari_under5 = round(runif(n(),indicator_count_ari_under5/1.25,indicator_count_ari_total_proxy*1.25),0),
         indicator_count_ari_over5 = round(runif(n(),indicator_count_ari_over5/1.25,indicator_count_ari_under5*1.25),0),
         indicator_denom = round(runif(n(),indicator_denom/1.25,indicator_count_ari_over5*1.25),0),
         indicator_denom_under5 = round(runif(n(),indicator_denom_under5/1.25,indicator_count_ari_ipd_death_under5*1.25),0),
         indicator_denom_over5 = round(runif(n(),indicator_denom_over5/1.25,indicator_count_ari_ipd_death_over5*1.25),0))


# Export new dummy data as an RDS into new repo
data_example_singlecounty <- saveRDS(dummy, file = "../../global_covid19_response/data/data_example_singlecounty.rds")