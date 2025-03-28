# Contact call convergence in vampire bats
# combine vocal distances and social data
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# set working directory
#setwd("/Users/gc1511/Dropbox (Personal)/Dropbox/_working/_ACTIVE/__students/Julia Vrtilek/2024/call convergence/2024-09") 

# load packages
library(tidyverse)

# function
matrix_to_df <- function(m1){
  data.frame(dyad = paste(rownames(m1)[col(m1)], colnames(m1)[row(m1)], sep="_"),
             value = c(t(m1)), stringsAsFactors = FALSE)
}

# get vocal distances
dist <- as.matrix(read_csv("vocal-distance-94-bats.csv"))
rownames(dist) <- colnames(dist)
distdf <- 
  matrix_to_df(dist) %>% 
  rename(dist= value) %>% 
  separate(dyad, into = c("bat1","bat2"), sep = "_", remove = F) %>% 
  filter(bat1 != bat2) %>% 
  mutate(dyad = ifelse(bat1<bat2, paste(bat1,bat2, sep="_"), paste(bat2,bat1, sep="_"))) %>% 
  group_by(dyad) %>% 
  summarize (dist= mean(dist, na.rm=T))

# get bat sexes and ages at recording
sex <- 
  read.csv("sex_age02.csv") %>% 
  mutate(birthdate = mdy(birthdate), date= mdy(date)) %>% 
  mutate(age.days = as.numeric(date- birthdate)) %>% 
  mutate(age.years = age.days/365) %>% 
  group_by(bat, sex) %>% 
  summarize(age= mean(age.years),
            n.recorded.dates=n(), .groups= 'drop') %>% 
  mutate(adult= is.na(age) | age>2) 

# get kinship
k <- 
  read_csv("kinship-all-bats01.csv") %>% 
  group_by(dyad) %>% 
  arrange(dyad, bat1) %>% 
  summarize(bat1= first(bat1), bat2= first(bat2), site1= first(site1), site2= first(site2), colony1= first(colony1), colony2= first(colony2), kinship= first(kinship)) %>% 
  ungroup()

# get directed social data
directed.social <- read_csv("2010-2019_vamp_social_interactions01.csv") 

# get undirected social data
social <- 
  directed.social %>% 
  mutate(dyad = ifelse(actor<receiver, paste(actor,receiver, sep= "_"), paste(receiver,actor, sep= "_"))) %>% 
  group_by(dyad) %>% 
  summarize(study= first(study),
            donation.rate = mean(donation.rate, na.rm=T),
            donation.lograte = mean(donation.lograte, na.rm=T),
            groom.rate = mean(groom.rate, na.rm=T),
            groom.lograte = mean(groom.lograte, na.rm=T),
            contact.rate = mean(contact.rate, na.rm=T),
            contact.lograte = mean(contact.lograte, na.rm=T)) 

# combine all data
d <- 
  full_join(k,social,by = "dyad") %>% 
  full_join(distdf,by = "dyad") %>% 
  # get bats with vocal distance
  filter(!is.na(dist)) %>% 
  # add sex 
  mutate(sex1= sex$sex[match(.$bat1, sex$bat)]) %>% 
  mutate(sex2= sex$sex[match(.$bat2, sex$bat)]) %>% 
  mutate(sex1= ifelse(bat1== "war-machine", "male", sex1)) %>% 
  mutate(sex2= ifelse(bat2== "war-machine", "male", sex2)) %>% 
  # label adults 
  mutate(adult1= sex$adult[match(.$bat1, sex$bat)]) %>% 
  mutate(adult2= sex$adult[match(.$bat2, sex$bat)]) %>% 
  mutate(adult1= ifelse(bat1== "war-machine", TRUE, adult1)) %>% 
  mutate(adult2= ifelse(bat2== "war-machine", TRUE, adult2)) %>% 
  mutate(both_adult = adult1 & adult2) %>% 
  mutate(both_young = !adult1 & !adult2) %>% 
  # label dyad for capture site
  mutate(site= case_when(
    site1 == "zoo" & site2 == "zoo" ~ "US",
    site1 == "chorrera" & site2 == "chorrera" ~ "CH",
    site1 == "tole" & site2 == "tole" ~ "TL",
    site1 == "las.pavas" & site2 == "las.pavas" ~ "LP",
    site1 == "lake.bayano" & site2 == "lake.bayano" ~ "LB",
    site1 == "chilibre" & site2 == "chilibre" ~ "CB",
    site1 != site2  ~ "different sites",
    TRUE ~ "other")) %>% 
  # label dyad by lab colony
  mutate(study= case_when(
    colony1 == "2014-colony" & colony2 == "2014-colony" ~ "2014",
    colony1 == "2017-colony" & colony2 == "2017-colony" ~ "2017",
    colony1 == "2019-colony" & colony2 == "2019-colony" ~ "2019",
    colony1 != colony2 & is.na(donation.rate)  ~ "different years",
    # bats d and scs were in both 2017 and 2019 datasets (but are labeled as 2019 bats)
    # so this line is necessary
    colony1 != colony2 & !is.na(donation.rate)  ~ "2017",
    TRUE ~ "other")) %>% 
  # add zoo vs wild
  mutate(both_zoo = site1 == "zoo" & site2 == "zoo") %>% 
  mutate(both_wild = site1 != "zoo" & site2 != "zoo") %>% 
  # label social experience treatments
  mutate(study_site= paste(study, site, sep="_")) %>% 
  # rank treatments 
  # 1 = never met
  # 2 = caught in same roost two years apart
  # 3 = introduced in captivity
  # 4 = caught together and caged together
  # 5 = long-term captive colony
  mutate(treatment_rank = case_when(
    study_site == "2019_CH" ~ 4,
    study_site == "2019_LB" ~ 4,
    study_site == "2017_LP" ~ 4,
    study_site == "2019_TL" ~ 4,
    study_site == "2017_TL" ~ 4,
    study_site == "2017_CB" ~ 4,
    study_site == "2014_US" ~ 5,
    study_site == "2017_different sites" ~ 3,
    study_site == "2019_different sites" ~ 3,
    study_site == "different years_TL" ~ 2,
    study_site == "different years_different sites" ~ 1,
    TRUE ~ 0)) %>% 
  mutate(treatment_label = case_when(
    treatment_rank == 1 ~ "1. never met",
    treatment_rank == 2 ~ "2. same wild roost\nbut different years",
    treatment_rank == 3 ~ "3. different wild roost\nthen caged together",
    treatment_rank == 4 ~ "4. same wild roost\nthen caged together",
    treatment_rank == 5 ~ "5. same long-term \ncaptive colony",
    TRUE ~ "other")) %>% 
  # add vocal similarity 
  # sim = inverse of distance on scale of 0 to 1
  mutate(sim = 1- (dist/max(dist))) %>% 
# label sexes of dyad
  mutate(dyad.sex = case_when(
    sex1 == "female" & sex2 == "female" ~ "female",
    sex1 == "male" & sex2 == "male" ~ "male",
    TRUE ~ "mixed_sex")) %>% 
  # add kinship categories
  mutate(kinship2= case_when(
    kinship == 0 ~ "0",
    kinship < 0.25 ~ "0.01-0.249",
    kinship >= 0.25 ~"0.25-0.5", 
    TRUE ~ NA_character_)) %>% 
  # beta models require response to be greater than zero
  # we have one observation with similarity = 0 (next lowest value is 0.111)
  # change that zero to 0.001
  mutate(sim = ifelse(sim ==0, 0.001, sim))

# get descriptive statistics for bats included in analysis
# check bats
(bats <- 
  tibble(bat= c(d$bat1, d$bat2), 
         colony= c(d$colony1, d$colony2), 
         sex= c(d$sex1, d$sex2), 
         adult= c(d$adult1, d$adult2), 
         site= c(d$site1, d$site2)) %>% 
  group_by(bat) %>% 
  summarize(colony= first(colony), sex= first(sex), adult= first(adult), site= first(site)))

# check number of bats by populations
bats %>% 
  mutate(colony= paste0(colony, site)) %>% 
  group_by(colony,sex, adult) %>% 
  summarize(n=n())

# check number of pairs by familiarity
d %>% 
  group_by(treatment_label) %>% 
  summarize(n=n())

# check pairs with known food sharing rates, allogrooming rate, or both
d %>% 
  filter(!is.na(donation.rate)) %>% 
  filter(is.na(groom.rate))
  
d %>% 
  filter(is.na(donation.rate)) %>% 
  filter(!is.na(groom.rate))

d %>% 
  filter(!is.na(donation.rate)) %>% 
  filter(!is.na(groom.rate))

# select variables
d2 <- 
  d %>% 
  dplyr::select(dyad, dyad.sex, bat1, bat2, both_zoo, both_wild, both_adult, both_young, study_site, treatment_rank, treatment_label, kinship, kinship2, donation.rate, donation.lograte, groom.rate, groom.lograte, contact.rate, contact.lograte, dist, sim)

# save data
write.csv(d2,"vocal_social_data.csv", row.names=F)

