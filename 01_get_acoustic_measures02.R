# Contact call convergence in vampire bats
# combine 27 acoustic measures from warblr with the fundamental frequency summary measures
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)

# set working directory
setwd(dirname(file.choose()))

# get 27 acoustic measures from warblr
specan <- readRDS("acoustic data/ds_completespecan_2023-10-02.RDS")

# get other acoustic measures of fundamental frequency
ff <- readRDS("acoustic data/ds_fundfreq_summary_2023-11-06.RDS")

# combine datasets by sound.files
d <- left_join(specan, ff, by = "sound.files")

# fix typo in label
d$sound.files <- gsub("2018-08-29","2017-09-03", d$sound.files)

# convert seconds to milliseconds
d$duration <- d$duration * 1000
d$time.median <- d$time.median * 1000
d$time.Q25 <- d$time.Q25 * 1000
d$time.Q75 <- d$time.Q75 * 1000
d$time.IQR <- d$time.IQR * 1000

# data cleaning-----------

# filter duration, peak frequency, and time variables to remove sounds that are not contact calls
d2 <- 
  d %>% 
  filter(duration > 3) %>% 
  filter(duration < 50) %>% 
  filter(peakf > 10) %>% 
  filter(peakf < 30) %>% 
  filter(time.Q25 > 0) %>% 
  filter(time.median > 0) %>% 
  filter(time.Q75 > 0) %>% 
  filter(time.IQR > 0) %>% 
  # add binary variable for whether fundamental frequency measures succeeded
  mutate(indicator = case_when(!is.na(meanslope) ~ T,
                               is.na(meanslope) ~ F,
                               is.nan(meanslope) ~ F))

# save measures
write.csv(d2, "vampire_call_measures.csv", row.names= F)
