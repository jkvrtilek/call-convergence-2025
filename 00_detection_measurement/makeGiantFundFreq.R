# Julia Vrtilek
# GOAL: read all fundamental frequency summary measures into one dataframe
# 9 Oct 2023

# load packages
library("plyr")

# set wd
setwd("/fs/ess/scratch/PAS1986/down_ff_bybat")

data <- ldply(list.files(pattern = "results*"), readRDS)

saveRDS(data, "/fs/ess/PAS1986/ds_fundfreq_summary.RDS")
write.csv(data, "/fs/ess/PAS1986/ds_fundfreq_summary.csv")