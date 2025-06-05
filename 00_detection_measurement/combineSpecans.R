# Julia Vrtilek
# combine giant and small specan tables
# 3 Sept 2023

library(tidyverse)

giant <- readRDS("/fs/ess/PAS1986/data/ds_giantspecan.RDS")

small <- readRDS("/fs/ess/PAS1986/data/ds_smallspecan.RDS")

all <- rbind(giant,small)

saveRDS(all, "/fs/ess/PAS1986/ds_completeSpecan.RDS")
write.csv(all, "/fs/ess/PAS1986/ds_completeSpecan.csv")
