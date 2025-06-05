# get summary measures from fundamental frequency contours:
# number of segments
# Y max slope across segments
# Y min slope across segments
# min (abs(slope))
# number of upward positive slopes
# number of downward negative slopes
# Y number of turns
# Y mean slope

# Julia Vrtilek
# Sep-Oct 2023

# load packages
library(tidyverse)

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# load data
raw <- readRDS(paste("/fs/scratch/PAS1986/down_ff_bybat/", args, sep = "")) %>% 
  select(-selec)

# remove all-NA rows
rawNoNA <- raw[rowSums(is.na(raw)) != ncol(raw) - 2, ]

# tidy data
ff <- rawNoNA %>% 
  pivot_longer(cols = 2:ncol(raw), names_to = "fundnum", values_to = "contour")

ff$sound.files <- as.factor(ff$sound.files)

d <- split(ff, f = ff$sound.files)

# create df to hold results
results <- ff %>% 
  select(sound.files) %>% 
  distinct() %>% 
  mutate(maxslope = NA) %>% 
  mutate(minslope = NA) %>% 
  mutate(abs_minslope = NA) %>%  
  mutate(pos_slopes = NA) %>% 
  mutate(neg_slopes = NA) %>% 
  mutate(turns = NA) %>% 
  mutate(meanslope = NA) %>% 
  mutate(segments = NA)

results[2:9] <- as.numeric(results$meanslope)

# function to count sign changes
nSignChanges <- function(x) {
  signs <- sign(x)
  sum(signs[-1] != signs[-length(x)])
}

# collect all slope-related measures
for (i in 1:length(d)) {
  d[[i]] <- d[[i]] %>% 
    mutate(slope = contour-lag(contour))
  slope <- d[[i]]$slope %>% 
    na.omit()
  results[i,"maxslope"] <- max(slope, na.rm = TRUE)
  results[i,"minslope"] <- min(slope, na.rm = TRUE)
  results[i,"abs_minslope"] <- min(abs(slope), na.rm = TRUE)
  results[i,"pos_slopes"] <- length(slope[slope >= 0])
  results[i,"neg_slopes"] <- length(slope[slope < 0])
  results[i,"turns"] <- nSignChanges(slope)
  results[i,"meanslope"] <- mean(slope, na.rm = TRUE)
}

# collect number of segments
for (i in 1:length(d)) {
  if (sum(!is.na(d[[i]]$contour)) <= 1) {
    results[i,"segments"] <- 0
  } else {
  fundnum <- d[[i]] %>% 
    filter(!is.na(contour)) %>% 
    separate(fundnum, into = c("trash","num")) %>% 
    mutate(num = as.integer(num)) %>% 
    mutate(test = num - lag(num))
  segments <- sum(fundnum$test > 1, na.rm = TRUE)
  results[i,"segments"] <- segments + 1
  }
}


saveRDS(results, file = paste("/fs/scratch/PAS1986/down_ff_bybat/results_", args, sep = ""))

