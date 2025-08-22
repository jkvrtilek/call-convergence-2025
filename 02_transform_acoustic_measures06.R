# Contact call convergence in vampire bats
# transform acoustic measures to change or remove units
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)

# function to make scale() function return a vector not a matrix
scale2 <- function(x){as.vector(scale(x, scale = FALSE))}

# Data for this script is available on Figshare: Vrtilek, Julia K.; Smith-Vidaurre, Grace; Carter, Gerald (2025). Data for "Vocal convergence during formation of cooperative relationships in vampire bats". figshare. Dataset. https://doi.org/10.6084/m9.figshare.29191334.v1

# get data
raw <- read.csv("vampire_call_measures.csv")

d <- 
  raw %>% 
  # convert time variables into percentages 
  mutate(time.Q25 = time.Q25/duration,
         time.median = time.median/duration,
         time.Q75 = time.Q75/duration,
         time.IQR = time.IQR/duration) %>% 
  # scale all numeric variables to remove units
  mutate(across(.cols=duration:peakf, .fns = scale2))

# save data
write.csv(d, "vampire_call_measures_transformed.csv", row.names= F)
