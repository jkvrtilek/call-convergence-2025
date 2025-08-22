# load packages
library(tidyverse)
library(warbleR)

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# read in data
extSelTable <- readRDS(paste("/users/PAS1918/vrtilek/downExtSelTables/", args[1], sep = ""))

# remove clipped selections
extSelTable <- extSelTable %>% 
  filter(clipped == "N")

# create wav.size variable for filtering
attrcheck <- attr(extSelTable, "check.results") %>% 
  select(sound.files, wav.size)

# filter out WAVs that are too small/short
extSelTable <- extSelTable %>% 
  left_join(attrcheck) %>% 
  filter(wav.size > 0.000001) %>% 
  filter(selection_length > 0.0005)

# Arguments
# X = one extended selection table; looping is in sbatch command
# bp = 5 to 100 kHz
# wl = 512 (default)
# threshold = 10, from Gerry's 2012 PLOS One paper
# fast = FALSE
# path = directory containing all extended selection tables
# ovlp = 50 (default), from Gerry's 2012 PLOS One paper
# wn = "blackman", from Gerry's 2012 PLOS One paper
specan <- spectro_analysis(extSelTable, path = "/users/PAS1918/vrtilek/downExtSelTables",
                           bp = c(5,100),
                           wl = 512,
                           threshold = 10,
                           fast = FALSE,
                           ovlp = 50,
                           wn = "blackman")

saveRDS(specan, file = paste("/users/PAS1918/vrtilek/fixed_Marcelo/spectro_analysis_", args[1], sep = ""))
