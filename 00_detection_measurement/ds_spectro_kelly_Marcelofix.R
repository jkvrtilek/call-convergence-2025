# for ds_2019-08-11_kelly.RDS, which didn't fail, it just got stuck running
# load packages
library(tidyverse)
library(warbleR)

args <- "ds_2019-08-11_kelly.RDS"

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

# run over subsets of the data, manually determined by when progress bar fails
specan1 <- spectro_analysis(extSelTable[1:20,], path = "/users/PAS1918/vrtilek/downExtSelTables",
                            bp = c(5,100),
                            wl = 512,
                            threshold = 10,
                            fast = FALSE,
                            ovlp = 50,
                            wn = "blackman")

specan2 <- spectro_analysis(extSelTable[25:190,], path = "/users/PAS1918/vrtilek/downExtSelTables",
                            bp = c(5,100),
                            wl = 512,
                            threshold = 10,
                            fast = FALSE,
                            ovlp = 50,
                            wn = "blackman")

specan3 <- spectro_analysis(extSelTable[195:225,], path = "/users/PAS1918/vrtilek/downExtSelTables",
                            bp = c(5,100),
                            wl = 512,
                            threshold = 10,
                            fast = FALSE,
                            ovlp = 50,
                            wn = "blackman")

specan4 <- spectro_analysis(extSelTable[230:240,], path = "/users/PAS1918/vrtilek/downExtSelTables",
                            bp = c(5,100),
                            wl = 512,
                            threshold = 10,
                            fast = FALSE,
                            ovlp = 50,
                            wn = "blackman")

specan5 <- spectro_analysis(extSelTable[250:391,], path = "/users/PAS1918/vrtilek/downExtSelTables",
                            bp = c(5,100),
                            wl = 512,
                            threshold = 10,
                            fast = FALSE,
                            ovlp = 50,
                            wn = "blackman")

# combine all data
kelly_specan <- rbind(specan1, specan2, specan3, specan4, specan5)

#check for duplicates
sum(duplicated(kelly_specan))

#save
saveRDS(kelly_specan, file = paste("/users/PAS1918/vrtilek/fixed_Marcelo/spectro_analysis_", args[1], sep = ""))
