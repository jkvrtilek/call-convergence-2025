# Julia Vrtilek
# GOAL: make an extended selection table for each file path
# POST downsampling and POST move to home directory, so no individual automatic_detections files; use whole selection table
# 28 Jan 2024

# load packages
X <- c("dplyr", "tidyverse", "data.table", "warbleR", "tuneR", "pbapply")
invisible(lapply(X, library, character.only = TRUE))

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# get bat name from args
temp <- strsplit(args[1],"/")
batname <- temp[[1]][6]
thisdate <- as.Date(temp[[1]][7])

# load whole selection table
table <- readRDS("/fs/ess/PAS1986/data/selectionTables_08-2023.RDS") %>% 
  separate(sound.files, into = c("date","bat","WAV"), sep = "_", remove = FALSE)

# keep only WAVs for target bat and filter out NA selections
x <- table %>% 
  filter(bat == batname) %>% 
  filter(date == thisdate | date == thisdate+1) %>% 
  filter(!is.na(selection_length)) %>% 
  filter(selection_length >= 0)

# make column listing downsampled WAVs
x2 <- x %>% 
  mutate(sound.files = paste("ds",sound.files,sep="_"))

# run function
st <- selection_table(x2, path = args, extended = TRUE, mar = 0.05)

# get name for ext sel table
name <- paste("ds",temp[[1]][7],temp[[1]][6],sep="_")

#save ext sel table
saveRDS(st, file = paste("/users/PAS1918/vrtilek/downExtSelTables/", name, ".RDS", sep = ""))
