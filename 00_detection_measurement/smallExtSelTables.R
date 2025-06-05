# Julia Vrtilek
# GOAL: make several extended selection tables with 100 selections each for the file paths that failed
# updated for that purpose 30 Aug 2023

# load packages
X <- c("plyr", "tidyverse", "data.table", "warbleR", "tuneR", "pbapply")
invisible(lapply(X, library, character.only = TRUE))

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# make list of selection tables in file path
x <- ldply(list.files(path = args[1], pattern = "^automatic_detections*", full.names = TRUE), read.csv, header = TRUE)

# filter out NA selections
x <- x %>% 
  filter(!is.na(selection_length)) %>% 
  filter(selection_length >= 0)

# make column listing downsampled WAVs
x <- x %>% 
  mutate(sound.files = paste("ds",sound.files,sep="_"))

temp <- strsplit(args[1],"/")
args[2] <- paste(temp[[1]][1],temp[[1]][2],temp[[1]][3],temp[[1]][4],"downsampled",temp[[1]][5],temp[[1]][6], sep = "/")

# get number of rounds the for loop will need to run
num <- ceiling(nrow(x)/100)
# initialize variable that IDs where each extended selection table should start
start <- 1

for (i in 1:num) {
  # set where extended selection table should end
  if (i < num) {
    end <- i*100
  } else {
    end <- nrow(x)
  }
  
  # select rows from x
  temp <- x[start:end,]
  
  # run function
  st <- selection_table(temp, path = args[2], extended = TRUE, mar = 0.05, confirm.extended = FALSE)
  
  # get name for ext sel table
  name <- str_split(st[[1]][1], "_")
  name <- paste(name[[1]][1], name[[1]][2], name[[1]][3], i, sep = "_")
  
  # save ext sel table
  saveRDS(st, file = paste("/fs/ess/scratch/PAS1986/downExtSelTables/small/", name, ".RDS", sep = ""))
  
  # index where extseltab should start
  start <- start + 100
}

