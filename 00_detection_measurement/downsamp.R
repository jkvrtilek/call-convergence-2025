# downsample all WAV files with a sampling rate above 250 kHz
# Julia Vrtilek
# Sep 25 2023, modified Jan 26 2024 for new WAV location in home directory

# load packages
library(tidyverse)
library(warbleR)

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# split args so I can use the pieces
x <- strsplit(args, split = "/")

# create directories in downsampled folder
dir.create(paste("/users/PAS1918/vrtilek/downsampled_calls",x[[1]][6],sep = "/"), showWarnings = FALSE)
dir.create(paste("/users/PAS1918/vrtilek/downsampled_calls",x[[1]][6],x[[1]][7],sep = "/"), showWarnings = FALSE)

# make list containing character vectors of wavs
wavs <- list.files(path = args[1], pattern = ".WAV$")

# downsampling loop, converts WAVs one at a time
for (i in 1:length(wavs)) {
  # read wave
  a <- readWave(paste(args,wavs[i],sep="/"))
  # if sampling rate is 250kHz, just copy the file
  # if > 250kHz, downsample to 250kHz
  # if < 250kHz, return error
  if (a@samp.rate == 250000) {
    a2 <- a
  } else if (a@samp.rate > 250000) {
    a2 <- resamp(a, g = 250000, output = "Wave")
  } else if (a@samp.rate < 250000) {
    write.table("x",paste("/users/PAS1918/vrtilek/downsampled_calls/",wavs[i],".txt",sep=""))
  }
  writeWave(a2, paste("/users/PAS1918/vrtilek/downsampled_calls/",x[[1]][6],"/",x[[1]][7],"/ds_",wavs[i],sep=""))
}




