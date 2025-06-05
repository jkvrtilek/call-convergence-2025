# Julia Vrtilek
# GOAL: run energy detector on wav files; spit out a CSV for each recording
# 27 Apr 2021
# Modified for OSC 5 Aug 2021

# load packages
X <- c("tidyverse", "data.table", "warbleR", "tuneR", "pbapply")
invisible(lapply(X, library, character.only = TRUE))

# source for function
source(file.path("/fs/ess/PAS1986/", "energyDetectorOSC.R"))

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

#make list containing character vectors of wavs
wavs <- list.files(path = args[1], pattern = ".WAV$")

#designate start date
startDate <- Sys.Date()

#make error log for path
file.create(paste("energyDetectorLog_", startDate, ".txt", sep = ""))

# run function
for (i in 1:length(wavs)) {
  energyDetector(wavs = wavs[i], path = args[1], outpath1 = "/fs/ess/scratch/PAS1986/selectionTables", outpath2 = args[1], lowerfreq = 40000, higherfreq = 100000, energyoffset = 2.5, amplitudethreshold = 0.01, spec_fft_size = 512, spec_fft_shift = 4, spec_energy_floor = -120, truncate = FALSE, shouldplot = FALSE, strategy = "hybrid", boundarystrategy = "amplitude", date = startDate)
}