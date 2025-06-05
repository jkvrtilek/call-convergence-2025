# adapting Grace's fundamental tracking code to run on OSC
# Julia Vrtilek
# 18 August 2023

# load packages
X <- c("plyr", "tidyverse", "data.table", "warbleR", "phonTools", "pbapply")
invisible(lapply(X, library, character.only = TRUE))

# load arguments from sbatch
args = commandArgs(trailingOnly=TRUE)

# load RDS
est <- readRDS(paste("/fs/ess/scratch/PAS1986/downExtSelTables", args, sep = "/"))

# for loop iterating over selections
for (i in 1:nrow(est)) {
  
  # make wave object
  r <- attr(est, "wave.objects")[[i]]
  samplingRate <- r@samp.rate
  windowLength <- (512/samplingRate) * 1000
  
  # get formant with soundgen
  sg <- soundgen::analyze(r, windowLength = windowLength, step = (windowLength/3), pitchFloor = 75, pitchCeiling = samplingRate/2, priorMean = NA, priorSD = NA, pitchMethods = c('autocor'), plot = TRUE, ylim = c(0, (samplingRate/2000)), interpol = NULL, pathfinding = "none", smooth = 0, nCands = 1, smoothVars = NULL, nFormants = 5, shortestSyl = 5, shortestPause = 20, formants = list(maxbw = 20000, coeffs = 20))
  formant_df <- sg$detailed
  
  saveRDS(formant_df, paste("/fs/ess/scratch/PAS1986/downfundfreq/ds_",est$sound.files[i],"_formant.RDS",sep=""))
  
  # Each row is a sound file, and each column is a value at each unique time point for ONE formant time series (cannot place multiple formants at the same time with warbleR::track_freq_contour)
  custom_contour_df <- formant_df %>%
    dplyr::mutate(
      sound.files = est$sound.files[i],
      selec = est$selec[i]
    ) %>%
    # Then select points within the start and end coordinates for the given call, otherwise the time series won't be plotted correctly below
    dplyr::filter(
      time >= (est$start[i])*1000 - windowLength & time <= (est$end[i])*1000 + windowLength
    ) %>%
    dplyr::rename(
      'f0' = "pitch"
    ) %>%
    dplyr::select(sound.files, selec, `f0`) %>%
    # Convert to kHz
    dplyr::mutate(
      `f0` = `f0`/1000
    ) %>%
    rowid_to_column() %>%
    dplyr::mutate(
      rowid = paste("f0-", rowid, sep = "")
    ) %>%
    # Make this data frame wider so that each row is a sound file
    pivot_wider(
      names_from = "rowid",
      values_from = `f0`
    )
  
  saveRDS(custom_contour_df, paste("/fs/ess/scratch/PAS1986/downfundfreq/ds_",est$sound.files[i],"_contour.RDS",sep=""))
  
}
