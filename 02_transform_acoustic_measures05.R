# Contact call convergence in vampire bats
# transform acoustic measures to change or remove units
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)

# function to make scale() function return a vector not a matrix
scale2 <- function(x){as.vector(scale(x, scale = FALSE))}

# get data
raw <- read.csv("vampire_call_measures.csv")

# convert time variables into percentages 
d <- 
  raw %>% 
  mutate(time.Q25 = time.Q25/duration,
         time.median = time.median/duration,
         time.Q75 = time.Q75/duration,
         time.IQR = time.IQR/duration) %>% 
  # scale all numeric variables to remove units
  mutate(across(.cols=duration:peakf, .fns = scale2)) %>% 
  # exclude calls without fundamental frequency measures
  filter(indicator == TRUE)

# save data
write.csv(d, "vampire_call_measures_transformed.csv", row.names= F)


# other possible transformations of variables below--------  

# # how to transform variables
# wsqrt <- c("dfrange","abs_minslope","pos_slopes","neg_slopes")
# wlog <- c("duration","time.median","time.Q25","time.Q75","time.IQR","skew","kurt","sfm","mindom","enddom","meanpeakf","turns")
# winv <- c("freq.Q25","meandom","modindx")
# rsqrt <- c("freq.median")
# rlog <- c("freq.Q75","freq.IQR","segments")
# rinv <- c("sp.ent","time.ent","entropy")
# none <- c("meanfreq","sd","maxdom","startdom","dfslope","peakf","maxslope","minslope","meanslope")
# 
# # do transformations
# # note: k = max + 1
# 
# # square root
# for (i in 1:length(wsqrt)) {
#   d[wsqrt[i]] <- sqrt(raw[wsqrt[i]])
# }
# 
# # log
# for (i in 1:length(wlog)) {
#   d[wlog[i]] <- log(raw[wlog[i]] +
#                     min(raw[wlog[i]][raw[wlog[i]] > 0], na.rm = T))
# }
# 
# # inverse
# for (i in 1:length(winv)) {
#   d[winv[i]] <- 1/raw[winv[i]]
# }
# 
# # reflect and root
# for (i in 1:length(rsqrt)) {
#   d[rsqrt[i]] <- sqrt(max(raw[rsqrt[i]]) +
#                       min(raw[rsqrt[i]][raw[rsqrt[i]] > 0], na.rm = T) -
#                       raw[rsqrt[i]])
# }
# 
# # reflect and log
# for (i in 1:length(rlog)) {
#   d[rlog[i]] <- log(max(raw[rlog[i]],na.rm=T) +
#                     min(raw[rlog[i]][raw[rlog[i]] > 0], na.rm = T) -
#                     raw[rlog[i]])
# }
# 
# # reflect and inverse
# for (i in 1:length(rinv)) {
#   d[rinv[i]] <- 1/(max(raw[rinv[i]]) +
#                    min(raw[rinv[i]][raw[rinv[i]] > 0], na.rm = T) -
#                    raw[rinv[i]])
# }

# # # where indicator is FALSE, turn all fund freq measurements to 0
# # ff <- colnames(d2)[30:37]
# # 
# # for (i in 1:nrow(d2)) {
# #   if (d2[i,"indicator"] == FALSE) {
# #     for (j in 1:length(ff)) {
# #       d2[i,ff[j]] <- 0
# #     }
# #     print(i)
# #   }
# # }
