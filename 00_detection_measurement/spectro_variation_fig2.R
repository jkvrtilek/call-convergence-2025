# get 5 spectrograms per session from 4 bats that were recorded 4+ times
# Julia Vrtilek
# assembled from earlier versions on Mar 1 2025

library(tidyverse)
library(warbleR)

# get list of bats with 4+ recording sessions ----

# load data
raw <- readRDS("/fs/ess/PAS1986/data/selectionTables_08-2023.RDS")

# load list of usable calls for filter
filt <- read.csv("/fs/ess/PAS1986/data/vampire_call_measures.csv")

st <- raw %>% 
  mutate(sel.filt = paste("ds", sound.files, selec, sep = "_")) %>% 
  filter(sel.filt %in% filt$sound.files) %>% 
  separate(sound.files, into=c('date', 'bat', 'file'), sep="_", remove = FALSE)

# get list of dates bats were recorded on
dates <- st %>% 
  dplyr::select(date,bat) %>% 
  distinct() %>% 
  arrange(bat)

# add "session" to catch bats that were recorded over midnight
dates$date <- as.Date(dates$date)
dates <- dates %>% 
  mutate(midnight = ifelse(bat != lag(bat), FALSE,
                           ifelse(lag(date) == date-1, TRUE,
                                  FALSE))) %>% 
  mutate(session = ifelse(midnight == TRUE, date - 1,
                          ifelse(midnight == FALSE, date,
                                 ifelse(is.na(midnight), date,
                                        NA)))) %>% 
  mutate(id = paste(bat, date, sep = "/")) %>% 
  select(!midnight)
dates$session <- as.Date(dates$session, origin="1970-01-01")

# fix first session row
dates$session[1] <- dates$session[2]

# fix veronica 2013-02-15 (seems like no calls were recorded before midnight so above code fails)
dates <- dates %>% 
  mutate(session = case_when(id == "veronica/2013-02-15" ~ as.Date("2013-02-14"),
                             TRUE ~ session))

# filter for bats with 4 or more sessions
session <- dates %>% 
  select(bat, session) %>% 
  distinct() %>% 
  group_by(bat) %>%
  filter(n()>4) %>% 
  filter(!is.na(session)) %>% 
  separate(session, into = c("year","month","day"), remove = FALSE) %>% 
  mutate(year = as.integer(year)) %>% 
  ungroup()

# get list of bats with 4 or more sessions
bats <- unique(session$bat)

# choose 4 bats: 1 from pre-2016, 1 from 2016, 2 from 2019 ----
set.seed(123)

pre.2016 <- session %>% 
  filter(year < 2016) %>% 
  slice_sample(n = 1)

postdoc <- session %>% 
  filter(year > 2015) %>% 
  filter(year < 2018) %>% 
  filter(bat != "scs") %>% # present in 2016-17 and 2019
  filter(bat != "ss") %>% # too few calls
  slice_sample(n = 1)

prof <- session %>% 
  filter(year == 2019) %>% 
  filter(bat != "scs") %>% # present in 2016-17 and 2019
  slice_sample(n = 2)

samp.bats <- c(pre.2016$bat, postdoc$bat, prof$bat)

# separate into a list with one df per bat session ----
d <- st %>% 
  filter(bat %in% samp.bats) %>% 
  mutate(sound.files = paste("ds", sound.files, sep = "_")) %>% 
  mutate(id = paste(bat, date, sep = "/")) %>% 
  select(!date:file) %>% 
  left_join(dates, by = "id")

# pull session data and pick 4 sessions per bat
bat_date <- session %>% 
  filter(bat %in% samp.bats) %>% 
  group_by(bat) %>% 
  slice_sample(n = 4)

# filter selection table by 4 chosen sessions
st_by_session <- d %>% 
  filter(session %in% bat_date$session)

# separate into a list with one df per session
sess_list <- split(st_by_session, f = st_by_session$session)

# choose 5 random files from each df in the list
sels <- data.frame(matrix(ncol = 11))
colnames(sels) <- colnames(st_by_session)

set.seed(125)

for (j in 1:length(sess_list)) {
  y <- sample(nrow(sess_list[[j]]), size = 3)
  sels <- rbind(sels, sess_list[[j]][y[1],], sess_list[[j]][y[2],], sess_list[[j]][y[3],])
}

# make duration of all spectrograms match
sels2 <- sels[2:nrow(sels),] %>% 
  mutate(duration = end - start)

longest <- max(sels2$duration)

sels3 <- sels2 %>% 
  mutate(adjustment = longest - duration) %>% 
  mutate(old.end = end) %>% 
  mutate(end = old.end + adjustment) %>% 
  mutate(WAVpath = "a")

# get path to WAV file for each selection
for (k in 1:nrow(sels2)) {
  sels3$WAVpath[k] <- paste("/users/PAS1918/vrtilek/downsampled_calls/",
                            sels3$bat[k], "/",
                            as.Date(sels3$session[k], origin="1970-01-01"),
                            sep = "")
}

# make spectrogram for each selection
for (l in 1:nrow(sels3)) {
  spectrograms(sels3[l,], flim = c(0,120), ovlp = 50, mar = 0.005, res = 400,
               line = FALSE, title.labels = NULL, gr = FALSE, osci = TRUE,
               path = sels3[l,15], 
               dest.path = "/fs/ess/PAS1986/randspecfig")
}


