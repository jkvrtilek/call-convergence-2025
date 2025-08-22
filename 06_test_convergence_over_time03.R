# Tests of call convergence over time
## figure showing recording dates
## linear discriminant function analysis to get vocal distances
## permuted DFA to get classification accuracy
# Julia Vrtilek, Gerry Carter

library(MASS)
library(tidyverse)
library(graph4lg)
library(vegan)

setwd(dirname(file.choose()))

# load and tidy data ----

#* call/recording data - large file, download from Figshare ----
raw <- read.csv("vampire_call_measures_transformed.csv")

# tidy data, filter for 2016-17 bats, and add sample sizes
d <- 
  raw %>% 
  separate(sound.files, into=c('ds','date', 'bat', 'file', 'call'), sep="_", remove = FALSE) %>% 
  separate(date, into=c("year","month","day"), sep = "-", remove = FALSE) %>% 
  filter(year == "2016" | year == "2017") %>% 
  group_by(bat) %>% 
  mutate(sample.size= n()) %>%
  ungroup()

min.sample.size <- 80

#* capture sites ----
capture.site <- read.csv("sex_age03.csv") %>% 
  dplyr::select(bat, capture.site) %>% 
  distinct

# filter calls by min sample size and capture site
d2 <- d %>% 
  filter(sample.size >= min.sample.size) %>% 
  left_join(capture.site, by = "bat") %>% 
  filter(capture.site == "tole" | capture.site == "las.pavas")

# get list of relevant bats
bats <- d2 %>%
  filter(date > "2016-09-01") %>% 
  dplyr::select(bat, capture.site) %>% 
  distinct()

#* first contact data from 2016 ----
fc16 <- read_csv("first_contact_events2016.csv")

# make df of 2016 first contact dates
cont_date16 <- fc16 %>% 
  mutate(time = as.character(time)) %>% 
  separate(col = time, into = c("date", "time"), sep = " ") %>% 
  distinct(date) %>% 
  mutate(date = as.Date(date, "%Y-%m-%d"))
# Feb 7 is a typo for Jul 2
# Jul 7 and Aug 26 carry over from introductions on Jul 6 and Aug 25
# legitimate dates: 2016-07-02, 2016-07-06, 2016-08-25
smallcage <- cont_date16$date[1]
largecage <- cont_date16$date[5]


# make recording dates figure ----
# prep data frame
rec_dates <- d2 %>% 
  dplyr::select(sound.files, date, bat, capture.site) %>% 
  filter(bat %in% bats$bat) %>% 
  mutate(date = as.Date(date, "%Y-%m-%d")) %>% 
  mutate(bat = factor(bat, levels = c("eve","ola","ivy","cat","six","dcd","cs","ccs","scs","rc","ccc","r"))) %>% 
  mutate(CaptureSite = case_when(capture.site == "las.pavas" ~ "Las Pavas",
                                 capture.site == "tole" ~ "Tole"))

# make labels
labels_lp <- rec_dates %>% 
  group_by(bat) %>%
  arrange(desc(date)) %>%
  slice(1L) %>% 
  filter(capture.site == "las.pavas")

labels_tole <- rec_dates %>% 
  group_by(bat) %>%
  arrange(desc(date)) %>%
  slice(1L) %>% 
  filter(capture.site == "tole")

# make plot
p <- rec_dates %>% 
  ggplot(aes(x=date, y=bat, group=bat)) +
  facet_wrap(~CaptureSite, scales = "free_y", ncol = 1) +
  geom_point(aes(color = capture.site), size = 3) + 
  geom_line() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  scale_color_manual(values = c("#0f2d59","#3eb7c7")) +
  guides(color=guide_legend("Bat Origin")) +
  geom_vline(xintercept = as.numeric(smallcage), color = "red", linetype = "dashed") +
  geom_vline(xintercept = as.numeric(largecage), color = "red") +
  ylab("individual bats") +
  xlab("recording dates") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 16)) +
  geom_text(data = labels_lp, aes(label = bat), nudge_x = +18) +
  geom_text(data = labels_tole, aes(label = bat), nudge_x = +18)

# save plot
ggsave(
  "06_test_convergence_over_time/recording_dates.pdf",
  plot = p,
  scale = 1,
  width = 8,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


# compare dyadic distances between two groups of bats pre- and post-introduction ----

#* make functions ----
# get mean and 95% CI of values x via bootstrapping
boot_ci <- function(x, perms=5000, bca=F) {
  get_mean <- function(x, d) {
    return(mean(x[d]))
  }
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  if(bca){
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="bca")
    low <- boot$bca[1,4]
    high <- boot$bca[1,5] 
  }else{
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="perc")
    low <- boot$perc[1,4]
    high <- boot$perc[1,5] 
  }
  c(low=low,mean=mean,high=high, N=round(length(x)))
}

# function to get mean and 95% CI via bootstrapping of values y within grouping variable x
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000, bca=F){
  df <- data.frame(effect=unique(x))
  df$low <- NA
  df$mean <- NA
  df$high <- NA
  df$n.obs <- NA
  for (i in 1:nrow(df)) {
    ys <- y[which(x==df$effect[i])]
    if (length(ys)>1 & var(ys)>0 ){
      b <- boot_ci(y[which(x==df$effect[i])], perms=perms, bca=bca) 
      df$low[i] <- b[1]
      df$mean[i] <- b[2]
      df$high[i] <- b[3]
      df$n.obs[i] <- b[4]
    }else{
      df$low[i] <- min(ys)
      df$mean[i] <- mean(ys)
      df$high[i] <- max(ys)
      df$n.obs[i] <- length(ys)
    }
  }
  df
}

#* dyadic DFA between Las Pavas bats and pre-intro Tole bats ----
# prep dataframe
d3 <- d2 %>% 
  filter(bat %in% bats$bat) %>% 
  mutate(
    group = case_when(capture.site == "las.pavas" ~ "lp",
                      capture.site == "tole" & date < "2016-08-01" ~ "tole1",
                      capture.site == "tole" & date > "2016-08-01" ~ "tole2")
  )

# make dataframe for each group
# tole pre-intro bats
tole_pre <- d3 %>% 
  filter(group == "tole1") %>% 
  mutate(individual = bat) %>% 
  dplyr::select(group, individual, duration:segments)
# tole post-intro bats
tole_post <- d3 %>% 
  filter(group == "tole2") %>% 
  mutate(individual = bat) %>% 
  dplyr::select(group, individual, duration:segments)
# las pavas bats
lp <- d3 %>% 
  filter(group == "lp") %>% 
  mutate(individual = bat) %>% 
  dplyr::select(group, individual, duration:segments)

# get all PRE-INTRO dyadic distances
pre <- rbind(tole_pre,lp)

# classify calls to bat using a single dfa without cross validation
dfa_pre <- lda(individual ~ 
                 duration+    
                 meanfreq+    
                 sd+          
                 freq.median+ 
                 freq.Q25+ 
                 freq.Q75+    
                 freq.IQR+    
                 time.median+
                 time.Q25+   
                 time.Q75+    
                 time.IQR+    
                 skew+        
                 kurt+
                 sp.ent+      
                 time.ent+  
                 entropy+     
                 sfm+         
                 meandom+     
                 mindom+     
                 maxdom+      
                 dfrange+    
                 modindx+     
                 startdom+    
                 enddom+      
                 dfslope+   
                 meanpeakf+   
                 peakf+
                 maxslope+
                 minslope+
                 abs_minslope+
                 pos_slopes+
                 neg_slopes+
                 turns+
                 meanslope+
                 segments,
               CV= F, 
               data=pre)

# get vocal distance as Mahalanobis distance between group centroids
pre_distance <- as.matrix(dist(dfa_pre$means %*% dfa_pre$scaling))

# get all POST-INTRO dyadic distances
post <- rbind(tole_post,lp)

# classify calls to bat using a single dfa without cross validation
dfa_post <- lda(individual ~ 
                  duration+    
                  meanfreq+    
                  sd+          
                  freq.median+ 
                  freq.Q25+ 
                  freq.Q75+    
                  freq.IQR+    
                  time.median+
                  time.Q25+   
                  time.Q75+    
                  time.IQR+    
                  skew+        
                  kurt+
                  sp.ent+      
                  time.ent+  
                  entropy+     
                  sfm+         
                  meandom+     
                  mindom+     
                  maxdom+      
                  dfrange+    
                  modindx+     
                  startdom+    
                  enddom+      
                  dfslope+   
                  meanpeakf+   
                  peakf+
                  maxslope+
                  minslope+
                  abs_minslope+
                  pos_slopes+
                  neg_slopes+
                  turns+
                  meanslope+
                  segments,
                CV= F, 
                data=post)

#* make matrix containing change in distance ----

# get vocal distance as Mahalanobis distance between group centroids
post_distance <- as.matrix(dist(dfa_post$means %*% dfa_post$scaling))

# reorder matrices
order <- c(unique(tole_pre$individual), unique(lp$individual))
pre_dist <- reorder_mat(pre_distance, order)
post_dist <- reorder_mat(post_distance, order)

# get call similarity before introduction
sim1 <- 1 - (pre_dist/max(pre_dist))

# get call similarity after introduction
sim2 <- 1 - (post_dist/max(post_dist))

# get increase in similarity
diff <- sim2 - sim1

#* make familiar vs introduced matrix ----
# 0 = originally familiar (same capture site)
# 1 = introduced (different capture site)
lp_bats <- unique(lp$individual)
tole_bats <- unique(tole_pre$individual)

a <- matrix(nrow = 12, ncol = 12)

col1 <- c(rep(0,times=7),rep(1,times=5))
col2 <- c(rep(1,times=7),rep(0,times=5))
v <- c(rep(col1,times=7),rep(col2,times=5))
a[,] <- v

diag(a) <- NA

rownames(a) <- c(tole_bats,lp_bats)
colnames(a) <- c(tole_bats,lp_bats)

# check names match
colnames(a) == colnames(diff)

#* Mantel test ----
set.seed(123)
mantel(diff, a, na.rm=T, method = "spearman")
# r = 0.28, p = 0.007

# compared to familiar pairs, pairs that were introduced have greater increases in similarity


# get within-dyad convergence ----
dd <-  
  tibble(bat1=rownames(diff)[row(diff)], 
         bat2=colnames(diff)[col(diff)], 
         diff=c(diff)) %>% 
  filter(bat1 != bat2) %>% 
  mutate(site1= ifelse(bat1 %in% lp_bats, "lp", "tole")) %>% 
  mutate(site2= ifelse(bat2 %in% lp_bats, "lp", "tole")) %>% 
  mutate(dyad= ifelse(bat1<bat2, paste(bat1,bat2, sep="_"), paste(bat2,bat1, sep="_"))) %>% 
  mutate(introduced= site1!= site2) %>% 
  group_by(dyad, introduced) %>% 
  summarize(diff= mean(diff)) %>% 
  ungroup() %>% 
  separate(dyad, into = c("bat1", "bat2"))

# how many introduced dyads had increased similarity? ----
t <- 
  dd %>% 
  mutate(increase= diff>0) %>% 
  mutate(group= case_when(
    introduced ~ "introduced",
    bat1 %in% lp_bats & bat2 %in% lp_bats ~ "Las Pavas",
    bat1 %in% tole_bats & bat2 %in% tole_bats ~ "Tole")) %>% 
  group_by(group) %>% 
  summarize(increase= sum(increase), n=n()) %>% 
  mutate(decrease = n-increase) %>% 
  dplyr::select(-n)
# 21 of 35 introduced pairs showed convergence (60%)
# 1 of 10 familiar Las Pavas pairs showed convergence (10%)
# 9 of 21 familiar Tole bats showed convergence (43%)


# multi-membership model ----
library(brms)
fit <- brm(diff ~ introduced + (1|mm(bat1,bat2)), 
           family= gaussian(), 
           data= dd, 
           seed = 123,
           iter = 6000,    
           warmup = 2000, 
           chains = 4)

summary(fit)
# Compared to familiar pairs, introduced bats showed greater increases in call similarity (coefficient = 0.06, Bayesian 95% credible interval = 0.02 to 0.11) 


# get misclassification rates for pre- and post- intro ----

# Mundry's program requires the following:
# one column for test factor (groups)
# one column for control factor (individuals)
# however many variable columns
# variables should be only numbers, group and individual can be named
# NO missing values
# will take smallest number of calls

#for nested design, data do not have to be balanced
#tests for difference between groups ('testfac')
#groups and subjects ('contrfac') have to numbered consecutively and with integers beginning with 1

# make into function
pDFA <- function(xdata, test_fac, contr_fac, variables, n.sel = 100, nperm = 1000) {
  if (is.factor(test_fac)==F) {test_fac=as.factor(test_fac)}
  model=paste("lda(test_fac~",variables,", prior=pr_prob, subset=sel.index, data=xdata)",sep="")
  f.table=as.data.frame(table(contr_fac))
  ncf.levels=nrow(f.table)
  ntf.levels=nrow(table(test_fac))
  pr_prob=rep(1/ntf.levels, ntf.levels)#define prior probabilities to be equal for either of two groups
  #get number of cases per subject (level of contrfac):
  #get number of calls to select per subject (subject)
  n.to.sel=min(f.table$Freq)
  #set number of random selections original classification rate should be based on:
  ur.c.val=0
  ur.sel=0
  number=(1:nrow(xdata))
  #get assignment of subjects to groups
  gr=c()
  subj.gr= table(contr_fac,test_fac)
  subject=rownames(subj.gr)
  for (i in 1:ncf.levels){
    for (k in 1:ntf.levels){
      if (subj.gr[i,k]>0) {gr=c(gr,colnames(subj.gr)[k])}
    }
  }
  
  for (k in 1:n.sel){
    #make random selection of same number of cases per subject
    sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
    for (i in 1:ncf.levels){
      sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
    }
    sel.index= number[sel==1]
    #do a DFA and store results in 'res':
    res=eval(parse(text=model))
    #get predictions and store them in 'pred':
    pred=predict(res,xdata,prior=pr_prob)$class
    ur.sel=ur.sel+sum((test_fac==pred)[sel==1])
    ur.c.val= ur.c.val+ sum((test_fac==pred)[sel==0])
  }
  #save number of correctly classified calls in variable 'orig.res':
  ur.sel= ur.sel/ n.sel
  ur.c.val= ur.c.val/ n.sel
  
  #set P-value to 1 (since original data should be treated as 1 permutation):
  p.sel=1
  p.c.val=1
  all.corr=matrix(NA, nrow=nperm, ncol=2)
  all.corr[1,1]=ur.sel
  all.corr[1,2]=ur.c.val
  
  if (length(gr)==ncf.levels){
    for (k in 1:(nperm-1)){
      #randomize subjects' assignments to groups:
      r.gr=sample(gr,length(gr), replace=F)
      for (i in 1:length(subject)){
        test_fac[contr_fac==subject[i]]=r.gr[i]
      }
      
      #make random selection or same number of cases per subject
      sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
      for (i in 1:ncf.levels){
        sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
      }
      sel.index= number[sel==1]
      #do a DFA and store results in 'res':
      res=eval(parse(text=model))
      #get predictions and store them in 'pred':
      pred=predict(res,xdata,prior=pr_prob)$class
      ran.sel= sum((test_fac==pred)[sel==1])
      ran.c.val= sum((test_fac==pred)[sel==0])
      if (ran.sel>=ur.sel){p.sel = p.sel + 1}
      if (ran.c.val>= ur.c.val){p.c.val= p.c.val + 1}
      all.corr[k+1,1]=ran.sel
      all.corr[k+1,2]=ran.c.val
    }
    what=c("N correctly assigned, original, selected", "P for selected", "N correctly assigned, original, cross-validated", "P for cross-validated", "N groups (levels of test factor)", "N subjects (levels of control factor)", "N cases total", "N cases selected per subject","N selected total", "N permutations","N random selections")
    value=c(ur.sel,p.sel/nperm,ur.c.val,p.c.val/nperm,ntf.levels,ncf.levels,nrow(xdata),n.to.sel,n.to.sel*ncf.levels,nperm,n.sel)
    result=data.frame(what,value)
  }else{
    result="at least one subject is member of two groups; no test done"
  }
  result
  write.table(result,file=paste("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence-2025/06_test_convergence_over_time/",levels(test_fac)[1],"_",levels(test_fac)[2],sep=""),row.names=F,col.names=T)
  all.corr[,1]#comprises the number correctly classified selected objects for all
  #permutations (with the first value being that for the original data)
  all.corr[,2]#same for the cross-validated objects
  hist(all.corr[,1])#shows the frequency distribution
  print(result)
}

#### pre-introduction Tole bats vs post-intro Las Pavas bats

xdata <- pre
test_fac <- as.factor(as.character(xdata$group))
contr_fac <- as.factor(as.character(xdata$individual))
variables <- paste(colnames(pre)[3:ncol(pre)], collapse = "+")

pDFA(xdata, test_fac, contr_fac, variables)

#### post-introduction Tole bats vs post-intro Las Pavas bats

xdata <- post
test_fac <- as.factor(as.character(xdata$group))
contr_fac <- as.factor(as.character(xdata$individual))
variables <- paste(colnames(post)[3:ncol(post)], collapse = "+")

pDFA(xdata, test_fac, contr_fac, variables)

