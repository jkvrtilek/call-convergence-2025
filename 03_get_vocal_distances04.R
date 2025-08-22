# Contact call convergence in vampire bats
# do linear discriminant function analysis (DFA) to classify calls to bat
# then use DFA get vocal distances from DFA
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# load packages
library(MASS)
library(patchwork)
library(tidyverse)

# Data for this script is available on Figshare: Vrtilek, Julia K.; Smith-Vidaurre, Grace; Carter, Gerald (2025). Data for "Vocal convergence during formation of cooperative relationships in vampire bats". figshare. Dataset. https://doi.org/10.6084/m9.figshare.29191334.v1

# get data
raw <- read.csv("vampire_call_measures_transformed.csv")

# tidy data and add sample sizes
d <- 
  raw %>% 
  separate(sound.files, into=c('ds','date', 'bat', 'file', 'call'), sep="_", remove = FALSE) %>% 
  group_by(bat) %>% 
  mutate(sample.size= n()) %>%
  ungroup()

# choose minimum number of calls per bat-----

# I chose this sample size based on the comparing accuracy per bat with and without cross-validation
# bats with more than 80 calls have correct classification rates similar to their cross-validated classification rates
# see plot "accuracy-with-and-without-cross-validation-5.pdf"
min.sample.size <- 100

# look at bats being excluded
d %>% 
  group_by(bat) %>% 
  summarize(n=n()) %>% 
  arrange(n) %>% 
  filter(n<min.sample.size) %>% 
  print(n=100)

# filter calls by min sample size
d2 <- 
  d %>% 
  filter(sample.size >= min.sample.size)

# count bats and calls per bat
d2 %>% 
  group_by(bat) %>% 
  summarize(n=n()) %>% 
  pull(n) %>% 
  mean() 


# plot n of bats by min number of calls per bat
# this takes about 2-3 minutes
if(TRUE){
  # gets sample sizes
  ss <-  1:1000
  n.bats <- rep(NA, length(ss)) 
  for (i in 1:length(ss)) {
    n.bats[i] <- 
      d %>% 
      filter(sample.size >= ss[i]) %>% 
      pull(bat) %>% 
      n_distinct()
    print(paste(i,"of",length(ss)))
  }
  # save
  tibble(ss, n.bats) %>% 
    write.csv("n_calls_by_n_bats.csv", row.names = F)  
}
t <- read.csv("n_calls_by_n_bats.csv")

(p <- 
  t %>%   
  ggplot(aes(x=ss, y=n.bats))+
    geom_point(size=1)+
    geom_line()+
    xlab("minimum number of calls necessary to include bat")+
  ylab("number of bats included in analysis")+
  geom_vline(xintercept= min.sample.size)+
  theme_bw())

# save plot
ggsave(
  filename = "n_calls_by_n_bats.pdf",
  plot = p,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# determine min/max/mean number of recording sessions per bat-------
# get list of dates bats were recorded on
rs <- d2 %>%
  dplyr::select(date, bat) %>% 
  distinct() %>% 
  arrange(bat)
rs$date <- as.Date(rs$date)

# add "session" column to account for recordings taking place over midnight
rs2 <- rs %>% 
  mutate(midnight = ifelse(bat != lag(bat), FALSE,
                           ifelse(lag(date) == date-1, TRUE,
                                  FALSE))) %>% 
  mutate(session = ifelse(midnight == TRUE, date - 1,
                          ifelse(midnight == FALSE, date,
                                 ifelse(is.na(midnight), date,
                                        NA)))) %>% 
  mutate(id = paste(bat, date, sep = "/"))
rs2$session <- as.Date(rs2$session, origin="1970-01-01")

# fix first session row
rs2$session[1] <- rs2$session[2]

# get number of sessions per bat
rs3 <- rs2 %>%
  # fix veronica 2013-02-15 (seems like no calls were recorded before midnight so above code fails)
  mutate(session = case_when(id == "veronica/2013-02-15" ~ as.Date("2013-02-14"),
                             TRUE ~ session)) %>% 
  dplyr::select(bat, session) %>% 
  group_by(bat) %>% 
  summarize(n = n()) %>% 
  arrange(n)
mean(rs3$n)

# classify calls to bat using dfa with cross-validation (leave one-out classification)--------
dfa <- lda(bat ~ 
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
           CV=T, 
           data=d2)

# get classification matrix
cm <- table(d2$bat, dfa$class)

# get overall correct classification rate (accuracy)
# this is the best accuracy estimate
correct.cases <- sum(diag(cm))
all.cases <- sum(cm)
accuracy <- correct.cases/all.cases
accuracy

# see median, mean, and range of correct assignments/all assignments to each bat
range(diag(cm)/colSums(cm), na.rm=T)
mean(diag(cm)/colSums(cm), na.rm=T)
median(diag(cm)/colSums(cm), na.rm=T)

# see median, mean and range of correct assignments/all calls from each bat
range(diag(cm)/rowSums(cm), na.rm=T)
mean(diag(cm)/rowSums(cm), na.rm=T)
median(diag(cm)/rowSums(cm), na.rm=T)

# save classification matrix
write.csv(cm, "cross-validated-dfa-classification.csv")

# plot accuracy as a function of sampling effort
tibble(bat= colnames(cm),
       accuracy= diag(cm)/rowSums(cm),
       n.cases= rowSums(cm)) %>%
  ggplot(aes(x=log(n.cases), y=accuracy))+
  geom_point(size=2)+
  geom_smooth(span=1)+
  xlab("log10 number of cases")+
  ylab('correct classification rate')+
  ggtitle('accuracy increase with number of calls per bat')+
  scale_colour_brewer(palette= "Dark2")+
  theme_bw()

# classify calls to bat using a single dfa without cross validation-------
# this method can give elevated accuracy (but doesn't matter because we are using this to get vocal distances)
dfa2 <- lda(bat ~ 
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
            CV=F, 
            data=d2)

# get classification rates
predictions <- predict(dfa2)
d2$prediction <- predictions$class
cm2 <- table(d2$bat, d2$prediction)

# get overall correct classification rate
correct.cases <- sum(diag(cm2))
all.cases <- sum(cm2)
accuracy2 <- correct.cases/all.cases
accuracy2

# compare accuracy of DFAs with and without CV------
# look at correlation between accuracy rates per bat (cross-validated or not)
n <- rowSums(cm)
c1 <- diag(cm)/n
c2 <- diag(cm2)/n

(p <- 
    tibble(c1,c2, n) %>% 
    mutate(sampling = ifelse(n>100, "over 100 calls", "less than 100 calls")) %>% 
    mutate(label= ifelse(sampling== "less than 100 calls", n, "")) %>% 
    ggplot(aes(x=c2, y=c1))+
    geom_smooth(method="lm")+
    geom_point(aes(color=sampling))+
    geom_text(aes(label=label), nudge_y=-0.02)+
    ylab("cross-validated accuracy assigning calls to bat")+
    xlab("accuracy assigning calls to bat\n(not cross-validated)")+
    scale_colour_brewer(palette= "Dark2")+
    theme_bw()+
    theme(legend.position =c(0.2,0.8), 
          legend.background= element_rect(linetype="solid", size= 0.1, colour ="black"))+
    ggtitle("correlation of accuracy by bat with and without cross validation", subtitle= paste("minimum calls per bat =", min.sample.size, "\ncorrelation =", round(cor(c1,c2),4))))
plotname <- paste0("accuracy-with-and-without-cross-validation-", min.sample.size, ".pdf")

# save plot
ggsave(
  filename = plotname,
  plot = p,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


# get vocal distances between bats---------

# get vocal distance as Mahalanobis distance between group centroids
meanGroup <- dfa2$means
distance <- as.matrix(dist(meanGroup %*% dfa2$scaling))

# get number of bats
n.bats <- n_distinct(d2$bat)
n.bats
# put number of bats in name of output
outname <- paste0("vocal-distance-",n.bats,"-bats.csv")

# save vocal distances
write.csv(distance, file= outname, row.names=F)  


# get DFA loadings---------------

# get proportion of variance explained by the discriminant functions
props.df <- dfa2$svd*100/sum(dfa2$svd) 
tibble(prop= props.df, x= 1:length(props.df)) %>% 
  ggplot(aes(x=x, y=prop))+
  geom_col()+
  xlab("discriminant function")+
  ylab("proportion of variance explained")

# get the coefficients of the linear discriminant functions (these tell you which variables were most important for identifying cases)
dfa2$scaling

# save all DFA loadings sorted by absolute value of DF1
loadings <- 
  dfa2$scaling %>% 
  as.data.frame() %>% 
  arrange(desc(abs(LD1)))

# put number of bats in name of output
outname2 <- paste0("dfa-loadings-",n.bats,"-bats.csv")

# save dfa loadings
write.csv(loadings, file= outname2, row.names=F)  


# OPTIONAL: compare observed accuracy to random classification rates------

# get expected
# re-run permutation test 
# this takes two hours
RUN_PERMUTATIONS <- FALSE
if(RUN_PERMUTATIONS){
  # get start time
  start <- Sys.time()
  
  # get observed accuracy
  obs <- accuracy2  
  
  # choose number of random datasets to test
  perms <- 1000
  
  # store results in vector 'exp' for expected rate (by chance)
  exp <-rep(NA, perms) 
  
  # repeat dfa with random data
  for (i in 1:perms) {
    
    # shuffle bat ids
    random.data <- 
      d2 %>% 
      mutate(bat= sample(bat))
    
    # classify calls to bat using dfa without cross validation
    dfa.rand <- lda(bat ~ 
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
                    data=random.data)
    
    # get classification rates
    predictions.rand <- predict(dfa.rand)
    d2$prediction.rand <- predictions.rand$class
    cm3 <- table(d2$prediction.rand, d2$bat)
    
    # get overall correct classification rate
    correct.cases <- sum(diag(cm3))
    all.cases <- sum(cm3)
    exp[i] <- correct.cases/all.cases
    
    # print progress
    print(paste(i, "of", perms))
  }
  
  # save results
  tibble(accuracy= c(obs,exp),
         label= c("observed", rep("expected", length(exp)))) %>%  
    write.csv("observed_v_expected_accuracy.csv")
  
  # runtime
  end <- Sys.time()
  end- start
  # perms = 1000 takes 2 hours
  
}

# get results
t <- read.csv("observed_v_expected_accuracy.csv")
exp <- t %>% filter(label== 'expected') %>% pull(accuracy)
obs <- t %>% filter(label== 'observed') %>% pull(accuracy)
  
# get 95% quantiles of expected values
exp.range <- round(quantile(exp, probs= c(0.025, 0.975), na.rm=T),5)

# plot expected
p1 <- 
  ggplot()+
  geom_histogram(aes(x=exp), color="black",fill="blue", binwidth = 0.00001)+
  xlab("expected values from null model")+
  theme_bw()+
  ggtitle("expected and observed accuracy", 
          subtitle = paste('observed = ',round(obs,3), ', 95% CI of expected: ', 
                           exp.range[1], ' to ', exp.range[2], "\np", ifelse(mean(exp>=obs)==0,
                                                                             paste("<",1/length(exp)), 
                                                                             paste("=",signif(mean(exp>=obs),digits=2))),", permutations = ",length(exp), sep=""))

# plot obs vs expected       
p2 <- 
  ggplot()+
  geom_histogram(aes(x=exp), color="black",fill="blue", binwidth = 0.0001)+
  geom_vline(aes(xintercept=obs), color="red", size=1)+
  xlab("expected values from null model (blue) and observed (red)")+
  theme_bw()

# combine
(p3 <- p1/p2)
# save plot
ggsave(
  "observed_vs_expected_accuracy.pdf",
  plot = p3,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


