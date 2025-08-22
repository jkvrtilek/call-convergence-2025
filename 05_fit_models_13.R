# Contact call convergence in vampire bats
# Is vocal distance predicted by familiarity (group membership), within-group affiliation, or food sharing?
# Julia Vrtilek and Gerry Carter

# clear workspace
rm(list=ls())

# set working directory
# setwd("/Users/gc1511/Dropbox (Personal)/Dropbox/_working/_ACTIVE/__students/Julia Vrtilek/2024/call convergence/2024-09") 

# load packages
library(tidyverse)
library(boot)
library(brms)
library(igraph)
library(bayesplot)
library(ggdist)
library(patchwork)
library(tidybayes)

# get data
d <- 
  read_csv('vocal_social_data.csv')

# look at social experience treatment
d %>% 
  group_by(treatment_rank) %>% 
  summarize(n.pairs= n())
# rank treatments 
# 1 = never met, 2540 pairs
# 2 = caught in same roost two years apart, 231 pairs
# 3 = introduced in captivity, 582 pairs
# 4 = caught together and caged together, 612 pairs
# 5 = long-term captive colony, 406 pairs


# make lists of bats that fit poorly into familiarity categories
keep.bats2to4 <- c("bat1","bat5","bat6","bat7")
rm.bats1to7 <- c("bat1","bat2","bat3","bat4","bat5","bat6","bat7")
# these 7 bats initially lived at the Chicago zoo and were moved to the University of Maryland
# bats 1,5,6,7 died in Maryland; bats 2-4 were then moved to OBC in Michigan

# because bats 1,5,6,7 never made it to the OBC, they only had the opportunity to share food with each other and bats 2-4
# but some of the bats at OBC were also born at the Chicago zoo, so bat 1,5,6,7 levels of familiarity with other zoo bats are uncertain
# we therefore exclude them from descriptive plots of kinship, familiarity, and food-sharing

# although bats 2-4 were added to the OBC colony, they were only there for a year
# for simplicity, we therefore exclude bats 2-4 in the plot of familiarity categories


# FUNCTIONS
# function to get 95% confidence interval of mean of vector x using classical bootstrapping
# argument 'bca = T' gives you bias-corrected and accelerated bootstrapping
boot_ci <- function(x, perms=5000, bca=F) {
  library(boot)
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

# get mean and 95% CIs via bootstrapping of values y within grouping variable x
# argument 'bca = T' gives you bias-corrected and accelerated bootstrapping
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000, bca=T){
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

# choose number of chains and chain length for Bayesian models in brms
# adjust as necessary
nchains <- 4
chain_length <- 6000
warmup_length <- 2000

# MODEL 1: Does kinship predict contact call similarity?--------

# get relevant data
# adult pairs with known kinship
t <- 
  d %>% 
  filter(both_adult) %>% 
  filter(!is.na(kinship)) 

# ***plot kinship effect----
# plot data
(p1 <- 
  t %>% 
   filter(!bat1 %in% keep.bats2to4) %>% 
   filter(!bat2 %in% keep.bats2to4) %>% 
  ggplot(aes(x=kinship, y=sim))+
  geom_point(size=1, alpha=0.5)+
  geom_smooth()+
  ylab("contact call similarity")+
  xlab("kinship")+
  theme_bw())

# plot means and bootstrapped 95% CIs
means <- 
  t %>% 
  filter(!bat1 %in% keep.bats2to4) %>% 
  filter(!bat2 %in% keep.bats2to4) %>% 
  boot_ci2(x = .$kinship2, y= .$sim)
(p2 <- 
  means %>% 
  ggplot(aes(x=effect, y=mean))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  ylab("contact call similarity")+
  xlab("kinship range")+
  theme_bw())

# save kinship plot
(kinship.plot <- p1+p2+plot_layout(widths= c(1,1)))
ggsave(
  "kinship_plot.pdf",
  plot = last_plot(),
  scale = 1,
  width = 8,
  height = 4,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# ***fit kinship model-----
fit.k <-
  brm(sim ~ 
        scale(kinship) +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# save posterior predictive check
(pp1 <- pp_check(fit.k, ndraws=100) + ggtitle("Model 1 posterior predictive check") +     
   theme(legend.position = "none", text = element_text(family = "Arial")))

# ***get model results----
(coefs.k <- 
  summary(fit.k)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "kinship") %>% 
  mutate(n.pairs = nrow(t)) %>%  
  mutate(sample= "pairs of adults with known kinship"))

# get variables
get_variables(fit.k)

# get samples from posterior distribution
pk <- 
  fit.k %>% 
  spread_draws(b_scalekinship) %>%
  mutate(model = "Model1") %>% 
  pivot_longer(b_scalekinship, names_to = 'term', values_to= 'coeff')

# remove temporary data
remove(t)

# MODEL 2: does kinship predict contact call similarity when conditioning on social experience?---------

# get relevant data
t <- 
  d %>% 
  filter(both_adult) %>% 
  # bats that could never have been seen sharing or grooming are labeled NA
  # for this analysis, convert NA for sharing and grooming to zeros
  mutate(affiliation.rate = ifelse(is.na(affiliation.rate), 0, affiliation.rate)) %>% 
  mutate(affiliation.lograte = ifelse(is.na(affiliation.lograte), 0, affiliation.lograte)) %>% 
  # label pairs that are familiar
  mutate(familiar = treatment_rank>1) %>%    
  # remove pairs with unknown kinship
  # note: this excludes all pairs from same roost on different years
  filter(!is.na(kinship))

# ***fit model for kinship after conditioning on social experience------------ 
fit.k2 <-
  brm(sim ~ 
        scale(affiliation.lograte) +
        scale(kinship)+
        familiar +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      seed =123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# save posterior predictive check
(pp2 <- pp_check(fit.k2, ndraws=100) + ggtitle("Model 2 posterior predictive check") +     
    theme(legend.position = "none", text = element_text(family = "Arial")))

# ***get model results-----
(coefs.k2 <- 
  summary(fit.k2)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= case_when(
    term == 'scaleaffiliation.lograte' ~ "affiliation rate (allogrooming + food sharing), conditioning on kinship and being caged together",
    term == 'scalekinship' ~ "kinship, conditioning on affiliation rate and being caged together",
    term == 'familiarTRUE' ~ "caged together, conditioning on affiliation rate and kinship")) %>% 
  mutate(n.pairs = nrow(t)) %>%  
  mutate(sample= "pairs of adults with known kinship"))

# plot model results
coefs.k2 %>% 
    # wrap text
    mutate(label= case_when(
      predictor == "affiliation rate (allogrooming + food sharing), conditioning on kinship and being caged together"~ 
        "affiliation rate\nconditioning on kinship\nand being caged together\n",
      predictor == "kinship, conditioning on affiliation rate and being caged together" ~ 
        "kinship\nconditioning on affiliation rate\nand being caged together",
      predictor == "caged together, conditioning on affiliation rate and kinship" ~
        "caged together\nconditioning on affiliation rate\nand kinship",
      TRUE ~ "error")) %>% 
    mutate(label= fct_rev(label)) %>%  
    ggplot(aes(x=Estimate, y=label))+
    geom_vline(xintercept = 0, linetype= "solid", color= 'black')+
    geom_point(size=2)+
    geom_errorbarh(aes(xmin=`l-95% CI`, xmax=`u-95% CI`, height=0.1), size=1)+
    xlab("effect on contact call similarity")+
    ylab("")+
    theme_bw()+
    theme(axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 0.5))

# get variables
get_variables(fit.k2)

# get samples from posterior distribution
pk2 <- 
  fit.k2 %>% 
  spread_draws(b_scalekinship, b_scaleaffiliation.lograte, b_familiarTRUE) %>%
  mutate(model = "Model2") %>% 
  pivot_longer(b_scalekinship:b_familiarTRUE, names_to = 'term', values_to= 'coeff')

remove(t)

# MODEL 3: Does time together among non-kin predict contact call similarity?----------

#***plot familiarity effect-----
unique(d$treatment_label)

# get relevant data to plot
t <- 
  d %>% 
  filter(!bat1 %in% rm.bats1to7) %>% 
  filter(!bat2 %in% rm.bats1to7) %>% 
  filter(both_adult) %>% 
  filter(dyad.sex== "female") %>% 
  filter(kinship <0.05| is.na(kinship)) %>% 
  # simplify treatment label and remove treatment 2, which is not really used in comparisons
  filter(treatment_rank != 2) %>% 
  mutate(treatment_label = substring(treatment_label, first= 4, last= 100)) %>% 
  # set order for plot
  mutate(treatment_label = factor(treatment_label, levels= c("never met", "different wild roost\nthen caged together", "same wild roost\nthen caged together","same long-term \ncaptive colony")))

# plot means and 95% CIs
means <- 
  t %>% 
  filter(!is.na(treatment_label)) %>% 
  boot_ci2(x = .$treatment_label, y= .$sim)

points <- 
  t %>% 
  filter(!is.na(treatment_label)) %>% 
  mutate(effect= treatment_label)

# plot with raw data
(p3 <- 
    means %>% 
    ggplot(aes(x=effect, y=mean))+
    geom_jitter(data= points, aes(y= sim), size=1, alpha=0.4, height=0, width=0.1, color= "darkgrey")+
    geom_boxplot(data= points, aes(y= sim), width=0.1, fill=NA, color="black", outlier.shape=NA)+
    geom_point(position = position_nudge(x = 0.25), size=1)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1)+
    ylab("contact call similarity")+
    xlab("familiarity level")+
    scale_colour_brewer(palette= "Dark2")+
    geom_hline(yintercept = means[1,"mean"], col = "black")+
    theme_bw()+
    theme(legend.position = "none"))

# plot without raw data
(p4 <- 
    means %>% 
    ggplot(aes(x=effect, y=mean))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
    ylab("contact call similarity")+
    xlab(" \nfamiliarity level")+
    scale_colour_brewer(palette= "Dark2")+
    geom_hline(yintercept = means[1,"mean"], col = "black")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 11))+
    theme(legend.position = "none"))

# combine plots
# remove x-axis from top plot
p3a <- 
  p3+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
(p5 <- p3a/p4)

# save plot
ggsave(
  "familiarity_plot.pdf",
  plot = p5,
  scale = 1,
  width = 7,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

remove(t)

#***fit model for introduction effect-------

# get relevant data to analyze
# get pairs of nonkin females that never met (treatment rank 1) and ones that were introduced (treatment rank 3)
# compare those two groups treatment rank 1 vs 3 
t <- 
  d %>% 
  filter(both_adult) %>% 
  # get female dyads
  filter(dyad.sex== "female") %>% 
  # get nonkin dyads
  filter(kinship == 0) %>% 
  # get bats that never met or were introduced
  filter(treatment_rank == 1 | treatment_rank == 3 ) %>% 
  # label familiar as TRUE/FALSE
  mutate(fam = treatment_rank > 1)

# fit model with nonkin females that were introduced together or not
fit.intro <-
  brm(sim ~ 
        fam +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# save posterior predictive check
(pp3 <- pp_check(fit.intro, ndraws=100) + ggtitle("Model 3 posterior predictive check") + 
    theme(legend.position = "none", text = element_text(family = "Arial")))

# ***get model results------
coefs.intro <- 
  summary(fit.intro)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "caged together") %>% 
  mutate(n.pairs= nrow(t)) %>% 
  mutate(sample= "pairs of unrelated female adults caught from different sites")

# get variables
get_variables(fit.intro)

# get samples from posterior distribution
pintro <- 
  fit.intro %>% 
  spread_draws(b_famTRUE) %>%
  mutate(model = "Model3") %>% 
  pivot_longer(b_famTRUE, names_to = 'term', values_to= 'coeff')

remove(t)

# MODEL 4: Does female within-group affiliation rate predict contact call similarity?------------

#***plot affiliation effect-------------
# get relevant data
# pairs of female bats in the same group
t <- 
  d %>% 
  # exclude pairs with bats that could never groom
  filter(!is.na(affiliation.lograte)) %>%
  filter(kinship <0.05) %>% 
  filter(both_adult) %>% 
  filter(dyad.sex== "female") %>% 
  mutate(colony= substring(study_site,1,4)) %>% 
  mutate(affiliation= affiliation.rate> 0) %>% 
  mutate(type= ifelse(colony== "2014", "food sharing in stable colony", "mostly allogrooming in merged colony"))

# plot
(p <- 
  t %>% 
    filter(!bat1 %in% keep.bats2to4) %>% 
    filter(!bat2 %in% keep.bats2to4) %>% 
  ggplot(aes(x=affiliation.lograte, y=sim, color= colony))+
    facet_wrap(~type, scales= "free_x")+
    geom_point(size=2)+
  geom_smooth(method="lm")+
  xlab("within-group affiliation log rate")+
  ylab("contact call similarity")+
  scale_color_brewer(palette= "Dark2")+
  theme_bw()+
  theme(legend.position = "top"))

# save plot
ggsave(
  "affiliation_plot.pdf",
  plot = p,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# ***fit affiliation model ------------ 
fit.aff <-
  brm(sim ~ 
        scale(affiliation.lograte) +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# save posterior predictive check
(pp4 <- pp_check(fit.aff, ndraws=100) + ggtitle("Model 4 posterior predictive check") + 
    theme(legend.position = "none", text = element_text(family = "Arial")))

# ***get model results----
(coefs.aff <- 
  summary(fit.aff)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
   
  mutate(predictor= "within-colony affiliation") %>%
    mutate(n.pairs= nrow(t)) %>%  
    mutate(sample= "pairs of unrelated female adults caged together"))

# get variables
get_variables(fit.aff)

# get samples from posterior distribution
pc <- 
  fit.aff %>% 
  spread_draws(b_scaleaffiliation.lograte) %>%
  mutate(model = "Model4") %>% 
  pivot_longer(b_scaleaffiliation.lograte, names_to = 'term', values_to= 'coeff')

remove(t)

# MODEL 5: Does female within-group non-kin food sharing predict contact call similarity?------------
# conditioning on kinship

# get relevant data
t <- 
  d %>% 
  filter(both_adult) %>% 
  filter(kinship <0.05) %>% 
  filter(dyad.sex== "female") %>% 
  mutate(colony= substring(study_site,1,4)) %>% 
  filter(!is.na(donation.rate)) %>% 
  mutate(share = ifelse(donation.rate>0, "food-sharing\nfemale pairs", "non-sharing\nfemale pairs")) %>% 
  mutate(share = factor(share)) %>% 
  mutate(share = fct_rev(share)) %>% 
  mutate(type= ifelse(colony== "2014", "long-term stable colony (2014)", "merged colony (2017)"))

# get data from all pairs (not just female pairs)
t2 <- 
  d %>% 
  filter(both_adult) %>% 
  filter(kinship <0.05) %>% 
  mutate(colony= substring(study_site,1,4)) %>% 
  filter(!is.na(donation.rate)) %>% 
  mutate(share = ifelse(donation.rate>0, "food-sharing\npairs", "non-sharing\npairs")) %>% 
  mutate(share = factor(share)) %>% 
  mutate(share = fct_rev(share)) %>% 
  mutate(type= ifelse(colony== "2014", "long-term stable colony (2014)", "merged colony (2017)"))

#***plot food sharing effect--------

# plot
(p1 <- 
  t %>% 
   filter(!bat1 %in% keep.bats2to4) %>% 
   filter(!bat2 %in% keep.bats2to4) %>% 
  ggplot(aes(x=donation.lograte, y=sim))+
    facet_wrap(~type, scales= "free_x")+
    geom_point(size=2, alpha=0.4)+
  geom_smooth(method="lm")+
  xlab("food sharing log rate")+
  ylab("contact call similarity")+
  scale_color_brewer(palette= "Dark2")+
  theme_bw()+
  theme(legend.position = "none"))

# plot means and 95% CIs
means <- 
  t %>% 
  filter(!bat1 %in% keep.bats2to4) %>% 
  filter(!bat2 %in% keep.bats2to4) %>% 
  boot_ci2(x = .$share, y= .$sim) 
points <- 
  t %>% 
  filter(!bat1 %in% keep.bats2to4) %>% 
  filter(!bat2 %in% keep.bats2to4) %>% 
  mutate(effect= share) 
(p2 <- 
    means %>% 
    mutate(effect= fct_rev(effect)) %>% 
    mutate(label= "female nonkin pairs") %>% 
    ggplot(aes(x=effect, y=mean))+
    facet_wrap(~ label)+
    geom_jitter(data= points, aes(y= sim), size=2, alpha=0.4, height=0, width=0.1)+
    geom_point(position = position_nudge(x = 0.25), size=3, color= "blue")+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1, color= "blue")+
    xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(size=11),
          axis.text.y= element_blank(),
          axis.title.y = element_blank()))

# same plot with both sexes
means2 <- 
  t2 %>% 
  filter(!bat1 %in% keep.bats2to4) %>% 
  filter(!bat2 %in% keep.bats2to4) %>% 
  boot_ci2(x = .$share, y= .$sim) 
points2 <- 
  t2 %>% 
  filter(!bat1 %in% keep.bats2to4) %>% 
  filter(!bat2 %in% keep.bats2to4) %>% 
  mutate(effect= share) 
(p3 <- 
    means2 %>% 
    mutate(label= "all nonkin pairs") %>% 
    ggplot(aes(x=effect, y=mean))+
    facet_wrap(~ label)+
    geom_jitter(data= points2, aes(y= sim), size=2, alpha=0.4, height=0, width=0.1)+
    geom_point(position = position_nudge(x = 0.25), size=3, color= "blue")+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1, color= "blue")+
    xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(size=11),
          axis.text.y= element_blank(),
          axis.title.y = element_blank()))


(p4 <- p1+p2+p3 +plot_layout(widths= c(2,1,1))+ plot_annotation(tag_levels = c("A", "B", "C")))

# save plot
ggsave(
  "sharing_plot.pdf",
  plot = p4,
  scale = 1,
  width = 10,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# ***fit sharing effect model------------ 
fit.s <-
  brm(sim ~ 
        scale(donation.lograte)+
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      seed = 123,
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# save posterior predictive check
(pp5 <- pp_check(fit.s, ndraws=100) + ggtitle("Model 5 posterior predictive check") + 
    theme(legend.position = "none", text = element_text(family = "Arial")))

# ***get model results-----
(coefs.s <- 
  summary(fit.s)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "within-colony food sharing") %>% 
  mutate(n.pairs= nrow(t)) %>% 
  mutate(sample= "pairs of unrelated female adults caged together"))

# get variables
get_variables(fit.s)

# get samples from posterior distribution
ps <- 
  fit.s %>% 
  spread_draws(b_scaledonation.lograte) %>%
  mutate(model = "Model5") %>% 
  pivot_longer(b_scaledonation.lograte, names_to = 'term', values_to= 'coeff')

remove(t)


# COMPILE all coefficients ------
(all_coefs <- rbind(coefs.k,coefs.k2, coefs.intro,coefs.aff, coefs.s))

# save coefficients as csv
all_coefs %>% 
  cbind(tibble(Model= c(1,2,2,2,3,4,5))) %>% 
  mutate(Sample= paste(n.pairs, sample)) %>% 
  mutate(Term= case_when(
    term == "scalekinship" ~ "kinship (scaled)",
    term == "scaleaffiliation.lograte" ~ "affiliation log rate (scaled)",
    term == "familiarTRUE" ~ "caged together (T/F)",
    term == "famTRUE" ~ "caged together (T/F)",
    term == "scaledonation.lograte" ~ "food sharing log rate (scaled)",
    TRUE ~ "error")) %>% 
  dplyr::select(Model, Sample, Term, Estimate:Tail_ESS) %>% 
  write.csv(file="model_results.csv")

# COMPILE all posterior distributions ------
(all_post <- rbind(pk,pk2, pintro,pc, ps))

# plot all posterior distributions--------
(models.plot <- 
   all_post %>% 
   filter(term != "b_captured.togetherTRUE") %>%  
   mutate(label= paste(model, term)) %>% 
   mutate(label= case_when(
     label == "Model1 b_scalekinship" ~
       "1. model 1: kinship",
     label == "Model2 b_scalekinship" ~ 
       "2. model 2: kinship\n(covariates: affiliation\nand co-housing)",
     label == "Model2 b_scaleaffiliation.lograte"~ 
       "3. model 2: affiliation\n(covariates: kinship\nand co-housing)",
     label == "Model2 b_familiarTRUE" ~
       "4. model 2: co-housing\n(covariates: affiliation\nand kinship)",
     label== "Model3 b_famTRUE" ~ 
       "5. model 3: co-housing\n(only nonkin females)",
     label == "Model4 b_scaleaffiliation.lograte" ~ 
       "6. model 4: affiliation \n(only co-housed\nnonkin females)",
     label == "Model5 b_scaledonation.lograte" ~ 
       "7. model 5: food sharing\n(only co-housed\nnonkin females)",
     TRUE ~ "error")) %>% 
   mutate(label= fct_rev(label)) %>%  
   ggplot(aes(y = label, x = coeff, fill=model)) +
   stat_halfeye(.width = c(0.95), linewidth= 5, size=5)+
   geom_vline(xintercept = 0)+
   ylab("")+
   xlab("estimated effect on vocal similarity (coefficient)")+
   theme_bw()+
   theme(legend.position= 'none', 
         axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5)))

# save plot
ggsave(
  "summary_plot.pdf",
  plot = models.plot,
  scale = 1,
  width = 5.5,
  height = 6.5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# make overall posterior predictive check figure-----
pp <- pp1 / pp2 / pp3 / pp4 / pp5
pp

# save output with timestamp-----
timestamp <- substr(gsub(x=gsub(":","",Sys.time()), 
                         pattern=" ", replace="_"), start=1, stop=15)
timestamp
save.image(file= paste("model_fit_workspace_", timestamp, ".Rdata", sep=""))
