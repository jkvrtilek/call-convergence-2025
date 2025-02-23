# Contact call convergence in vampire bats
# Is vocal distance predicted by familiarity (group membership), within-group contact, or food sharing?
# using brms
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
# update this link after data are public
d <- read_csv('vocal_social_data.csv')

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
  ggplot(aes(x=kinship, y=sim))+
  geom_point(size=1, alpha=0.5)+
  geom_smooth(method= "lm")+
  ylab("contact call similarity")+
  xlab("kinship")+
  theme_bw())

# plot categories as boxplot
t %>% 
  ggplot(aes(x=kinship2, y=sim))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  ylab("contact call similarity")+
  xlab("kinship")+
  theme_bw()

# plot means and bootstrapped 95% CIs
means <- 
  t %>% 
  boot_ci2(x = .$kinship2, y= .$sim)
points <- 
  t %>% 
  mutate(effect= kinship2)
(p2 <- 
  means %>% 
  ggplot(aes(x=effect, y=mean))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  ylab("contact call similarity")+
  xlab("kinship range")+
  theme_bw())

# ***save kinship plot----
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
fitk <-
  brm(sim ~ 
        scale(kinship) +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fitk, ndraws=100)
summary(fitk)

# ***get model results----
(coefs.k <- 
  summary(fitk)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "kinship") %>% 
  mutate(n.pairs = nrow(t)) %>%  
  mutate(sample= "pairs of adults with known kinship"))

# get variables
get_variables(fitk)

# get samples from posterior distribution
pk <- 
  fitk %>% 
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
  mutate(contact.rate = ifelse(is.na(contact.rate), 0, contact.rate)) %>% 
  mutate(contact.lograte = ifelse(is.na(contact.lograte), 0, contact.lograte)) %>% 
  # label pairs that are familiar
  mutate(familiar = treatment_rank>1) %>%    
  # remove pairs with unknown kinship
  # note: this excludes all pairs from same roost on different years
  filter(!is.na(kinship))

# ***fit model for kinship after conditioning on social experience------------ 
fit.k2 <-
  brm(sim ~ 
        scale(contact.lograte) +
        scale(kinship)+
        familiar +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fit.k2, ndraws=100)
summary(fit.k2)

# ***get model results-----
(coefs.k2 <- 
  summary(fit.k2)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= case_when(
    term == 'scalecontact.lograte' ~ "contact rate (allogrooming + food sharing), conditioning on kinship and being caged together",
    term == 'scalekinship' ~ "kinship, conditioning on contact rate and being caged together",
    term == 'familiarTRUE' ~ "caged together, conditioning on contact rate and kinship")) %>% 
  mutate(n.pairs = nrow(t)) %>%  
  mutate(sample= "pairs of adults with known kinship"))

# plot model results
coefs.k2 %>% 
    # wrap text
    mutate(label= case_when(
      predictor == "contact rate (allogrooming + food sharing), conditioning on kinship and being caged together"~ 
        "contact rate\nconditioning on kinship\nand being caged together\n",
      predictor == "kinship, conditioning on contact rate and being caged together" ~ 
        "kinship\nconditioning on contact rate\nand being caged together",
      predictor == "caged together, conditioning on contact rate and kinship" ~
        "caged together\nconditioning on contact rate\nand kinship",
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
  spread_draws(b_scalekinship, b_scalecontact.lograte, b_familiarTRUE) %>%
  mutate(model = "Model2") %>% 
  pivot_longer(b_scalekinship:b_familiarTRUE, names_to = 'term', values_to= 'coeff')

remove(t)

# MODEL 3: Does time together among non-kin predict contact call similarity?----------

#***plot familiarity effect-----

# get relevant data to plot
t <- 
  d %>% 
  filter(both_adult) %>% 
  filter(dyad.sex== "female") %>% 
  filter(kinship <0.05| is.na(kinship))

# plot treatment categories for non-kin
t %>% 
  ggplot(aes(x=treatment_label, y=sim, color= treatment_label))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  geom_jitter(alpha=0.1, height=0, width=0.1)+
  ylab("contact call similarity")+
  xlab("familiarity level")+
  scale_colour_brewer(palette= "Dark2")+
  theme_bw()+
  theme(legend.position= "none")

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
    ggplot(aes(x=effect, y=mean, color=effect))+
    geom_jitter(data= points, aes(y= sim), size=1, alpha=0.5, height=0, width=0.1)+
    geom_boxplot(data= points, aes(y= sim), width=0.1, fill=NA, color="black", outlier.shape=NA)+
    geom_point(position = position_nudge(x = 0.25), size=1)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1)+
    ylab("contact call similarity")+
    xlab("familiarity level")+
    scale_colour_brewer(palette= "Dark2")+
    theme_bw()+
    theme(legend.position = "none"))

# plot without raw data
(p4 <- 
    means %>% 
    ggplot(aes(x=effect, y=mean, color=effect))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
    ylab("contact call similarity")+
    xlab(" \nfamiliarity level")+
    scale_colour_brewer(palette= "Dark2")+
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
  width = 9,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


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
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fit.intro, ndraws=100)
summary(fit.intro)

# get model results
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

#***fit same-roost effect for bats that were caged together-----

# are co-housed bats more similar when captured from the same roost?
# get pairs of nonkin females from different roost (treatment rank 3) or same roost Treatment rank 4) 
# compare those two groups treatment rank 3 vs 4 
# get relevant data to analyze
t <- 
  d %>% 
  filter(both_adult) %>% 
  # get female dyads
  filter(dyad.sex== "female") %>% 
  # get nonkin dyads
  filter(kinship <0.05) %>% 
  # get bats that never met or were introduced
  filter(treatment_rank == 3 | treatment_rank == 4 ) %>% 
  # label familiar as TRUE/FALSE
  mutate(same.roost = treatment_rank == 4)

# fit model with nonkin females that were from same roost or not
fit.roost <-
  brm(sim ~ 
        same.roost +
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fit.roost, ndraws=100)
summary(fit.roost)

# get model results
coefs.roost <- 
  summary(fit.roost)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
  mutate(predictor= "captured together") %>% 
  mutate(n.pairs= nrow(t)) %>% 
  mutate(sample= "pairs of unrelated female adults caged together")

# this effect is hard to interpret because some of the introduced bats were recorded prior to introduction which can bias the estimate
# the same bias cannot explain difference in similarity between familiarity levels 1 and 3

remove(t)

# MODEL 4: Does within-group contact rate predict contact call similarity?------------

#***plot contact effect-------------
# get relevant data
# pairs of female bats in the same group
t <- 
  d %>% 
  # exclude pairs with bats that could never groom
  filter(!is.na(contact.lograte)) %>%
  filter(kinship <0.05) %>% 
  filter(both_adult) %>% 
  filter(dyad.sex== "female") %>% 
  mutate(colony= substring(study_site,1,4)) %>% 
  mutate(captured.together= case_when(
    study_site == "2014_US" ~ TRUE,
    study_site == "2017_different sites" ~ FALSE,
    study_site == "2019_different sites" ~ FALSE,
    study_site == "2017_TL" ~ TRUE,
    study_site == "2019_TL" ~ TRUE,
    study_site == "2019_LB" ~ TRUE,
    study_site == "2019_CH" ~ TRUE,
    study_site == "2017_LP" ~ TRUE)) %>% 
  mutate(contact= contact.rate> 0)
 
# plot
(p <- 
  t %>% 
  mutate(`captured together` = captured.together) %>% 
  ggplot(aes(x=contact.lograte, y=sim,
             color=`captured together`, shape=`captured together`, linetype= `captured together`))+
    facet_wrap(~colony, scales= "free_x")+
    geom_point(size=2)+
  geom_smooth(method="lm")+
  xlab("within-group contact log rate")+
  ylab("contact call similarity")+
  scale_color_brewer(palette= "Dark2")+
  theme_bw()+
  theme(legend.position = "top"))

# save plot
ggsave(
  "contact_plot.pdf",
  plot = p,
  scale = 1,
  width = 6,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

  
# ***fit social grooming effect model ------------ 
fit.c <-
  brm(sim ~ 
        scale(contact.lograte) +
        captured.together+
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fit.c, ndraws=100)
summary(fit.c)

# get model results
(coefs.c <- 
  summary(fit.c)$fixed %>% 
  rownames_to_column(var= "term") %>% 
  filter(term!= "Intercept") %>% 
   
  mutate(predictor= "within-colony contact") %>%
    mutate(n.pairs= nrow(t)) %>%  
    mutate(sample= "pairs of unrelated female adults caged together"))

# get variables
get_variables(fit.c)

# get samples from posterior distribution
pc <- 
  fit.c %>% 
  spread_draws(b_scalecontact.lograte, b_captured.togetherTRUE) %>%
  mutate(model = "Model4") %>% 
  pivot_longer(b_scalecontact.lograte:b_captured.togetherTRUE, names_to = 'term', values_to= 'coeff')

remove(t)

# MODEL 5: Does within-group non-kin food sharing predict contact call similarity?------------
# conditioning on kinship

# get relevant data
t <- 
  d %>% 
  filter(both_adult) %>% 
  filter(kinship <0.05) %>% 
  filter(dyad.sex== "female") %>% 
  filter(!is.na(donation.rate)) %>% 
  mutate(captured.together= study_site != "2017_different sites") %>% 
  mutate(share = ifelse(donation.rate>0, "food-sharing\npairs", "non-sharing\npairs"))

#*** plot food sharing effect--------

# plot with captured together or not
(p1 <- 
  t %>% 
  ggplot(aes(x=donation.lograte, y=sim,
             color=treatment_label))+
    facet_wrap(~treatment_label, scales= "free_x")+
    geom_point(size=2)+
  geom_smooth(method="lm")+
  xlab("food sharing log rate")+
  ylab("contact call similarity")+
  scale_color_brewer(palette= "Dark2")+
  theme_bw()+
  theme(legend.position = "none"))

# plot means and 95% CIs
means <- 
  t %>% 
  boot_ci2(x = .$share, y= .$sim) 

points <- 
  t %>% 
  mutate(effect= share) 
(p2 <- 
    means %>% 
    mutate(effect= fct_rev(effect)) %>% 
    ggplot(aes(x=effect, y=mean))+
    geom_jitter(data= points, aes(y= sim, color= treatment_label), size=2, alpha=1, height=0, width=0.1)+
    geom_point(position = position_nudge(x = 0.25), size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), position = position_nudge(x = 0.25), size=1)+
    ylab("contact call similarity")+
    xlab("")+
    scale_color_brewer(palette= "Dark2")+
    theme_bw()+
    theme(axis.text.x = element_text(size=11))+
    theme(legend.position = "none"))

(p3 <- p1/p2 +plot_layout(widths= c(3,1)))

# save plot
ggsave(
  "sharing_plot.pdf",
  plot = p3,
  scale = 1,
  width = 5,
  height = 8,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# ***fit sharing effect model------------ 
fit.s <-
  brm(sim ~ 
        scale(donation.lograte)+
        captured.together+
        (1|mm(bat1,bat2)),
      data = t, 
      family = "beta",
      cores = nchains,
      chains = nchains,
      iter = chain_length,
      warmup = warmup_length)

# check model
pp_check(fit.s, ndraws=100)
summary(fit.s)

# get model results
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
  spread_draws(b_scaledonation.lograte, b_captured.togetherTRUE) %>%
  mutate(model = "Model5") %>% 
  pivot_longer(b_scaledonation.lograte:b_captured.togetherTRUE, names_to = 'term', values_to= 'coeff')

remove(t)

# COMPILE all coefficients ------
(all_coefs <- rbind(coefs.k,coefs.k2, coefs.intro,coefs.c, coefs.s))

# save coefficients as csv
all_coefs %>% 
  cbind(tibble(Model= c(1,2,2,2,3,4,4,5,5))) %>% 
  mutate(Sample= paste(n.pairs, sample)) %>% 
  mutate(Term= case_when(
    term == "scalekinship" ~ "kinship (scaled)",
    term == "scalecontact.lograte" ~ "contact log rate (scaled)",
    term == "familiarTRUE" ~ "caged together (T/F)",
    term == "famTRUE" ~ "caged together (T/F)",
    term == "scaledonation.lograte" ~ "food sharing log rate (scaled)",
    term== "captured.togetherTRUE" ~ "captured together (T/F)",
    TRUE ~ "error")) %>% 
  dplyr::select(Model, Sample, Term, Estimate:Tail_ESS) %>% 
  write.csv(file="model_results.csv")


# plot all coefficients--------
(models.plot <- 
  all_coefs %>% 
  filter(term != "captured.togetherTRUE") %>%  
  # wrap text
  mutate(label= case_when(
    predictor == "caged together" ~ 
      "5. caged together (model 3)\n(nonkin females only)",
    predictor == "within-colony contact" ~ 
      "6. contact rate\namong nonkin females\nconditioning on social history\n(model 4)",
    predictor == "within-colony food sharing"~ 
      "7. food-sharing rate\namong nonkin females\nconditioning on social history\n(model 5)",
    predictor == "contact rate (allogrooming + food sharing), conditioning on kinship and being caged together"~ 
      "3. contact rate (model 2)\nconditioning on kinship\nand being caged together\n",
    predictor == "kinship, conditioning on contact rate and being caged together" ~ 
      "2. kinship (model 2)\nconditioning on contact rate\nand being caged together",
    predictor == "caged together, conditioning on contact rate and kinship" ~
      "4. caged together (model 2)\nconditioning on contact rate\nand kinship",
    predictor == "kinship" ~
      "1. kinship (model 1)",
    TRUE ~ "error")) %>% 
  mutate(label= fct_rev(label)) %>%  
  ggplot(aes(x=Estimate, y=label))+
  geom_vline(xintercept = 0, linetype= "solid", color= 'black')+
  geom_point(size=2)+
  geom_errorbarh(aes(xmin=`l-95% CI`, xmax=`u-95% CI`, height=0.1), size=1)+
  xlab("effect on contact call similarity")+
  ylab("")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0)))

# save plot
ggsave(
  "summary_plot.pdf",
  plot = models.plot,
  scale = 1,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# COMPILE all posterior distributions ------
(all_post <- rbind(pk,pk2, pintro,pc, ps))

# plot all posterior distributions--------
 (models.plot2 <- 
   all_post %>% 
   filter(term != "b_captured.togetherTRUE") %>%  
   mutate(label= paste(model, term)) %>% 
   mutate(label= case_when(
     label == "Model1 b_scalekinship" ~
       "1. kinship (model 1)",
     label == "Model2 b_scalekinship" ~ 
       "2. kinship (model 2)\nconditioning on contact rate\nand co-housing",
     label == "Model2 b_scalecontact.lograte"~ 
       "3. contact rate (model 2)\nconditioning on kinship\nand co-housing\n",
     label == "Model2 b_familiarTRUE" ~
       "4. co-housing (model 2)\nconditioning on contact rate\nand kinship",
     label== "Model3 b_famTRUE" ~ 
       "5. co-housing (model 3)\n(nonkin females only)",
     label == "Model4 b_scalecontact.lograte" ~ 
       "6. contact rate (model 4) \nconditioning on capture site\n(familiar nonkin females only)",
     label == "Model5 b_scaledonation.lograte" ~ 
       "7. food-sharing rate (model 5) \nconditioning on capture site\n(familiar nonkin females only)",
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
  "summary_plot2.pdf",
  plot = models.plot2,
  scale = 1,
  width = 6,
  height = 7,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)


# save output with timestamp-----
timestamp <- substr(gsub(x=gsub(":","",Sys.time()), 
                         pattern=" ", replace="_"), start=1, stop=15)
timestamp
save.image(file= paste("model_fit_workspace_", timestamp, ".Rdata", sep=""))
