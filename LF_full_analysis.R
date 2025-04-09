
#last modified 28/02/2022 Monday, model comparision to write in supplementary materials
####### analysis with cleaned-up data file Fitness_final_analysis.xlsx######
packageVersion("loo")
library(cmdstanr)
set_cmdstan_path("C:/Users/lulu/R/R-Library/cmdstan-2.27.0/")
check_cmdstan_toolchain(fix=TRUE)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(16))
library(readxl)
library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayestestR)
library(tidybayes)
library(bayesplot)
library(logspline)
library(ggpubr)
library(BiocManager)
library(gridExtra)
library(dplyr)

setwd("C:/R_data/BayesianWorkshop")
d <- read_xlsx("Fitness_final_analysis.xlsx")

names(d)

d <- d %>% mutate(
                  colony = factor(colony),
                  source_colony = factor(source_colony),
                  treatment = as.factor(ifelse(treatment==0,"control","restricted")),
                  queen_fecundity = brood + w_young,
                  osz = (Initial_colony_size - mean(Initial_colony_size))/sd(Initial_colony_size),
                  )
#m4.1.nos is the final model to compare with.
# regular model, starting model

###### Model with additive overdispersion, bernoulli, binomial, possion ######################################################

#==============================================
#model 4, a smooth term used for queen fecundity, modified based on model3
#==============================================
bf_survQ <- bf(queen_survival ~ treatment*osz + 
                 (1|p|colony)) + bernoulli()
bf_survC <- bf(colony_survival ~ treatment*osz + 
                 (1|p|colony)) + bernoulli()
bf_survW <- bf(w_old | trials(Initial_colony_size) ~  treatment*osz + 
                 (1|p|colony)) + binomial()
bf_fe <- bf(queen_fecundity ~ s(osz,by=treatment) + treatment +
            (1|p|colony)) + poisson()
bf_al <- bf(alates ~ treatment*osz + 
              (1|p|colony)) + poisson()


bf_tot <- bf_survQ + bf_survC + bf_survW + bf_fe + bf_al + set_rescor(rescor = F)


## priors
prior_survQ <- prior(normal(0,10), class = Intercept, resp= "queensurvival") +
  prior(normal(0,1), class = b, resp= "queensurvival")
prior_survC <- prior(normal(0,10), class = Intercept, resp= "colonysurvival") +
  prior(normal(0,1), class = b, resp= "colonysurvival")
prior_survW <- prior(normal(0,10), class = Intercept, resp= "wold") +
  prior(normal(0,1), class = b, resp= "wold")
prior_fe <- prior(normal(0,10), class = Intercept, resp= "queenfecundity") +
  prior(normal(0,1), class = b, resp= "queenfecundity")
prior_al <- prior(normal(0,10), class = Intercept, resp= "alates") +
  prior(normal(0,1), class = b, resp= "alates")

my_priors <- prior_survQ + prior_survC + prior_survW +  prior_fe + prior_al


bf_tot <- bf_survQ + bf_survC + bf_survW + bf_fe + bf_al + set_rescor(rescor = F)
m4.smooth <- brm(bf_tot,
          prior = my_priors,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          backend = "cmdstanr",
          threads = threading(2),
          control = list(adapt_delta = 0.999,
                         max_treedepth=10),
          seed = 666,
          file = "m4.smooth",
          data=d)

waic(m4.1.nos, m4.smooth)# The model without smooth term is better


#=======================
#zero_inflated Poisson
#=======================

bf_fe <- bf(queen_fecundity ~ treatment*osz + 
               (1|p|colony)) + zero_inflated_poisson()

bf_al <- bf(alates ~ treatment*osz + 
               (1|p|colony)) + zero_inflated_poisson()

bf_tot <- bf_survQ + bf_survC + bf_survW + bf_fe + bf_al + set_rescor(rescor = F)

m4_0inflat <- brm(bf_tot,
          prior = my_priors,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          backend = "cmdstanr",
          threads = threading(2),
          control = list(adapt_delta = 0.999,
                         max_treedepth=10),
          seed = 666,
          file = "m4_0inflat",
          data=d)

waic(m4.1.nos, m4_0inflat, m4.smooth)


#===========================================
### Additional random effect source colony
#===========================================
names(d)
bf_survQ <- bf(queen_survival ~ treatment*osz + 
                 (1|p|colony)+(1|q|source_colony)) + bernoulli()
bf_survC <- bf(colony_survival ~ treatment*osz + 
                 (1|p|colony)+(1|q|source_colony)) + bernoulli() 
bf_survW <- bf(w_old | trials(Initial_colony_size) ~  treatment*osz + 
                 (1|p|colony)+(1|q|source_colony)) + binomial()
bf_fe <- bf(queen_fecundity ~ osz*treatment  +
              (1|p|colony)+(1|q|source_colony)) + poisson()
bf_al <- bf(alates ~ treatment*osz + 
              (1|p|colony)+(1|q|source_colony)) + poisson()
bf_tot <- bf_survQ + bf_survC + bf_survW + bf_fe + bf_al + set_rescor(rescor = F)

m4_add_colonyEf <- brm(bf_tot,
          prior = my_priors,
          warmup = 1000,
          iter = 3500,
          chains = 4,
          cores = 4,
          backend = "cmdstanr",
          threads = threading(2),
          control = list(adapt_delta = 0.999,
                         max_treedepth=10),
          seed = 666,
          file = "m4_add_colonyEf",
          data=d)
describe_posterior(m4_add_colonyEf)
# one<-add_criterion(m4.1.nos, "loo",moment_match = TRUE)
# two<-add_criterion(m4_0inflat, "loo",moment_match = TRUE)
# three<-add_criterion(m4.smooth, "loo",moment_match = TRUE)
# four<-add_criterion(m4_add_colonyEf, "loo",moment_match = TRUE)

one<-add_criterion(model, "waic")
two<-add_criterion(m4_0inflat, "waic")
three<-add_criterion(m4.smooth, "waic")
four<-add_criterion(m4_add_colonyEf, "waic")
waic(one)#waic=353.8
waic(two)#waic=366.9
waic(three)#waic=354.3
waic(four)#waic=352.6
loo_compare(one, two, three, four, criterion="waic")

#the model with Poisson, without smooth and additional source colony effect is the best

