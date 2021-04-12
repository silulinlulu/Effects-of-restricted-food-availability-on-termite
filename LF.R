#============================================================================
#installation
#install.packages("brms")
#install.packages("cowplot")
#install.packages("bayestestR")
library(bayestestR)
#install.packages("viridis")
#if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
#BiocManager::install("cmdstanr")
library(BiocManager)
#install.packages("cmdstan")
#install.packages("devtools")
d#evtools::install_github("stan-dev/cmdstanr")
library(cmdstanr)
#check_cmdstan_toolchain(fix = T)
#install_cmdstan(dir="C:/.cmdstanr")
#set_cmdstan_path()

#dotR <- file.path("C:/R_installment/R-3.5.2/library", ".R")
#if (!file.exists(dotR))dir.create(dotR)
#M <- file.path(dotR, "Makevars.win")
#if (!file.exists(M)) file.create(M)
#cat("\nCXX14FLAGS=-O3 -march=native",
   # "CXX14 = C:/RBuildTools/3.5/mingw_64/bin/g++.exe -m$(WIN) -std=c++1y",
    #"CXX11FLAGS=-O3 -march=native",
    #file = M, sep = "\n", append = TRUE)
#install.packages(c("StanHeaders","rstan"),type="source")
#install.packages("car")
#install.packages("mcglm")
# we recommend running this is a fresh R session or restarting your current session
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#Sys.setenv(MAKEFLAGS = paste0("-j",parallel::detectCores()))
#example(stan_model, package = "rstan", run.dontrun = TRUE)

#============================================================================

setwd("C:/R_data/BayesianWorkshop")
library(tidyverse)
library(cowplot)
library(brms)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(readxl)
library(viridis)
library(cmdstanr)
library(dplyr)
library(reshape2)
d <- read_excel("fitness_LF vsHF_sent _Lulu.xls")
head(as.data.frame(d))
#Col_ID: colony id, source of genetic effects
#treatment: 0 is control and 2 is low food availability
#fitness_total_ok: all workers+alates
#survival: whether a colony is alive a the end of teh experiment
#alates_total: self-explainatory
#worker_nymphl: young workers+old workers
#worker_only: young workers+old workers-nymphs
#nymphs
#soldiers
#brood: eggs+larvaes
#ker_young: workers born during experiment=queen fecundity
#worker_old: large workers+nymphs
#transcriptome: whether we have transcriptome data from this colony or not

d$treatment[d$treatment==0]<-'control'
d$treatment[d$treatment==2]<-'LF' #low food
names(d)

#ggplot to have an idea of teh correlations
ggplot(d,aes(worker_young,brood,color=treatment)) +
  theme_cowplot(16) +
  geom_point() +
  facet_grid(~Col_ID) +
  background_grid(major = "xy", minor = "xy") +
  labs(x="worker_ol/colony_size") +
  geom_smooth(method = "gam", se=F, formula = y~s(x,k=3))
str(d)

dat<-melt(d, id=c("colony","Col_ID","treatment","survival","transcriptome" ))
ggplot(dat, aes(x=variable, y=value, fill=treatment)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free")
ggplot(d, aes(x=worker_old, y=brood, color=treatment)) + 
  geom_smooth()
ggplot(d, aes(x=worker_old, y=worker_young, color=treatment)) + 
  geom_smooth()

names(d)
summary(d)
d$fecundity<-d$brood+d$worker_young
d$sexuals<-d$worker_old+d$alates_total #old workers+nymphs and alates are sexuals
#Questions1: how did the treatment affect different fitness measurements
d <- d %>% mutate(Col_ID = factor(Col_ID),
                  treatment=factor(treatment))




car::scatterplotMatrix(~ survival+fitness_total_ok+alates_total+ worker_nymphl+sexuals+worker_only+brood+worker_young+worker_old| treatment ,
                       #ellipse = T,
                       col = viridis(3),
                       upper.panel = NULL,
                       gap = 0,
                       data = d)

car::scatterplotMatrix(~ fitness_total_ok+alates_total+ worker_nymphl+sexuals+worker_only+fecundity+worker_old| treatment ,
                       #ellipse = T,
                       col = viridis(3),
                       upper.panel = NULL,
                       gap = 0,
                       data = d)

car::scatterplotMatrix(~ survival+fitness_total_ok+alates_total+ worker_nymphl+sexuals+worker_only+brood+worker_young+worker_old| Col_ID ,
                       #ellipse = T,
                       col = viridis(3),
                       upper.panel = NULL,
                       gap = 0,
                       data = d)

car::scatterplotMatrix(~ fitness_total_ok+alates_total+ worker_nymphl+sexuals+worker_only+fecundity+worker_old| Col_ID ,
                       #ellipse = T,
                       col = viridis(3),
                       upper.panel = NULL,
                       gap = 0,
                       data = d)

View(d)

#difference/crediable interval of the difference
#2. waic to compare the model with one without thet treatment curve, to see whether 

#=====================model 2 is colony maturing earlier?  sexuals~colony size+treatment
bf_1 <- bf(fecundity ~ s(worker_nymphl, by=treatment)+treatment + (1|p|Col_ID)) + poisson()
bf_2 <- bf(alates_total ~ s(worker_nymphl, by=treatment) +treatment+ (1|p|Col_ID)) + poisson()



#bf1 <- bf_yw + bf_brood + bf_vi + set_rescor(rescor = F)
bf <- bf_1 + bf_2

prior_1 <- prior(normal(0,100), class = Intercept, resp = "fecundity") +
  prior(normal(0,1), class = b, resp = "fecundity")

prior_2 <- prior(normal(0,100), class = Intercept, resp = "alatestotal") +
  prior(normal(0,1), class = b, resp = "alatestotal")

my_priors <- prior_1+prior_2

m2 <- brm(bf,
          prior = my_priors,
          warmup = 2000,
          iter = 8000,
          chains = 4,
          cores = 4,
          #backend = "cmdstanr",
          #threads = threading(2),
          control = list(adapt_delta = 0.9999,
                         max_treedepth=15),
          seed = 123,
          file = "Csec_LF_3",
          data=d)

summary(m2)
bayestestR::describe_posterior(m2)
hypothesis(m2,"fecundity_Intercept>fecundity_treatmentLF")# fecundity=brood_young workers produced signif
#significantly lower in low food condition
hypothesis(m2,"alatestotal_Intercept<alatestotal_treatmentLF")# low food colonies produced significantly more
#alates=sexuals
hypothesis(m2,"alatestotal_sworker_nymphl:treatmentcontrol_1 < alatestotal_sworker_nymphl:treatmentLF_1",class = "bs")
hypothesis(m2,"alatestotal_sworker_nymphl.treatmentcontrol_1>0")
hypothesis(m2,"alatestotal_sworker_nymphl:treatmentLF_1>0",class = "bs")
parnames(m2)

## correlations all over the place
bayes_R2(m2)

plot(m2)
conditional_effects(m2)
brms::conditional_smooths(m2)

#conditional_effects(m1, "worker_old", resp = "workeryoung")

pp_check(m2, resp = "alatestotal")
pp_check(m2, resp = "fecundity")

difference/crediable interval of the difference

2. waic to compare the model with one without thet treatment curve, to see whether 
#====================================================================
#model3, cseclf4 changed sexual measures
#====================================================================
bf_1 <- bf(brood ~ s(worker_young, by=treatment)+treatment + (treatment|p|Col_ID)) + poisson()
bf_2 <- bf(sexuals ~ s(worker_young, by=treatment) +treatment+ (treatment|p|Col_ID)) + poisson()
#bf_ <- bf(viablepeas | trials(totalpeas) ~ pot*water + (1|p|block)) + binomial()

#bf1 <- bf_yw + bf_brood + bf_vi + set_rescor(rescor = F)
bf <- bf_1 + bf_2

prior_1 <- prior(normal(0,100), class = Intercept, resp = "brood") +
  prior(normal(0,1), class = b, resp = "brood")

prior_2 <- prior(normal(0,100), class = Intercept, resp = "sexuals") +
  prior(normal(0,1), class = b, resp = "sexuals")

#prior_vi <- prior(normal(0,100), class = Intercept, resp = "viablepeas") +
#prior(normal(0,1), class = b, resp = "viablepeas")

my_priors <- prior_1+prior_2

m3 <- brm(bf,
          prior = my_priors,
          warmup = 2000,
          iter = 8000,
          chains = 4,
          cores = 4,
          #backend = "cmdstanr",
          #threads = threading(2),
          control = list(adapt_delta = 0.99999,
                         max_treedepth=18),
          seed = 123,
          file = "Csec_LF_4",
          data=d)
## about 50 sec

summary(m3)
bayestestR::describe_posterior(m3)
hypothesis(m3,"sexuals_Intercept>sexuals_treatmentLF") #number of sexuals significantly lower in Low food
hypothesis(m3,"brood_Intercept<brood_treatmentLF")# number of brood did not differ between treatment
## correlations all over the place
bayes_R2(m3)

plot(m3)
conditional_effects(m3)
brms::conditional_smooths(m3)

#conditional_effects(m1, "worker_old", resp = "workeryoung")

pp_check(m2, resp = "alatestotal")
pp_check(m2, resp = "fecundity")
