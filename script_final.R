#last modified 30/11/2021 Tuesday, describeposteriors, centrality=mean
#change figures to the same scale
#refer to https://mc-stan.org/cmdstanr/ for installing cmdstanr
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#BiocManager::install("cmdstanr")
library(cmdstanr)
set_cmdstan_path("C:/Users/lulul/Documents/R/win-library/4.1/.cmdstanr/cmdstan-2.27.0/")
check_cmdstan_toolchain(fix=FALSE)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(16))
library(readxl)
library(brms)
library(loo)

###########################

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
#install.packages("rstan")
#setwd("C:/R_data/BayesianWorkshop") do not chnag eworking direcotry
wd="C:/Users/lulul/Documents/phd Thesis/LF/"
setwd(wd)
d <- read_xlsx(paste0(wd,"Fitness_final_analysis.xlsx"))
#View(d)
#include colony/individual level random effect, correlations between fitness, add source colony
#as genetic effect.
## Define factors
names(d)


d <- d %>% mutate(colony = factor(colony),
                  source_colony = factor(source_colony),
                  treatment = as.factor(ifelse(treatment==0,"control","restricted")),
                  osz = (Initial_colony_size - mean(Initial_colony_size))/sd(Initial_colony_size),
                  queen_fecundity = brood + w_young)

###### Model with additive overdispersion ######################################################

bf_survQ <- bf(queen_survival ~ treatment*osz +
                 (1|p|colony)) + bernoulli()
bf_survC <- bf(colony_survival ~ treatment*osz +
                 (1|p|colony)) + bernoulli()
bf_survW <- bf(w_old | trials(Initial_colony_size) ~  treatment*osz +
                 (1|p|colony)) + binomial()
bf_al <- bf(alates ~ treatment*osz +
              (1|p|colony)) + poisson()

bf_fe <- bf(queen_fecundity ~ osz*treatment +
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


model <- brm(bf_tot,
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
                file = paste0(wd,"model"),
                save_pars = save_pars(all = TRUE),
                data=d)
#takes about 10 minutes
sink("output_model_final.txt")
summary(model) ## Some significant positive correlations (due to shared fates if colony as a whole dies)
describe_posterior(model,ci=0.95,ci_method = "ETI",centrality = 'mean')
sink()


conditional_effects(model)

pp_check(model, resp = "queensurvival") # looks good, the black line in the middle is the real data
#others are predicted line
pp_check(model,resp = "colonysurvival")
pp_check(model, resp = "wold") # looks good
pp_check(model, resp = "queenfecundity") # looks good
pp_check(model, resp = "alates") # looks fine, regression predicted more 0 counts than actual

bayes_R2(model) ## pretty high 



#predict returns the fitted values before the inverse of the link function is applied (to return the data to the same
#scale as the response variable), and fitted shows it after it is applied.



#extract the fiitted regression line from  the model
#me=predict(m4.1.nos,re_formula = NA, scale = "response", robust=FALSE,#return mean estimate and sd
           #summary = T)%>% as_tibble() %>% bind_cols(d)
get_variables(model)#know what were estimaed and could be extracted draws


d.epred.c <- expand.grid(osz = seq(min(d$osz,na.rm=T),max(d$osz,na.rm = T),
                                   length.out = 100),
                         Initial_colony_size = d$Initial_colony_size,
                         treatment = levels(d$treatment)[1],
                         colony = levels(d$colony)[1]
)

d.epred.r <- expand.grid(osz = seq(min(d$osz,na.rm=T),max(d$osz,na.rm = T),
                                   length.out = 100),
                         Initial_colony_size = d$Initial_colony_size,
                         treatment = levels(d$treatment)[2],
                         colony = levels(d$colony)[1]
)




d_fit_qs_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                robust = FALSE,#return mean estimate and sd
               summary = T, resp = 'queensurvival')  %>% as_tibble() %>% bind_cols(d.epred.c) %>%
  mutate(variable=rep('queen survival',nrow(d.epred.c)))

d_fit_qs_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                   robust = FALSE,#return mean estimate and sd
                   summary = T, resp = 'queensurvival')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('queen survival',nrow(d.epred.r)))
  
d_fit_cs_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                   robust = FALSE,#return mean estimate and sd
                   summary = T, resp = 'colonysurvival')  %>% as_tibble()%>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('colony survival',nrow(d.epred.c)))

d_fit_cs_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'colonysurvival')  %>% as_tibble()%>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('colony survival',nrow(d.epred.r)))

d_fit_ws_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                   robust = FALSE,#return mean estimate and sd
                   summary = T, resp = 'wold')  %>% as_tibble() %>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('worker survival',nrow(d.epred.c)))


d_fit_ws_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'wold')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('worker survival',nrow(d.epred.r)))

d_fit_qf_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                   robust = FALSE,#return mean estimate and sd
                   summary = T, resp = 'queenfecundity')  %>% as_tibble() %>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('queen fecundity',nrow(d.epred.c)))

d_fit_qf_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'queenfecundity')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('queen fecundity',nrow(d.epred.r)))

d_fit_alates_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                   robust = FALSE,#return mean estimate and sd
                   summary = T, resp = 'alates')  %>% as_tibble()%>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('colony fecundity',nrow(d.epred.c)))

d_fit_alates_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                         robust = FALSE,#return mean estimate and sd
                         summary = T, resp = 'alates')  %>% as_tibble()%>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('colony fecundity',nrow(d.epred.r)))

dt= bind_rows(d_fit_alates_c,d_fit_alates_r, d_fit_cs_c,d_fit_cs_r, d_fit_qf_c,d_fit_qf_r,
              d_fit_qs_c,d_fit_qs_r,d_fit_ws_c,d_fit_ws_r)
#Valid response variables are: 'queensurvival', 'colonysurvival', 'wold', 'queenfecundity', 'alates'


##mee=fitted(m4.1.nos,re_formula = NA, scale = "response", robust=FALSE,
           #summary = T)%>% bind_cols(d)
#summary=T returns predicted response values, If summary = FALSE the 
#output resembles those of posterior_pred.brmsfit

#me=fitted(m4.1.nos,re_formula = NA, scale = "response",
#summary = F)%>% as_tibble() #results are returned on the scale of the response variable
#when fitted(summary=false), the returned values are the same as posterior_epred
#10000 rows= 10000 iterations, 34 observations for each variable



#plot worker survival:treatment
p5=ggplot(data=filter(dt,variable=='worker survival'),aes(x=osz, y=Estimate/Initial_colony_size,color=treatment))+
  geom_jitter(data=d,aes(x = osz, y = w_old/Initial_colony_size, color=treatment),width = 0.1, height = 0.1
              ,alpha=0.4)+
  theme_bw()+xlab("incipient colony size (z-score)")+ylab("worker survival")+ggtitle("f")+
  geom_line(lwd=1.0) + 
  geom_ribbon(aes(ymin=Q2.5/Initial_colony_size,ymax=Q97.5/Initial_colony_size,
                  fill=treatment),alpha=0.08,
              linetype='blank')+scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue")) + 
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)


#queen survival:treatment:osz
p3=ggplot(data=filter(dt,variable=='queen survival'),aes(x=osz, y=Estimate,color=treatment))+
  geom_jitter(data=d,aes(x = osz, y = queen_survival, color=treatment),width = 0.1, height = 0.1
              ,alpha=0.4)+
  theme_bw()+xlab("incipient colony size (z-score)")+ylab("queen survival")+ggtitle("d")+
  geom_line(lwd=1.0) + 
  geom_ribbon(aes(ymin=Q2.5,ymax=Q97.5,
                  fill=treatment),alpha=0.08,
              linetype='blank')+scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)

p4=ggplot(data=filter(dt,variable=='colony survival'),aes(x=osz, y=Estimate,color=treatment))+
  geom_jitter(data=d,aes(x = osz, y = colony_survival, color=treatment),width = 0.1, height = 0.1
              ,alpha=0.4)+
  theme_bw()+xlab("incipient colony size (z-score)")+ylab("colony survival")+ggtitle("e")+
  geom_line(lwd=1.0) + 
  geom_ribbon(aes(ymin=Q2.5,ymax=Q97.5,
                  fill=treatment),alpha=0.08,
              linetype='blank')+scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)
  

p8=ggplot(data=filter(dt,variable=='queen fecundity'),aes(x=osz, y=Estimate,color=treatment))+
  geom_jitter(data=d,aes(x = osz, y = queen_fecundity, color=treatment),width = 0.1, height = 0.1
              ,alpha=0.4)+
  theme_bw()+xlab("incipient colony size (z-score)")+ylab("queen fecundity")+ggtitle("c")+
  geom_line(lwd=1.0) + 
  geom_ribbon(aes(ymin=Q2.5,ymax=Q97.5,
                  fill=treatment),alpha=0.08,
              linetype='blank')+scale_color_manual(values=c("red", "blue"))+
  scale_y_continuous(n.breaks = 5)

p9=ggplot(data=filter(dt,variable=='colony fecundity'),aes(x=osz, y=Estimate,color=treatment))+
  geom_jitter(data=d,aes(x = osz, y = alates, color=treatment),width = 0.1, height = 0.1
              ,alpha=0.4)+
  theme_bw()+xlab("incipient colony size (z-score)")+ylab("colony fecundity (alates)")+ggtitle("d")+
  geom_line(lwd=1.0) + 
  geom_ribbon(aes(ymin=Q2.5,ymax=Q97.5,
                  fill=treatment),alpha=0.08,
              linetype='blank')+scale_color_manual(values=c("red", "blue"))+
  scale_y_continuous(n.breaks = 5)


#here use the the value when osz=0, coony is the first colony rather than the mean of all colonies and 
#all colony size
#treatment effect is when osz=0
d.epred.c <- expand.grid(osz = 0,
                         Initial_colony_size = d$Initial_colony_size,
                         treatment = levels(d$treatment)[1],
                         colony = levels(d$colony)[1]
)

d.epred.r <- expand.grid(osz = 0,
                         Initial_colony_size = d$Initial_colony_size,
                         treatment = levels(d$treatment)[2],
                         colony = levels(d$colony)[1]
)




d_fit_qs_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'queensurvival')  %>% as_tibble() %>% bind_cols(d.epred.c) %>%
  mutate(variable=rep('queen survival',nrow(d.epred.c)))

d_fit_qs_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'queensurvival')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('queen survival',nrow(d.epred.r)))

d_fit_cs_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'colonysurvival')  %>% as_tibble()%>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('colony survival',nrow(d.epred.c)))

d_fit_cs_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'colonysurvival')  %>% as_tibble()%>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('colony survival',nrow(d.epred.r)))

d_fit_ws_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'wold')  %>% as_tibble() %>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('worker survival',nrow(d.epred.c)))


d_fit_ws_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'wold')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('worker survival',nrow(d.epred.r)))

d_fit_qf_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'queenfecundity')  %>% as_tibble() %>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('queen fecundity',nrow(d.epred.c)))

d_fit_qf_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                     robust = FALSE,#return mean estimate and sd
                     summary = T, resp = 'queenfecundity')  %>% as_tibble() %>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('queen fecundity',nrow(d.epred.r)))

d_fit_alates_c <- fitted(model, newdata = d.epred.c, re_formula = NA, scale = "response",
                         robust = FALSE,#return mean estimate and sd
                         summary = T, resp = 'alates')  %>% as_tibble()%>% bind_cols(d.epred.c)%>%
  mutate(variable=rep('colony fecundity',nrow(d.epred.c)))

d_fit_alates_r <- fitted(model, newdata = d.epred.r, re_formula = NA, scale = "response",
                         robust = FALSE,#return mean estimate and sd
                         summary = T, resp = 'alates')  %>% as_tibble()%>% bind_cols(d.epred.r)%>%
  mutate(variable=rep('colony fecundity',nrow(d.epred.r)))

dt_treat= bind_rows(d_fit_alates_c,d_fit_alates_r, d_fit_cs_c,d_fit_cs_r, d_fit_qf_c,d_fit_qf_r,
              d_fit_qs_c,d_fit_qs_r,d_fit_ws_c,d_fit_ws_r)



p0=ggplot()+geom_jitter(data=d,aes(x = treatment, y = queen_survival),width = 0.1, height = 0.1
                        ,alpha=0.4)+
  theme_bw()+xlab("treatment")+ylab("queen survival")+ggtitle("a")+
  geom_pointrange(data=filter(dt_treat,variable=='queen survival'),
                  mapping=aes(x=treatment,y=Estimate,ymin=Q2.5,ymax=Q97.5),color='red')+
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)


p1=ggplot()+geom_jitter(data=d,aes(x = treatment, y = colony_survival),width = 0.1, height = 0.1
                        ,alpha=0.4)+
  theme_bw()+xlab("treatment")+ylab("colony survival")+ggtitle("b")+
  geom_pointrange(data=filter(dt_treat,variable=='colony survival'),
                  mapping=aes(x=treatment,y=Estimate,ymin=Q2.5,ymax=Q97.5),color='red')+
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)

p2=ggplot()+geom_jitter(data=d,aes(x = treatment, y =w_old/Initial_colony_size),width = 0.1, height = 0.1
                        ,alpha=0.4)+
  theme_bw()+xlab("treatment")+ylab("worker survival")+ggtitle("c")+
  geom_pointrange(data=filter(dt_treat,variable=='worker survival'),
                  mapping=aes(x=treatment,y=Estimate/Initial_colony_size,ymin=Q2.5/Initial_colony_size,ymax=Q97.5/Initial_colony_size),color='red')+
  scale_y_continuous(limits=c(-0.1, 1.1),n.breaks = 10)

p6=ggplot()+geom_jitter(data=d,aes(x = treatment, y = queen_fecundity),width = 0.1, height = 0.1
                        ,alpha=0.4)+
  theme_bw()+xlab("treatment")+ylab("queen fecundity")+ggtitle("a")+
  geom_pointrange(data=filter(dt_treat,variable=='queen fecundity'),
                  mapping=aes(x=treatment,y=Estimate,ymin=Q2.5,ymax=Q97.5),color='red')+
  scale_y_continuous(n.breaks = 5)

p7=ggplot()+geom_jitter(data=d,aes(x = treatment, y = alates),width = 0.1, height = 0.1
                        ,alpha=0.4)+
  theme_bw()+xlab("treatment")+ylab("colony fecundity (alates)")+ggtitle("b")+
  geom_pointrange(data=filter(dt_treat,variable=='colony fecundity'),
                  mapping=aes(x=treatment,y=Estimate,ymin=Q2.5,ymax=Q97.5),color='red')+
  scale_y_continuous(n.breaks = 5)
#using the colony size in the end generate significant results
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p3)

plot1=grid.arrange(arrangeGrob(p0,p1,p2,p3 + theme(legend.position="none"),
                         p4 + theme(legend.position="none"),
                         p5 + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(5, 1))
mylegend<-g_legend(p8)
plot2=grid.arrange(arrangeGrob(p6,p7,
                         p8 + theme(legend.position="none"),
                         p9 + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(5, 1))
ggsave("Fig2.png", plot = plot1, device = "png", dpi = 300, width = 10, height = 8, units = "in")
ggsave("Fig3.png", plot = plot2, device = "png", dpi = 300, width = 6, height = 8, units = "in")

