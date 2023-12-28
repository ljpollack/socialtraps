library(bayesplot)
library(brms)
library(broom)
library(flextable)
library(GGally)
library(ggmcmc)
library(ggrepel)
library(gtools)
library(loo)
library(patchwork)
library(psych)
library(survival)
library(cmdstanr)
library(tidyverse)
library(tidybayes)

#These models were prepared closely following a discrete time survival analysis tutorial for bayesian inference (Applied Longitudinal Data Analysis 2023, A.S. Kurz)



#clean and reorganize data for analysis ----

#load data frame
appendix <- read_csv("~/Appendix_Social_Traps_102623.csv")

appendix_updated <- appendix %>%
  mutate(Trap_120_binary = case_when(
    Trapped_within == "10" ~ 1,
    Trapped_within == "20" ~ 1,
    Trapped_within == "30" ~ 1,
    Trapped_within == "40" ~ 1,
    Trapped_within == "50" ~ 1,
    Trapped_within == "60" ~ 1,
    Trapped_within == "70" ~ 1,
    Trapped_within == "80" ~ 1,
    Trapped_within == "90" ~ 1,
    Trapped_within == "100" ~ 1,
    Trapped_within == "110" ~ 1,
    Trapped_within == "120" ~ 1,
    Trapped_within == "130" ~ 1,
    is.na(Trapped_within) ~ 0))

#make new trapped within column with anyone over 130 min dying at min 140
appendix_status_change <- appendix_updated %>%
  mutate(Trapped_within_2 = case_when(
    Trapped_within == "10" ~ 10,
    Trapped_within == "20" ~ 20,
    Trapped_within == "30" ~ 30,
    Trapped_within == "40" ~ 40,
    Trapped_within == "50" ~ 50,
    Trapped_within == "60" ~ 60,
    Trapped_within == "70" ~ 70,
    Trapped_within == "80" ~ 80,
    Trapped_within == "90" ~ 90,
    Trapped_within == "100" ~ 100,
    Trapped_within == "110" ~ 110,
    Trapped_within == "120" ~ 120,
    Trapped_within == "130" ~ 130,
    is.na(Trapped_within) ~ 140))

#make this a new status
appendix_status <- appendix_status_change %>%
  mutate(Status = case_when(
    Trapped_within_2 == "10" ~ 1,
    Trapped_within_2 == "20" ~ 1,
    Trapped_within_2 == "30" ~ 1,
    Trapped_within_2 == "40" ~ 1,
    Trapped_within_2 == "50" ~ 1,
    Trapped_within_2 == "60" ~ 1,
    Trapped_within_2 == "70" ~ 1,
    Trapped_within_2 == "80" ~ 1,
    Trapped_within_2 == "90" ~ 1,
    Trapped_within_2 == "100" ~ 1,
    Trapped_within_2 == "110" ~ 1,
    Trapped_within_2 == "120" ~ 1,
    Trapped_within_2 == "130" ~ 1,
    Trapped_within_2 == "140" ~ 1))

#use survSplit to make these into distinct time points
appendix_surv <- survSplit(Surv(Trapped_within_2, Status)~., data=appendix_status, cut = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), start='tstart', end='tstop')

#filter out the times with point 140 since these were just dummy
appendix_surv_filtered <-appendix_surv %>% filter(tstop != "140")

#add tstop columns
appendix_hazard <-appendix_surv_filtered %>% 
  mutate(
    d10 = case_when(
      tstop == "10" ~ 1,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d20 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 1,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d30 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 1,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d40 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 1,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d50 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 1,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d60 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 1,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d70 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 1,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d80 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 1,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d90 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 1,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d100 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 1,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d110 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 1,
      tstop == "120" ~ 0,
      tstop == "130" ~ 0),
    d120 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 1,
      tstop == "130" ~ 0),
    d130 = case_when(
      tstop == "10" ~ 0,
      tstop == "20" ~ 0,
      tstop == "30" ~ 0,
      tstop == "40" ~ 0,
      tstop == "50" ~ 0,
      tstop == "60" ~ 0,
      tstop == "70" ~ 0,
      tstop == "80" ~ 0,
      tstop == "90" ~ 0,
      tstop == "100" ~ 0,
      tstop == "110" ~ 0,
      tstop == "120" ~ 0,
      tstop == "130" ~ 1)
      )


#fitting models for comparison -----

#just time periods
fitk.1 <-
  brm(data = appendix_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d10 + d20 + d30 + d40 + d50 + d60 + d70 + d80 + d90 + d100 + d110 + d120 + d130,
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      seed = 11)

#include treatment
fitk.2 <-
  brm(data = appendix_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d10 + d20 + d30 + d40 + d50 + d60 + d70 + d80 + d90 + d100 + d110 + d120 + d130 + Treatment,
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 4000, warmup = 1000,
      seed = 11)

#include sex
fitk.3 <-
  brm(data = appendix_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d10 + d20 + d30 + d40 + d50 + d60 + d70 + d80 + d90 + d100 + d110 + d120 + d130 + Sex,
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 4000, warmup = 1000,
      seed = 11)

#include both treatment and sex
fitk.4 <-
  brm(data = appendix_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d10 + d20 + d30 + d40 + d50 + d60 + d70 + d80 + d90 + d100 + d110 + d120 + d130 + Treatment + Sex,
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 4000, warmup = 1000,
      seed = 11)

#compare models
model_weights(fitk.1, fitk.2, fitk.3, fitk.4, weights = "loo") %>% round(digits = 3) 

fitk.1  <- add_criterion(fitk.7, c("loo", "waic"))
fitk.2  <- add_criterion(fitk.8, c("loo", "waic"))
fitk.3  <- add_criterion(fitk.9, c("loo", "waic"))
fitk.4 <- add_criterion(fitk.10, c("loo", "waic"))

loo_compare(fitk.1,fitk.2, fitk.3, fitk.4, criterion = "loo") %>% print(simplify = F)
loo_compare(fitk.1,fitk.2, fitk.3, fitk.4, criterion = "waic") %>% print(simplify = F)
#treatment matters, sex seems not to have an impact on fit

#add in date to main model
fitk.11 <-
  brm(data = appendix_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d10 + d20 + d30 + d40 + d50 + d60 + d70 + d80 + d90 + d100 + d110 + d120 + d130 + Treatment + Sex + (1| Date),
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 4000, warmup = 1000,
      seed = 11)

model_weights(fitk.8, fitk.10, fitk.11, weights = "loo") %>% round(digits = 3)


#exploring main model outcomes ----

#hazard ratio of dying 
fixef(fitk.11)%>% exp() 

#in order to estimate median odds ratio, must filtering some variables
rg11 <- fitk.11 %>%
  emmeans::ref_grid(nuisance = c("d20", "d30", "d40", "d60", "d70", "d80", "d90", "d110", "d120"))

rg11@grid

#differences between treatments
fitk.11_emm <- fitk.11 %>% 
  emmeans::emmeans(~ Treatment,
                   nuisance = c("d20", "d30", "d40", "d60", "d70", "d80", "d90", "d110", "d120")) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
fitk.11_emm  %>% 
  median_hdci()

fitk.11 %>%
  emmeans::ref_grid(rg.limit = 50000)

#differences between sexes

fitk.11_emm <- fitk.11 %>% 
  emmeans::emmeans(~ Sex,
                   nuisance = c("d20", "d30", "d40", "d60", "d70", "d80", "d90", "d110", "d120")) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
fitk.11_emm  %>% 
  median_hdci()

##count trials 
appendix%>%
  count(Treatment, Sex)

##count number that died in all time periods
appendix %>%
  count(Trapped_within)


