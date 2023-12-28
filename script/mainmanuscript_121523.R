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

#prep data----

#load data frame
st <- read_csv("~/Social_Traps_092223.csv",col_types = cols(Treatment = col_character()))

#add binary
st_updated <- st %>%
  mutate(Trap_24_binary = case_when(
    Trapped_within == "1" ~ 1,
    Trapped_within == "3" ~ 1,
    Trapped_within == "5" ~ 1,
    Trapped_within == "24" ~ 1,
    is.na(Trapped_within) ~ 0))

#make treatment linear
st_linearTreament <- st_updated %>%
  mutate(Treatment_linear = as.numeric(Treatment))

#first convert the data to the long format

#make new trapped within column with anyone over 24hours dying at hour 25
st_status_change <- st_linearTreament %>%
  mutate(Trapped_within_2 = case_when(
    Trapped_within == "1" ~ 1,
    Trapped_within == "3" ~ 3,
    Trapped_within == "5" ~ 5,
    Trapped_within == "24" ~ 24,
    is.na(Trapped_within) ~ 25))

#make this a new status
st_status <- st_status_change %>%
  mutate(Status = case_when(
    Trapped_within_2 == "1" ~ 1,
    Trapped_within_2 == "3" ~ 1,
    Trapped_within_2 == "5" ~ 1,
    Trapped_within_2 == "24" ~ 1,
    Trapped_within_2 == "25" ~ 1))

#use survSplit to make these into distinct time points
St_surv <- survSplit(Surv(Trapped_within_2, Status)~., data=st_status, cut = c(1, 3, 5, 24), start='tstart', end='tstop')

#filter out the times with point 25 since these were just dummy variables
St_surv_filtered <-St_surv %>% filter(tstop != "25")

#add tstop columns
St_hazard <-St_surv_filtered %>% 
  mutate(
    d1 = case_when(
      tstop == "1" ~ 1,
      tstop == "3" ~ 0,
      tstop == "5" ~ 0,
      tstop == "24" ~ 0),
    d3 = case_when(
      tstop == "1" ~ 0,
      tstop == "3" ~ 1,
      tstop == "5" ~ 0,
      tstop == "24" ~ 0),
    d5 = case_when(
      tstop == "1" ~ 0,
      tstop == "3" ~ 0,
      tstop == "5" ~ 1,
      tstop == "24" ~ 0),
    d24 = case_when(
      tstop == "1" ~ 0,
      tstop == "3" ~ 0,
      tstop == "5" ~ 0,
      tstop == "24" ~ 1))


# fitting models for comparison -----

#basic model structure
fit11.11 <-
  brm(data = St_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d1 + d3 + d5 + d24 + Treatment_linear + Focal,
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      seed = 11)

fit11.11 <- add_criterion(fit11.11, c("loo", "waic"))


#add in random effect of date
fit11.12 <-
  brm(data = St_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d1 + d3 + d5 + d24 + Treatment_linear + Focal + (1|Date),
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      seed = 11)

fit11.12 <- add_criterion(fit11.12, c("loo", "waic"))

#add in arena
fit11.13 <-
  brm(data = St_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d1 + d3 + d5 + d24 + Treatment_linear + Focal + (1|Date) + (1|Arena),
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      seed = 11)

fit11.13 <- add_criterion(fit11.13, c("loo", "waic"))

#add in trap side
fit11.14 <-
  brm(data = St_hazard,
      family = binomial,
      Status | trials(1) ~ 0 + d1 + d3 + d5 + d24 + Treatment + Focal + Trap_side + (1|Date) + (1|Arena),
      prior(normal(0, 4), class = b),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      seed = 11)

fit11.14 <- add_criterion(fit11.14, c("loo", "waic"))

loo_compare(fit11.11, fit11.14, fit11.12,fit11.13, criterion = "loo") %>% print(simplify = F)
loo_compare(fit11.11,fit11.14, fit11.12, fit11.13, criterion = "waic") %>% print(simplify = F)

##date helps fit a little, arena and trap side does not, although all are very similar

#exploring model fit11.12
conditional_effects(fit11.12)

#estimating the differences between focals
fit11.12_emm <- fit11.12 %>% 
  emmeans::emmeans(~ Focal) %>% 
  emmeans::contrast(method = "pairwise") %>% 
  gather_emmeans_draws() %>% 
  mutate(.value = exp(.value)) 
###values
fit11.12_emm  %>% 
  median_hdci()

#summary stats ------

#count trials 
St_hazard %>%
  filter(d1 == "1") %>%
  count(Focal, Treatment)

#count number that died by 24
St_hazard %>%
  filter(d1 == "1") %>%
  count(Trapped_within == "24") 

#hazard ratio of dying
fixef(fit11.12)%>%exp()

#for discussion - from raw data - mean trap differences between treatments
st %>% 
  filter(Trapped_within != "NA") %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(Trapped_within))


#figure 1 ----- 
#conditional effects of focal and treatment on hazard

#extract conditional effects
ce <- conditional_effects(fit11.12, effects = "Treatment_linear:Focal")

#check basic plot
p2 <- plot(ce, plot = F)
p2

#make cleaner for publication
p2$Treatment_linear$data %>%
  mutate(across(estimate__:upper__, inv_logit_scaled)) %>%
  ggplot(aes(x = Treatment_linear)) +
  geom_line(aes(y = estimate__, color = Focal)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Focal), alpha = 0.1)+
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(color = "grey92"),
        legend.position = c(.225, .825)) +
  xlab("Number of drowned conspecifics") +
  scale_color_discrete(labels=c("red: mated female", "green: virgin female", "blue:mated male"))+
  scale_fill_discrete(labels=c("red: mated female", "green: virgin female", "blue:mated male"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab("Hazard of becoming trapped") 
