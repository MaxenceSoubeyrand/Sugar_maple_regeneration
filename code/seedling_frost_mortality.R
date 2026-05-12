#Script that model the mortality and frost injuries. 

rm(list=ls())

library(AICcmodavg)
library(tidyverse)
theme_set(theme_bw())
library(readxl)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggpubr)
library(binom)
library(ggeffects)

site <- data.frame(site=c("MUS", "HED", "MON", "KEK", "KIP"), 
                   latitude=c(49.932, 49.243, 48.460, 48.191, 46.740),
                   longitude=c(-78.698, -78.311, -79.418, -79.112, -78.905))

#Growth data
height <- read_excel("data/height.xlsx")

#c as cold injuries, and d and m as mortality. 

height_mort_cold <- height %>% 
  pivot_longer(cols=c(`HT 2005`, `HT 2006`, `HT 2007`), 
               names_to="year", values_to = "height") %>% 
  mutate(year=gsub("HT ", "", year)) %>% 
  mutate(height2=case_when(height=="c"~"cold injury",
                          height=="d"~"dead",
                          height=="m"~"dead",
                          .default = height)) %>% 
  group_by(site, seedling_ID, light, block) %>% 
  mutate(height3=case_when(height2=="cold injury"~lag(height2),
                          .default = height2)) %>%
  mutate(height3 = case_when(
    height2 == "cold injury" & any(height2 != "cold injury") ~ lag(height2),
    height2 == "cold injury" & all(height2 == "cold injury") ~ NA_character_, #If 3 cold injuries, NA
    TRUE ~ height2
  )) %>%
  # propagate dead to the following years
  mutate(height3 = {
    h <- height3
    for(i in 2:length(h)) {
      if(!is.na(h[i-1]) && h[i-1] == "dead") h[i] <- "dead"
    }
    h
  }) %>%
  mutate(height3 = case_when(
    height3 == "cold injury" & any(height3 != "cold injury") ~ lag(height3),
    TRUE ~ height3
  )) %>%
  mutate(height4 = case_when(
    height3 == "dead" ~ NA_character_, 
    TRUE ~ height3)) %>% 
  dplyr::select(site, seedling_ID, light, block, year, mort_cold=height2, height=height3) %>% 
  mutate(mort_cold=case_when(
    mort_cold %in% c("dead", "cold injury") ~ mort_cold,
    TRUE ~ "alive"
  )) %>% 
  filter(year!="2005") #Remowe 2005 because no mortality at the transplantation

#Summary
summary_ci_dead <- height_mort_cold %>%
  group_by(site, year, mort_cold) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(mort_cold!="alive") %>% 
  mutate(mort_cold=case_when(mort_cold=="cold injury"~"dead",
                             .default = mort_cold))

summary_ci_dead$site <- factor(summary_ci_dead$site, 
                               levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Figure 5.a
status_plot <- ggplot(summary_ci_dead, aes(x = year, y = n/360, fill = mort_cold)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~site, ncol=5) +
  scale_fill_manual(values = c("dead" = "black", "cold injury" = "deepskyblue")) +
  labs(x = "Year", y = "Mortality", fill = "Status") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

df_glm <- height_mort_cold %>%
  mutate(dead_bin=ifelse(mort_cold=="dead", 1, 0),
         dead_bin=ifelse(mort_cold=="cold injury", 1, dead_bin)) %>%
  group_by(site, seedling_ID, light, block) %>%
  mutate(
    dead_cumsum = cumsum(dead_bin),   
    keep = dead_cumsum <= 1         
  ) %>%
  filter(keep) %>%          
  dplyr::select(-dead_cumsum, -keep) %>%
  ungroup()

#Mean and condident intervals
mort_summary <- df_glm %>%
  group_by(light) %>%
  summarise(
    n = n(),
    dead = sum(dead_bin == 1),
    prop_dead = mean(dead_bin == 1)
  ) %>%
  ungroup()%>%
  rowwise() %>%
  mutate(
    ci = list(binom.confint(dead, n, conf.level = 0.90, methods = "exact"))
  ) %>%
  unnest(ci, names_sep = "_")

mort_summary

#Adding climate variables
clim <- readRDS("data/clim_growth.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

df_glm <- df_glm %>% 
  left_join(clim)

df_glm$light  <- case_when(
  df_glm$light  == "SH-1" ~ "Open-canopy",
  df_glm$light  == "SH-2" ~ "Intermediate",
  df_glm$light == "SH-3" ~ "Dense-canopy")
df_glm$light  <- factor(df_glm$light, levels = c("Dense-canopy", "Intermediate", "Open-canopy"))

df_glm <- df_glm %>% left_join(site)

#Random effect selection
mod_re1 <- glmer(dead_bin ~ latitude + light + 
                scale(temperature) + 
                scale(precipitation) + 
                scale(nb_late_frost_event) + 
                (1 | site/block) + (1 | seedling_ID) + (1 | year),
              data = df_glm, family = binomial) 

mod_re2 <- glmer(dead_bin ~ latitude +light + 
                   scale(temperature) + 
                   scale(precipitation) + 
                   scale(nb_late_frost_event) + 
                   (1 | site/block) + (1 | seedling_ID),
                 data = df_glm, family = binomial) 


mod_re3 <- glmer(dead_bin ~ latitude + light +
                   scale(temperature) + 
                   scale(precipitation) + 
                   scale(nb_late_frost_event) + 
                   (1 | site/block),
                 data = df_glm, family = binomial) 


mod_re4 <- glm(dead_bin ~ latitude + light +
                   scale(temperature) + 
                   scale(precipitation) + 
                   scale(nb_late_frost_event),
                 data = df_glm, family = binomial) 


AIC_re <- data.frame(mod=1:4, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3), AIC(mod_re4)))
AIC_re

#Fixed effect selection
mod1 <- glmer(dead_bin ~ latitude +light +
                     scale(temperature) + 
                     scale(precipitation) + 
                     scale(nb_late_frost_event) + 
                     light:scale(temperature) +
                     light:scale(precipitation) +
                     light:scale(nb_late_frost_event) +
                     scale(temperature):scale(precipitation) + 
                     (1|site/block),
                   data = df_glm, family = binomial)  

mod2 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                scale(precipitation) + 
                scale(nb_late_frost_event) + 
                light:scale(temperature) +
                light:scale(precipitation) +
                scale(temperature):scale(precipitation) + 
                (1|site/block),
              data = df_glm, family = binomial) 

mod3 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                scale(precipitation) + 
                scale(nb_late_frost_event) + 
                light:scale(temperature) +
                scale(temperature):scale(precipitation) + 
                (1|site/block),
              data = df_glm, family = binomial) 

mod4 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                scale(precipitation) + 
                scale(nb_late_frost_event) +
                scale(temperature):scale(precipitation) + 
                (1|site/block),
              data = df_glm, family = binomial) 

mod5 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                scale(precipitation) + 
                scale(nb_late_frost_event) +
                (1|site/block),
              data = df_glm, family = binomial) 

mod6 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                scale(precipitation) + 
                (1|site/block),
              data = df_glm, family = binomial) 

mod7 <- glmer(dead_bin ~ latitude +light +
                scale(temperature) + 
                (1|site/block),
              data = df_glm, family = binomial) 

mod8 <- glmer(dead_bin ~ latitude +light +
                (1|site/block),
              data = df_glm, family = binomial) 

mod9 <- glmer(dead_bin ~ latitude +
                (1|site/block),
              data = df_glm, family = binomial) 

mod10 <- glmer(dead_bin ~ 1 +
                (1|site/block),
              data = df_glm, family = binomial) 



liste_mod <- list(
  mod1 = mod1,
  mod2 = mod2,
  mod3 = mod3,
  mod4 = mod4,
  mod5 = mod5,
  mod6 = mod6,
  mod7 = mod7,
  mod8 = mod8,
  mod9 = mod9,
  mod10 = mod10)
aictab(liste_mod)

sel_mod <- mod2
#Model 2 is not convergiong but very close to the threshold, so it is acceptable. 

source("diagonistic_plot_function.R")

diagnostic_plots(sel_mod, df_glm)
ggsave(plot=diagnostic_plots(sel_mod, df_glm), filename="figures/SI_diag_plot_mortality.png", 
       width=8, height=10)

summary(sel_mod)
car::Anova(sel_mod, type=2)
# Conditional and marginal R²
r2 <- r.squaredGLMM(sel_mod)
r2

emm_light <- emmeans(
  sel_mod,
  ~ light,
  type = "response"
)

#  light * temperature predictions (figure 5.b)
gge_light_temp <- ggpredict(sel_mod, terms = c("temperature [all]"))

temp_dead <- ggplot(gge_light_temp, aes(x = x, y = predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Temperature (in °C)",
    y = "Probability of mortality and frost injury") +
  theme_minimal()+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

# light * precipitations predictions (figure 5.c)
gge_light_precip <- ggpredict(sel_mod, terms = c("precipitation [all]", "light"))

prec_dead <- ggplot(gge_light_precip, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Precipitation (in mm)",
    y = NULL,
    color = "Light treatment",
    fill = "Light treatment"
  ) +
  theme_minimal()+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

# light * late frost predictions (figure 5.d)
gge_light_frost <- ggpredict(sel_mod, terms = c("nb_late_frost_event [all]"))

frost_dead <- ggplot(gge_light_frost, aes(x = x, y = predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "# late frost events",
    y = NULL) +
  theme_minimal()+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))


pred_mort <- ggarrange(temp_dead, prec_dead, frost_dead, nrow=1, ncol=3, common.legend = T,
          legend="bottom", labels=c("b", "c", "d"), label.x = 0.1,
          label.y = 1.05, widths=c(1.07,1,1))

frost_dead_plot <- ggarrange(ggarrange(status_plot, labels=c("a"), label.x = 0.04), pred_mort, 
          nrow=2, ncol=1,
          heights=c(0.8,1))+
  bgcolor("white") +
  border("white")

ggsave(plot=frost_dead_plot, 
       filename="figures/frost_dead_plot.png", 
       width=8, height=6.2)