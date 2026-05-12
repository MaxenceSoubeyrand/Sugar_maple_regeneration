#Script that model the height growth

rm(list=ls())

library(tidyverse)
theme_set(theme_bw())
library(readxl)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggpubr)
library(AICcmodavg)
library(ggeffects)

site <- data.frame(site=c("MUS", "HED", "MON", "KEK", "KIP"), 
                   latitude=c(49.932, 49.243, 48.460, 48.191, 46.740),
                   longitude=c(-78.698, -78.311, -79.418, -79.112, -78.905))

#Growth data
height <- read_excel("data/height.xlsx")

head(height)

#Manipulation to have 0 cm growth when cold injuries.
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
    height2 == "cold injury" & all(height2 == "cold injury") ~ NA_character_, 
    TRUE ~ height2
  )) %>%
  # propager "dead" aux années suivantes, en vérifiant NA
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
  dplyr::select(site, seedling_ID, light, block, year, mort_cold=height2, height=height4) %>%
  na.omit() 

#Growth calcul

growth <- height_mort_cold %>% 
  group_by(site, seedling_ID, block, light) %>% 
  mutate(RGR=log(as.numeric(height))-log(lag(as.numeric(height))),
         RGR=case_when(year=="2006"~RGR/2,
                          year=="2007"~RGR)) %>% 
  mutate(growth=as.numeric(height)-lag(as.numeric(height)),
         growth=case_when(year=="2006"~growth/2,
                       year=="2007"~growth)) %>% 
  na.omit()


growth$site <- factor(growth$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Climate data
clim <- readRDS("data/clim_growth.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

#Join climate and height data
growth <- growth %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(seedling_ID=paste(site, seedling_ID, light, block, sep="_")) %>% 
  mutate(site_block=paste0(site, "_", block)) %>% 
  mutate(height=as.numeric(height)) %>% 
  na.omit()

growth$light  <- case_when(
  growth$light  == "SH-1" ~ "Open-canopy",
  growth$light  == "SH-2" ~ "Intermediate",
  growth$light == "SH-3" ~ "Dense-canopy")
growth$light  <- factor(growth$light, levels = c("Dense-canopy", "Intermediate", "Open-canopy"))

growth <- growth %>% left_join(site)

#Random effect model selection
mod_re1 <- lmer(log(growth+1)~ latitude + scale(height) + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) +
                (1 | site_block) + (1 | seedling_ID) + (1 | year), data=growth)

mod_re2 <- lmer(log(growth+1)~ latitude + scale(height) + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) + 
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_re3 <- lmer(log(growth+1)~ latitude + scale(height) + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) + 
                  (1 | site_block) , data=growth)

mod_re4 <- lm(log(growth+1)~ latitude + scale(height) + light + 
                scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) , data=growth)

AIC_re <- data.frame(mod=1:4, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3), AIC(mod_re4)))
AIC_re

#Fixed effect model selection
mod_me1 <- lmer(log(growth+1)~ latitude + scale(height) + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) +
                  
                  scale(temperature):scale(precipitation) +
                  light:scale(height) + light:scale(temperature) + 
                  light:scale(precipitation) + light:scale(nb_late_frost_event)+
                  (1 | site_block) + (1 | seedling_ID) , data=growth)

mod_me2 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) +
                  
                  scale(temperature):scale(precipitation) +
                  light:scale(height) + light:scale(temperature)+
                  light:scale(precipitation) +
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me3 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) +
                  scale(temperature):scale(precipitation) +
                light:scale(height) + light:scale(temperature)+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me4 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(nb_late_frost_event) + scale(precipitation) +
                  scale(temperature):scale(precipitation) +
                  light:scale(height)+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me5 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(precipitation) + scale(nb_late_frost_event)+
                  scale(temperature):scale(precipitation) +
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me6 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(precipitation) + scale(nb_late_frost_event)+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me7 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature) + scale(precipitation)+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me8 <- lmer(log(growth+1)~ latitude + light + 
                  scale(temperature)+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me9 <- lmer(log(growth+1)~ latitude + light+
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me10 <- lmer(log(growth+1)~ latitude +
                  (1 | site_block) + (1 | seedling_ID), data=growth)

mod_me11 <- lmer(log(growth+1)~1+
                   (1 | site_block) + (1 | seedling_ID), data=growth)

liste_mod <- list(
  mod_me1 = mod_me1,
  mod_me2 = mod_me2,
  mod_me3 = mod_me3,
  mod_me4 = mod_me4,
  mod_me5 = mod_me5,
  mod_me6 = mod_me6,
  mod_me7 = mod_me7,
  mod_me8 = mod_me8,
  mod_me9 = mod_me9,
  mod_me10 = mod_me10,
  mod_me11 = mod_me11)

aictab(liste_mod)

sel_mod <- mod_me1

summary(sel_mod)


source("diagonistic_plot_function.R")

#plot diagnostics
ggsave(plot=diagnostic_plots(sel_mod, growth), filename="figures/SI_diag_plot_growth.png", 
       width=8, height=10)


#anova on the model with most effect
anova(sel_mod)

#Marginal and conditionnal R²
r2 <- r.squaredGLMM(sel_mod)
r2

#Figures 6
#Figures 6.b precipitation × light

gge_precip <- ggpredict(sel_mod, terms = c("precipitation [all]"))

p_precip <- ggplot(gge_precip, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Precipitation (in mm)",
    y = NULL,
    color = "Light treatment",
    fill = "Light treatment"
  ) +
  theme_minimal()+ 
  coord_cartesian(ylim = c(0, 80)) +
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

#Figures 6.a temperature × light
gge_temp <- ggpredict(sel_mod, terms = c("temperature [17.7:21.1 by=0.1]", "light"))

p_temp <- ggplot(gge_temp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Temperature (in °C)",
    y = "Height growth (cm.year⁻¹)",
    color = "Light treatment",
    fill = "Light treatment"
  ) +
  theme_minimal()+ 
  coord_cartesian(ylim = c(0, 80)) +
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

#Figures 6.c late frost events × light
gge_frost <- ggpredict(sel_mod, terms = c("nb_late_frost_event [all]", "light"))

p_frost <- ggplot(gge_frost, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "# late frost events",
    y = NULL,
    color = "Light treatment",
    fill = "Light treatment"
  ) +
  coord_cartesian(ylim = c(0, 80)) +
  theme_minimal()+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

quantile(growth$precipitation, prob=c(0, 0.1, 0.5, 0.9, 1))

min(growth$temperature)
max(growth$temperature)
gge_prec_temp <- ggpredict(
  sel_mod,
  terms = c(
    "temperature [17.7:21.1 by=0.1]",
    "precipitation [316, 330, 383]"
  )
)

p_temp_prec <- ggplot(gge_prec_temp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Temperature",
    y = NULL) +
  coord_cartesian(ylim = c(0, 80)) +
  scale_color_viridis_d(name = "Precipitation (in mm)") +
  scale_fill_viridis_d(name = "Precipitation (in mm)") +
  theme_minimal()+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))+
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  )

clim_height <- ggarrange(
  p_temp, p_precip, p_frost,
  ncol = 3, nrow = 1,
  common.legend = TRUE, legend = "bottom",
  labels = c("a","b","c"),
  label.x = c(0, -0.04, -0.04),  
  label.y = c(1, 1, 1)
) +
  bgcolor("white") +
  border("white")

clim_height2 <- ggarrange(clim_height, ggarrange(p_temp_prec, labels=c("d"),
                                                 label.x = -0.04,  
                                                 label.y =  1), ncol=2, nrow=1, widths=c(0.75,0.25))

ggsave(plot=clim_height2, filename="figures/growth_clim_effect.png", 
       width=11, height=4)


#Figure 7 Height growth differrences between sites
growth$fitted <-  unname(exp(predict(sel_mod, newdata = growth, type = "response")) - 1)

site_effect <- growth %>% 
  group_by(site, light) %>% 
  summarize(mean = mean(fitted),
            sd = sd(fitted)) %>% 
  rename(`Light level` = light)

site_effect$site <- factor(site_effect$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("within or", "north of the current sugar maple range")) %>% 
  rename(`Site`=labs)

rect_df$`Site` <- factor(rect_df$`Site`, levels = c("within or", "north of the current sugar maple range"))

growth_obs <- growth %>% 
  group_by(site, light) %>% 
  summarize(med_growth=median(growth)) %>% 
  rename(`Light level`=light)

growth_obs$site <- factor(growth_obs$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

site_plot <- ggplot()+
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE)+
  scale_fill_viridis_d() +
  geom_point(data=site_effect, aes(x=site, y= mean)) +
  geom_errorbar(data=site_effect, aes(x=site, ymin=mean-sd, ymax=mean+sd)) +
  geom_point(data=growth_obs, aes(x=site, y= med_growth), color="red", shape=3, size=3) +
  facet_wrap(~`Light level`) +
  ylab("Height growth (cm.year⁻¹)") + 
  xlab("Sites") + 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

site_plot

ggsave(plot=site_plot, filename="figures/growth_site_effect.png", 
       width=8, height=4)

