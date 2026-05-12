#Script taht model the germination rate of sugar maple. 

rm(list=ls())

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(readxl)
library(lme4)
library(lmerTest)
library(car)
library(AICcmodavg)
library(viridis)
library(sjPlot)
library(MuMIn)
library(ggeffects)
library(DHARMa)
library(emmeans)

site <- data.frame(site=c("MUS", "HED", "MON", "KEK", "KIP"), 
                   latitude=c(49.932, 49.243, 48.460, 48.191, 46.740),
                   longitude=c(-78.698, -78.311, -79.418, -79.112, -78.905))

germ <- read_excel("data/germination.xlsx") %>% 
  mutate(year=as.factor(year),
         light=as.factor(light))
head(germ)

germ$site <- factor(germ$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Climate
clim <- readRDS("data/clim_germination.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

#Join climate and germination
germ <- germ %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(site_block=paste0(site, "_", block))

#Mean germiantion rate accroding to light
germ %>%
  group_by(light) %>%
  summarize(
    mean = mean(germination, na.rm = TRUE),
    sd   = sd(germination, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),  # erreur standard
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se
  )

germ$site <- relevel(as.factor(germ$site), ref = "KIP")

germ$light <- factor(germ$light, levels = c("10", "30", "100"), 
                          labels = c("Dense-canopy", "Intermediate", "Open-canopy"))

germ <- germ %>% left_join(site)

#Random effect selection model
mod_re1 <- glmer(cbind(germination,100-germination) ~ latitude + light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) + 
                   (1 | site_block) + (1 | year), family = binomial, data=germ)

mod_re2 <- glmer(cbind(germination,100-germination)~ latitude + light  + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                   (1 | site_block), family = binomial, data=germ)

mod_re3 <- glm(cbind(germination,100-germination)~ latitude + light + scale(temperature) + 
                 scale(nb_late_frost_event) + scale(precipitation), family = binomial, data=germ)

AIC_re <- data.frame(mod=1:3, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3)))
AIC_re

#Fixed effect selection
mod_me1 <- glmer(cbind(germination,100-germination)~ latitude + light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation)  +
                   scale(temperature):scale(precipitation) +
                  light:scale(temperature) + light:scale(precipitation) + light:scale(nb_late_frost_event)+
                   (1 | site_block), 
                 family = binomial, data=germ)

mod_me2 <- glmer(cbind(germination,100-germination)~ latitude +light +scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                   scale(temperature):scale(precipitation) +
                  light:scale(temperature) + light:scale(precipitation)+
                 (1 | site_block), family = binomial, data=germ)

mod_me3 <- glmer(cbind(germination,100-germination)~ latitude +light + scale(temperature)+ 
                  scale(nb_late_frost_event) + scale(precipitation)+
                   scale(temperature):scale(precipitation) +
                  light:scale(temperature)+
                 (1 | site_block), family = binomial, data=germ)

mod_me4 <- glmer(cbind(germination,100-germination)~ latitude +light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                   scale(temperature):scale(precipitation)+
                 (1 | site_block), family = binomial, data=germ)

mod_me5 <- glmer(cbind(germination,100-germination)~ latitude +light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation)+
                 (1 | site_block), family = binomial, data=germ)

mod_me6 <- glmer(cbind(germination,100-germination)~ latitude+ light + scale(temperature) + 
                   scale(precipitation)+
                 (1 | site_block), family = binomial, data=germ)

mod_me7 <- glmer(cbind(germination,100-germination)~ latitude +light + scale(temperature)+
                 (1 | site_block), family = binomial, data=germ)

mod_me8 <- glmer(cbind(germination,100-germination)~ latitude +light+
                 (1 | site_block), family = binomial, data=germ)

mod_me9 <- glmer(cbind(germination,100-germination)~ latitude+
                 (1 | site_block), family = binomial, data=germ)

mod_me10 <- glmer(cbind(germination,100-germination)~ 1+
                  (1 | site_block), family = binomial, data=germ)

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
  mod_me10 = mod_me10)
aictab(liste_mod)

sel_mod <- mod_me1




emm_light <- emmeans(
  sel_mod,
  ~ light,
  type = "response"
)

source("diagonistic_plot_function.R")

#plot diagnostics
ggsave(plot=diagnostic_plots(sel_mod, germ), filename="figures/SI_diag_plot_germ.png", 
       width=8, height=10)



summary(sel_mod)
car::Anova(sel_mod, type=2)
anova(sel_mod)

#Marginal and conditional R²
library(MuMIn)
r.squaredGLMM(sel_mod)


#Prédiction for each climate effect with interaction (Figure 3)
pred_temp  <- ggpredict(sel_mod, terms = c("temperature [all]", "light"),
                        back_transform = TRUE)

pred_prec  <- ggpredict(sel_mod, terms = c("precipitation [all]", "light"))

pred_frost <- ggpredict(sel_mod, terms = c("nb_late_frost_event [all]", "light"))


make_plot <- function(pred, xlab) {
  ggplot(pred, aes(x = x, y = predicted, color = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.2, color = NA) +
    labs(x = xlab, y = "Predicted germination probability", color = "Light treatment", fill = "Light treatment") +
    theme_minimal()+ 
    coord_cartesian(ylim = c(0, 0.49)) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5),
          strip.text.x = element_text(size = 14),
          axis.text=element_text(size=14),
          axis.title=element_text(size=15),
          legend.title = element_text(size = 16),
          legend.text = element_text(size=14),
          legend.position = "bottom") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))+
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
  
}

p1 <- make_plot(pred_temp,  "Temperature")
p2 <- make_plot(pred_prec,  "Precipitation (in mm)")+
  ylab(NULL)
p3 <- make_plot(pred_frost, "# late frost events")+
  ylab(NULL)

quantile(germ$precipitation, prob=c(0, 0.1, 0.5, 0.9, 1))

pred_prec_temp  <- ggpredict(sel_mod, terms = c("temperature [all]", "precipitation [118, 141, 180]"))
str(pred_prec_temp)

pred_prec_temp <- as.data.frame(pred_prec_temp) %>% 
  select(temperature=x, germination=predicted,
         std.error, conf.low, conf.high, precipitation=group)




prec_temp <- ggplot(pred_prec_temp,
       aes(x = temperature,
           y = germination,
           color = precipitation,
           group = precipitation)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = precipitation),
    alpha = 0.25,
    color = NA
  ) +
  
  scale_color_viridis_d(name = "Precipitation (in mm)") +
  scale_fill_viridis_d(name = "Precipitation (in mm)") +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted germination probability"
  ) +
  
  coord_cartesian(ylim = c(0, 0.49)) +
  theme_minimal()+
  theme(
    legend.position = "right"
  )+
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5),
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.title = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.position = "bottom")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
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





#Combine plots
plot_germ_effect <- ggarrange(p1, p2, p3, nrow=1, ncol=3, 
                              common.legend = T, legend="bottom",
                              labels=c("a","b", "c"),
                              label.x = c(-0.02),  # décalage spécifique pour chaque
                              label.y = c(1.02)) + 
  bgcolor("white")+
  border("white")


ggsave(plot=plot_germ_effect, filename="figures/germination_interaction_effect.png", 
       width=9, height=4.5)


#Prediction on the model by site (Figure 4)
germ$fitted <-  predict(sel_mod, newdata = germ, type = "response")

site_effect <- germ %>% 
  group_by(site, light) %>% 
  summarize(mean = mean(fitted),
            sd = sd(fitted)) %>% 
  rename(`Light level` = light)

site_effect$site <- factor(site_effect$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("within or", "north of the current sugar maple range")) %>% 
  rename(`Site`=labs)

rect_df$`Site` <- factor(rect_df$`Site`, levels = c("within or", "north of the current sugar maple range"))

germ_obs <- germ %>% 
  group_by(site, light) %>% 
  summarize(med_germ=median(germination)) %>% 
  rename(`Light level`=light)

germ_obs$site <- factor(germ_obs$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

site_plot <- ggplot()+
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE)+
  scale_fill_viridis_d() +
  geom_point(data=site_effect, aes(x=site, y= mean*100)) +
  geom_errorbar(data=site_effect, aes(x=site, ymin=mean*100-sd*100, ymax=mean*100+sd*100)) +
  geom_point(data=germ_obs, aes(x=site, y= med_germ), color="red", shape=3, size=3) +
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

ggsave(plot=site_plot, filename="figures/germination_site_effect.png", 
       width=8, height=4)

#Annual germination rate
germ_year <- germ %>% 
  group_by(year, site, light) %>% 
  summarize(germ=mean(germination),
            sd_germ=sd(germination))

plot_year <- ggplot(germ_year, aes(x=site, y=germ, color=light, ymin=germ-sd_germ, ymax=germ+sd_germ))+
  geom_errorbar(position = position_dodge(width = 0.6), linewidth=1)+
  facet_wrap(~year) +
  labs(
    x = "Site",
    y = "Observed germination probability",
    color = "Light treatment"
  ) + 
  theme(strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.title = element_text(size = 16),
        legend.text = element_text(size=14))

ggsave(plot=plot_year, filename="figures/SI_germ_year.png", 
       width=12, height=4.5)
  
########