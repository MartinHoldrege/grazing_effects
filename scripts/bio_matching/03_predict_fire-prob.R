# Purpose: Predict fire probabiliy, using long term means from stepwat
# and see how closely it actually matches fire probability simulated
# in stepwat (to see if using the equation to predict fire
# probability this way is a good enough approximation
# of how changing (biomass) inputs will change fire probability in stepwat)

# Author: Martin Holdrege

# Started: Nov 8, 2024


# dependencies ------------------------------------------------------------

library(tidyverse)
source('src/fire_functions.R')
source('src/general_functions.R')
source('src/fig_params.R')
source('src/fig_functions.R')
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

# created in scripts/02_summarize_bio.R
bio <- readRDS('data_processed/site_means/summarize_bio.RDS')

path <- 'data_processed/temp_rds/qm_funs.rds'
if(file.exists(path)) {
  qm_l <- readRDS(path)
} else {
  source('scripts/bio_matching/02_quantile_matching.R')
  qm_l <- readRDS(path)
}

# params ------------------------------------------------------------------

run <- "fire1_eind1_c4grass0_co20_2311"
dataset <- qm_l$dataset

# functions ---------------------------------------------------------------

# functions for matching from stepwat biomass to rap biomass
# (via quantile matching)
qm_Aherb <- qm_l$Aherb
qm_Pherb <- qm_l$Pherb

# now predict fire probability, but after stepwat biomass 
# has been put on the scale of RAP biomass
predict_fire_qm <- function(mat, map, psp, afg, pfg) {
  afg_adj  <- qm_Aherb(afg)
  pfg_adj <- qm_Pherb(pfg)
  
  predict_fire(mat, map, psp, afg_adj, pfg_adj)
}

# prepare dataframes ------------------------------------------------------
clim0 <- bio$clim_all2
pft5_bio2 <- bio$pft5_bio2
fire_med1 <- bio$fire_med1

fire_med2 <- fire_med1 %>% 
  select(-n_fires, -fire_return)

# climate

clim_med <- clim0 %>% 
  group_by(site, years, RCP) %>% 
  summarize(across(where(is.numeric), .fns = median),
            .groups = 'drop')

# wide 
bio_w1 <- pft5_bio2 %>% 
  filter(run == !!run) %>% 
  calc_aherb(group_cols = c('run', 'years', 'RCP', 'graze', 'id', 
                            'site')) %>%  # add a total annual herbacious category
  filter(PFT %in% c('Pherb', 'Aherb')) %>% 
  select(run,id, site, graze, biomass, PFT, RCP, years) %>% 
  pivot_wider(names_from = 'PFT',
              values_from = 'biomass')

bio_w2 <- bio_w1 %>% 
  left_join(clim_med, by = c('site', 'RCP', 'years'))

comb1 <- bio_w2 %>% 
  left_join(fire_med2) %>% 
  rename(sim_fire_prob = fire_prob)
  
# predict fire probability --------------------------------------------------

comb2 <- comb1 %>% 
  mutate(pred_fire_prob_raw = predict_fire(mat = MAT,
                                      map = MAP,
                                      psp = psp,
                                      afg = Aherb,
                                      pfg = Pherb)*100, # converting to %
         pred_fire_prob_qm = predict_fire_qm(mat = MAT,
                                           map = MAP,
                                           psp = psp,
                                           afg = Aherb,
                                           pfg = Pherb)*100 # converting to %
         )

# * grazing effect ----------------------------------------------------------
# the simulated and predicted effects of grazing on fire probability

tmp1 <- comb2 %>% 
  select(run, site, graze, RCP, years, matches('fire')) %>% 
  pivot_longer(cols = matches('fire'))

# difference in fire probability between light grazing
# and another level of grazing
delta_fire_grazeL <- tmp1 %>% 
  filter(graze == 'Light') %>% 
  select(-graze) %>% 
  rename(value_light_graze = value) %>% 
  right_join(tmp1) %>% 
  mutate(delta_fire_prob = value - value_light_graze) %>% 
  filter(graze != 'Light') %>% 
  select(-value, -value_light_graze) %>% 
  pivot_wider(values_from = 'delta_fire_prob',
              names_from = 'name')

# figure params ---------------------------------------------------------------

# axis labels
lab_pred_qm <- 'Predicted w/ quantile matched STEPWAT2 biomass'
lab_pred_raw <- 'Predicted w/ raw STEPWAT2 biomass'
lab_sim <- 'Simulated by STEPWAT2'
# captions
cap_qm <- paste('Quantile matching done between', dataset, 'STEPWAT2 and RAP')

cap1 <- paste0('Fire equation used to predict fire probability using raw', 
               ' (not quantile matched) STEPWAT biomass.',
               '\n', run)
cap2 <- paste0('Fire equation used to predict fire probability using both raw',
               ' and quantile matched STEPWAT2 biomass', 
               '\n', cap_qm, '\n', run)
# figures--simulated vs raw predicted -----------------------------------------



pdf('figures/bio_matching/fire-prob_predicted_vs_simulated.pdf')
ggplot(comb2, aes(sim_fire_prob, pred_fire_prob_raw)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(subtitle = 'Fire probability (%)', 
       x = lab_sim,
       y = lab_pred_raw,
       caption = cap1) +
  facet_wrap(RCP~years)

g1 <- ggplot(delta_fire_grazeL, aes(sim_fire_prob, pred_fire_prob_raw, 
                                    color = graze)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(subtitle = 'Change in fire probability (# fires/100 yrs)  relative to light grazing',
       x = lab_sim,
       y = lab_pred_raw,
       caption = cap1) +
  scale_color_manual(values = cols_graze)

g1 + facet_wrap(RCP~years)

g1 + 
  facet_wrap(~graze, ncol = 2) +
  theme(legend.position = 'none') +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5)
  
dev.off()  

# figures--qm vs raw predicted -----------------------------------------


pdf('figures/bio_matching/fire-prob_qm_vs_raw.pdf')

ggplot(comb2, aes(pred_fire_prob_raw, pred_fire_prob_qm)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(subtitle = 'Fire probability (%)', 
       x = lab_pred_raw,
       y = lab_pred_qm,
       caption = cap2) +
  facet_wrap(RCP~years)

g1 <- ggplot(delta_fire_grazeL, aes(pred_fire_prob_raw, pred_fire_prob_qm, 
                                    color = graze)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(subtitle = 'Change in fire probability (# fires/100 yrs)  relative to light grazing',
       x = lab_pred_raw,
       y = lab_pred_qm,
       caption = cap2) +
  scale_color_manual(values = cols_graze)

g1 + facet_wrap(RCP~years)

g1 + 
  facet_wrap(~graze, ncol = 2) +
  theme(legend.position = 'none') +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.5)

dev.off()  

