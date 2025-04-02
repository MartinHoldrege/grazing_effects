# Purpose: figures of site level fire probability, and change in probabilities
# etc.

# Some of this code was formerly in the 03_bio_figs.R script

# Started: March 14, 2025

# Author Martin Holdrege


# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)

path_rds <- 'data_processed/site_means/summarize_bio.RDS'

if(!file.exists(path_rds)) {
  source("scripts/02_summarize_bio.R") # where needed dataframes are created
}

source("src/fig_params.R")
source("src/fig_functions.R")
source('src/general_functions.R')

# params ------------------------------------------------------------------
# misc. parameters and settings for figure creation

theme_set(theme_custom1())
ref_graze <- 'Light'

# read in data ------------------------------------------------------------

l <- readRDS(path_rds) # created in 02_sommarize_bio.R

fire_med1 <- l$fire_med1
fire_d_wgraze <- l$fire_d_wgraze
clim1 <- l$clim1
fire_d_wgcm <- l$fire_d_wgcm
runs_graze <- l$runs_graze

stopifnot(
  # check that correct reference grazing level is assumed
  ! ref_graze %in% fire_d_wgcm$graze
)

# prepare df's ------------------------------------------------------------

fire_med2 <- fire_med1 %>% 
  left_join(clim1, by = 'site')


# boxplot of fire return interval by RCP, and grazing intensity

map(runs_graze, function(x) {
  df <- fire_med2 %>% 
    filter(run == x) 
  g1 <- ggplot(df, aes(x = id, y = fire_prob, fill = RCP)) +
    geom_text(data = ~box_anno(., var = "fire_prob",
                               group_by = c("graze")),
              aes(x, y, label = graze, fill = NULL),
              size = 2.5) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_rcp, name = "Scenario") +
    scale_x_discrete(labels = id2year) +
    geom_vline(xintercept = line_loc, linetype = 2) +
    theme(legend.position = "top",
          axis.text = element_text(size = 7))  +
    labs(x = lab_yrs,
         y = lab_firep0) 
  
  g2 <- ggplot(df, aes(x = rcp_label(RCP, years, include_parenth = FALSE), 
                       y = fire_prob, fill = graze))+
    geom_boxplot() +
    scale_fill_graze() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    labs(y = lab_firep0,
         x = NULL)
  
  # fire difference boxplot
  g3 <- fire_d_wgraze %>% 
    filter(run == x) %>% 
    group_by(id) %>% 
    #remove_outliers(var = "fire_prob_diff") %>% 
    ggplot(aes(id2, fire_prob_diff, fill = graze)) +
    geom_text(data = ~box_anno(., var = "fire_prob_diff", 
                               group_by = c("RCP"), id = "id2"),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_boxplot(position = "dodge",
                 outlier.color = NA) +
    scale_fill_graze() +
    # so just display the RCP
    scale_x_discrete(labels = id2year) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 2.5, linetype = 2) +
    labs(x = lab_yrs,
         y = lab_firep1)
  
  jpeg(paste0("figures/fire/fire_return_boxplots_", x, ".jpeg"),
       height = 7, width = 7, res = 600, units = "in")
  print(g1 + g2 + g3 + guide_area() + patchwork::plot_layout(guides = 'collect'))
  dev.off()
})


# * within gcm grazing effect -----------------------------------------------
# compare fire probability in light vs very  heavy grazing
fire_med3 <- fire_med1 %>% 
  filter(run == runs_graze['default'])

fire_med_ref0 <- fire_med3 %>% 
  filter(graze == ref_graze) %>% 
  rename(fire_prob_ref = fire_prob) %>% 
  select(-graze, -id, -n_fires, -fire_return) %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE))

# compare reference (light) grazing prob vs delta prob when going to higher graze
fire_d_wgcm2 <-  fire_d_wgcm %>% 
  filter(run == runs_graze['default']) %>% 
  left_join(fire_med_ref0)


# compare ref grazing prob to other grazing levels
fire_med_ref1 <- fire_med3 %>% 
  filter(graze != ref_graze) %>% 
  left_join(fire_med_ref0)
fire_med_ref1 %>% 
  mutate(diff = fire_prob - fire_prob_ref) %>% 
  pull(diff) %>% 
  summary()

pdf(paste0("figures/fire/compare_graz_levs_v1_", runs_graze['default'],".pdf"),
    width = 8, height = 7)

ggplot(fire_d_wgcm2, aes(fire_prob_ref, fire_prob_diff,
                         color = graze)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'loess', se = FALSE) +
  facet_wrap(~rcp_year) +
  scale_color_graze() +
  labs(x = paste0('Fires (#/century) under ', ref_graze, ' grazing'),
       y = lab_firep1,
       subtitle = 'Change in fire probability with change in grazing intensity')

ggplot(fire_med_ref1, aes(fire_prob_ref, fire_prob,
                         color = graze)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  facet_wrap(~rcp_year) +
  geom_abline(slope = 1) +
  scale_color_graze() +
  labs(x = paste0('Fires (#/century) under ', ref_graze, ' grazing'),
       y = paste0('Fires (#/century) under elevated grazing intensity'),
       subtitle = 'Comparing fire probabilties')

dev.off()

# * fire vs climate -------------------------------------------------------

rn <- runs_graze['default']
pdf(paste0("figures/fire/fire_vs_clim_", rn, ".pdf"),
    width = 8, height = 7)

g1 <- fire_med2 %>% 
  filter(run == rn) %>% 
  ggplot(aes(y = fire_prob, color = graze)) +
  lemon::facet_rep_wrap(~RCP + years) + 
  scale_color_graze() +
  labs(y = lab_firep0, 
       caption = paste(rn, '\n x-axis shows current climate in all panels')) +
  theme(legend.position = c(0.85, 0.15),
        axis.text = element_text(size = 7))

climate_scatter(g1, include_psp = TRUE)
dev.off() 
