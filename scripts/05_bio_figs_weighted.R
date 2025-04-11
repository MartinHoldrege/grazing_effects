# Purpose: TBD

# Author: Martin Holdrege

# Started: April 11, 2025

# params ------------------------------------------------------------------

v_interp <- 'v4' # interpolation version
run <- "fire1_eind1_c4grass1_co20_2503"
# dependencies ------------------------------------------------------------

library(tidyverse)
theme_set(theme_custom1())
library(patchwork)
source("src/fig_params.R")
source("src/fig_functions.R")
source('src/general_functions.R')

# read in data ------------------------------------------------------------

path_rds <- 'data_processed/site_means/summarize_bio.RDS'

if(!file.exists(path_rds)) {
  source("scripts/02_summarize_bio.R") # where needed dataframes are created
}
l <- readRDS(path_rds) # created in 02_sommarize_bio.R

# weight (based on area of interpolation) corresponding to each site
# in each ecoregion
# created in 04_interpolation_weights.R
w1 <- read_csv(paste0('data_processed/interpolation_data/interpolation_weights_', 
                       v_interp, '.csv'),
               show_col_types = FALSE)

bio2 <- l$pft5_bio2

# combine site level and weights ------------------------------------------

bio3 <- bio2 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  filter_clim_extremes() %>% 
  filter(run == !!run) %>% 
  df_factor() 


# longer format
bio4 <- bio3 %>% 
  select(-matches('utilization'), -matches('indivs')) %>% 
  pivot_longer(matches('biomass'),
               names_to = 'summary',
               values_to = 'biomass') %>% 
  mutate(summary = str_extract(summary, '(?<=_).*'),
         summary = ifelse(is.na(summary), "median", summary),
         summary = factor(summary, levels = c('low', 'median', 'high')),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE))


# boxplots ----------------------------------------------------------------

pfts <- c('Sagebrush', 'Pherb', 'Aherb')
# continue here

plots <- map(pfts, function(pft) {
  bio4 %>% 
    filter(PFT == pft) %>% 
    ggplot(aes(graze, biomass, fill = summary)) +
    geom_boxplot(aes(weight = weight), position = position_dodge2(preserve = 'single'),
                 outlier.size = 0.25, outlier.alpha = 0.5) + 
    facet_grid(region~rcp_year) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x.top = element_text(size = rel(0.7))) +
    labs(x = lab_graze,
         y = lab_bio0,
         subtitle = pft,
         fill = 'Summary across GCMs')
})


g <- wrap_plots(plots, nrow = 1) +
  plot_layout(guides = 'collect', axis_titles = 'collect')
g2 <- g&theme(legend.position = 'bottom')

png(paste0("figures/biomass/bio_weighted_by-region_3pft_boxplot_", run, ".png"),
    width = 8, height = 8, units = 'in', res = 600)
g2
dev.off()
