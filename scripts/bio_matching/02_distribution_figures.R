# Purpose: figures comparing RAP and stepwat biomass distributions

# Author: Martin Holdrege


# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
source('src/fig_params.R')
source('src/fig_functions.R')
theme_set(theme_custom1())

# read in data -----------------------------------------------------------
# list created in 01_combine_data.R
m <- readRDS('data_processed/temp_rds/rap_sw_matching.rds') 

# partial dependence plot values from fire equation
# file created in cheatgrass_fire/scripts/06_figures_pdp_vip_quant.R
df_pdp2 <- readRDS("../cheatgrass_fire/data_processed/df_pdp2.rds");

# params ------------------------------------------------------------------
run <- m$runs[1]
RCP <- 'Current'

# dataframes etc --------------------------------------------------------------

clim <- m$clim
sw_comb <- m$sw_comb
rap_comb <- m$rap_comb
qual_cutoff <- m$qual_cutoff

# *prepare pdp predictions -------------------------------------------------
pft_lookup <- c('pfgAGB' = 'Pherb','afgAGB' = 'Aherb')

df_pdp3 <- df_pdp2 %>% 
  filter(inter_var == 'Mean prediction',
         variable %in% c("pfgAGB", 'afgAGB')) %>% 
  mutate(PFT = pft_lookup[variable]) %>% 
  select(-variable)

# climate figures -------------------------------------------------------------

g <- ggplot(clim, aes(color = dataset_overlap)) +
  theme_bw() +
  scale_color_manual(name = 'STEPWAT site falls on \nHoldrege et al 2024 1km cell',
                     values = c('darkgrey', 'blue' ))

g1 <- g + geom_point(aes(MAT_C, MAP_mm))

g2 <- g + geom_point(aes(CorrTP2, MAP_mm)) +
  labs(x = 'T-P correlation')
g3 <- g + geom_point(aes(CorrTP2, MAT_C)) +
  labs(x = 'T-P correlation')
p <- g1 + g2 + g3 + plot_layout(guides = 'collect',
                                ncol = 2) 
p&theme(legend.position = 'top')

linetype <- function() {
  scale_linetype_manual(name = NULL, 
                        values = c("fire1" = 1, "fire0" = 2))
}

# density figures ---------------------------------------------------------

pdf(paste0('figures/bio_matching/RAP-vs-sw_hists_', run, '.pdf'),
    width = 8)

# interpolated data comparison
ggplot() +
  geom_histogram(data = rap_comb[rap_comb$dataset == 'interpolated',], 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = sw_comb[sw_comb$dataset == 'interpolated', ], 
               aes(biomass, color = graze, linetype = fire)) +
  geom_line(data = df_pdp3, 
            aes(x = x_value, y = yhat*2)) +
  scale_color_manual(values = cols_graze,
                     name = 'STEPWAT2 grazing') +
  facet_wrap(~PFT, ncol = 1, scales = 'free') +
  scale_y_continuous(sec.axis = sec_axis(~./2, name = 'Mean predicted fire probability')) +
  expand_limits(x = 0) +
  labs(x = lab_bio0,
       subtitle = 'Interpolated stepwat data',
       caption = paste0(
         'Histogram shows RAP data, colored lines show stepwat (', RCP, ' climate), black is fire model',
         '\nData from 1km gridcells where interpolation matching quality < ', 
         qual_cutoff,
         '\nMatching critera from Palmquist used.',
         '\nFor given pixels RAP data includes 3-year averages for all years in Holdrege et al.'
       )) +
  theme(legend.position = 'top') +
  linetype()

# stepwat site level comparison

ggplot() +
  geom_histogram(data = rap_comb[rap_comb$dataset != 'interpolated',], 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = sw_comb[sw_comb$dataset != 'interpolated', ], 
               aes(biomass, color = graze, linetype = fire)) +
  geom_line(data = df_pdp3,
            aes(x = x_value, y = yhat*2)) +
  scale_color_manual(values = cols_graze,
                     name = 'STEPWAT2 grazing') +
  facet_wrap(~PFT, ncol = 1, scales = 'free') +
  scale_y_continuous(sec.axis = sec_axis(~. /2, name = 'Mean predicted fire probability')) +
  expand_limits(x = 0) +
  labs(x = lab_bio0,
       subtitle = 'Stepwat site level data (not interpolated)',
       caption = paste0(
         'Stepwat data from ', length(unique(rap_comb$site)), ' sites',  
         "\n(other sites didn't overlap with Holdrege et al. study area).",
         '\nRAP data from Holdrege et al. only from grid-cells in those site locations',
         "Black line is the mean fire model prediction")) +
  theme(legend.position = 'top') +
  linetype()

# histograms of just stepwat data 
# showing these seperately because density curves may be hiding details

for (pft in pft_lookup) {
  g <- sw_comb %>% 
    filter(PFT == pft) %>% 
    ggplot(aes(x = biomass, y = after_stat(density), fill = fire)) +
    geom_histogram() +
    facet_grid(dataset~graze) +
    labs(title = pft,
         subtitle = 'Comparing site level and interpolated stepwat biomass')
  print(g)
}

dev.off()

