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

# functions for quantile matching
path <- 'data_processed/temp_rds/qm_funs.rds'
if(file.exists(path)) {
  qm_l <- readRDS(path)
} else {
  source('scripts/bio_matching/02_quantile_matching.R')
  qm_l <- readRDS(path)
}

# params ------------------------------------------------------------------

run <- qm_l$run
run2 <- qm_l$run2
RCP <- 'Current'
stopifnot(run == unique(m$sw_comb$run)) # there may be cases where this doesn't need to be true
qual_cutoff <- qm_l$qual_cutoff

# dataframes etc --------------------------------------------------------------

sw_comb <- m$sw_comb
rap_comb <- m$rap_comb

# *prepare pdp predictions -------------------------------------------------
pft_lookup <- c('pfgAGB' = 'Pherb','afgAGB' = 'Aherb')
pfts <- pft_lookup
names(pfts) <- pfts
df_pdp3 <- df_pdp2 %>% 
  filter(inter_var == 'Mean prediction',
         variable %in% c("pfgAGB", 'afgAGB')) %>% 
  mutate(PFT = pft_lookup[variable]) %>% 
  select(-variable)

# calculate quantile matched biomass ------------------------------------------

sw_comb$biomass_qm <- NA

for (pft in pfts) {
  f <- qm_l[[pft]]
  sw_comb$biomass_qm[sw_comb$PFT ==pft] <- f(sw_comb$biomass[sw_comb$PFT ==pft])
}

# Pherb:total herb
ratio <- sw_comb %>% 
  filter(dataset != 'site level') %>% 
  pivot_wider(id_cols = c('run', 'fire', 'graze', 'dataset', 'site', "cell_num"),
              values_from = c('biomass', 'biomass_qm'),
              names_from = 'PFT') %>% 
  mutate(P_ratio_qm = biomass_qm_Pherb/(biomass_qm_Aherb + biomass_qm_Pherb),
         P_ratio_raw = biomass_Pherb/(biomass_Aherb + biomass_Pherb))


# fig params --------------------------------------------------------------



# climate figures -------------------------------------------------------------

# need to use December 2024 version of repo for these figures

# density figures ---------------------------------------------------------

# captions
cap_match <- paste0('\nMatching critera from Palmquist used.',
                    '\nFor given pixels RAP data includes 3-year averages for all years in Holdrege et al.')
cap_interp <- paste0('\nData from 1km gridcells where interpolation matching quality < ', 
                     qual_cutoff)

cap_qm <- qm_l$caption

cap_sites <- paste0(
  'Stepwat data from ', length(unique(rap_comb$site)), ' sites',  
  "\n(other sites didn't overlap with Holdrege et al. study area).",
  '\nRAP data from Holdrege et al. only from grid-cells in those site locations',
  "Black line is the mean fire model prediction")

pdp_line <- function() {
  geom_line(data = df_pdp3, 
            aes(x = x_value, y = yhat*2))
}

pdf(paste0('figures/bio_matching/RAP-vs-sw_hists_', run, '_qm', qual_cutoff, 
           '.pdf'), width = 8)

# interpolated data comparison
g1 <- ggplot() +
  geom_histogram(data = rap_comb[rap_comb$dataset == 'interpolated',], 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  pdp_line() +
  scale_color_manual(values = cols_graze,
                     name = 'STEPWAT2 grazing') +
  facet_wrap(~PFT, ncol = 1, scales = 'free') +
  scale_y_continuous(sec.axis = sec_axis(~./2, name = 'Mean predicted fire probability')) +
  expand_limits(x = 0) +
  theme(legend.position = 'top') +
  labs(x = lab_bio0)

print(g1 +
  geom_density(data = sw_comb[sw_comb$dataset == 'interpolated', ], 
               aes(biomass, color = graze)) +
  labs(subtitle = 'Interpolated  stepwat data',
       caption = paste('Histogram shows RAP data, colored lines show stepwat (', 
                       RCP, ' climate), black is fire model',
                       cap_interp, cap_match)))

print(g1 +
  geom_density(data = sw_comb[sw_comb$dataset == 'interpolated', ], 
               aes(biomass_qm, color = graze)) +
  labs( subtitle = 'Interpolated stepwat (quantile mapped) data',
       caption = paste('Histogram shows RAP data, colored lines show stepwat quantile mapped values\n',
                       cap_qm,
                       cap_interp, cap_match)))

# stepwat site level comparison
# see 12/2024 version of repo



dev.off()

# Pherb/total herb ratio -------------------------------------------------------

# how does the perenial fraction change with quantile matching?
set.seed(1234)
g <- ratio %>% 
  slice_sample(n = 1000) %>% 
  filter(run == !!run) %>% 
  ggplot(aes(P_ratio_raw, P_ratio_qm)) +
  geom_point(aes(color = graze)) +
  scale_color_manual(values = cols_graze) +
  geom_abline(color = 'darkgray') +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(0, 1)) +
  labs(x = 'Pherb/(Pherb + Aherb) [Raw stepwat biomass]',
       y = 'Pherb/(Pherb + Aherb) [quantile mapped biomass]',
       subtitle = 'Proportion perennials before and after quantile mapping',
       caption = paste('Current climate conditions.\n',
                       run2, "\n Sample of pixels shown"))

pdf(paste0('figures/bio_matching/Pherb_ratio_qm', qual_cutoff, "_", run2, '.pdf'))
print(g)
dev.off()
