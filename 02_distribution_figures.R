# Purpose: figures comparing RAP and stepwat biomass distributions

# Author: Martin Holdrege


# dependencies ------------------------------------------------------------


# climate figures -------------------------------------------------------------

g <- ggplot(clim2, aes(color = dataset_overlap)) +
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

pdf(paste0('figures/RAP/RAP-vs-sw_hists_', run, '.pdf'),
    width = 8)

# interpolated data comparison
ggplot() +
  geom_histogram(data = rap_df_m, 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = df_sw2, 
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
  geom_histogram(data = rap_sub2, 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = sw_site_bio1, 
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
         'Stepwat data from ', length(unique(site_ids)), ' sites',  
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

