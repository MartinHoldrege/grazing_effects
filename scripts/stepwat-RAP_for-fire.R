# Purpose: Compare stepwat and RAP biomass, at the 200 simulation sites 
# so that some sort of quantile matching relationship can determined
# so that stepwat biomass is translated to biomass levels used for 
# fitting the Holdrege et al. 2024 fire equation

# Author: Martin Holdrege

# Oct 3, 2024


# params ------------------------------------------------------------------

version_interp <- 'v3' # interpolation version (see 03_interpolate.R)
run <- "fire1_eind1_c4grass1_co20_2311"

# matching quality cutoff, this decides
# the cells from which both stepwat and rap data will be pulled from
qual_cutoff <- 0.5 
RCP <- 'Current'
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(patchwork)
source("src/general_functions.R")
source('src/fig_params.R')
source('src/fig_functions.R')
theme_set(theme_custom1())
# read in data ------------------------------------------------------------

# data from the Holdrege et al. 2024 publication
rap_df1 <- read_csv("../cheatgrass_fire/data_processed/data_publication/Wildfire_Climate_Biomass_Data.csv")

cell_nums <- rast("../cheatgrass_fire/data_processed/data_publication/cell_nums.tif")

# created in grazing_effects/scripts/02_summarize_bio.R
pft5_bio2 <- readRDS('data_processed/site_means/summarize_bio.RDS')$pft5_bio2

sites1 <- read_csv("data_raw/site_locations.csv")
clim1 <- read_csv("data_processed/site_means/dbWeather_200sites.csv")

# partial dependence plot values from fire equation
# file created in cheatgrass_fire/scripts/06_figures_pdp_vip_quant.R
df_pdp2 <- readRDS("../cheatgrass_fire/data_processed/df_pdp2.rds");

# * interpolated stepwat --------------------------------------------------
PFTs <- c('Aherb', 'Pherb')

paths <- list.files("data_processed/interpolated_rasters/biomass/",
                    pattern = paste0(run, version_interp, ".*", RCP, ",*.tif"),
                    full.names = TRUE)

r_sw1 <- rast(paths)

r_qual <- rast(paste0("data_processed/interpolation_quality/matching_quality",
                 version_interp, ".tif"))


# vectors etc -------------------------------------------------------------

pft_lookup <- c('pfgAGB' = 'Pherb','afgAGB' = 'Aherb')

# prepare pdp predictions -------------------------------------------------

df_pdp3 <- df_pdp2 %>% 
  filter(inter_var == 'Mean prediction',
         variable %in% c("pfgAGB", 'afgAGB')) %>% 
  mutate(PFT = pft_lookup[variable]) %>% 
  select(-variable)

# stepwat site level analysis ---------------------------------------------

# * select gridcells of rap data  ------------------------------------------

sites2 <- sf::st_as_sf(sites1, coords = c("X_WGS84", "Y_WGS84"),
                       crs = 'EPSG:4326') %>% 
  st_transform(crs = crs(cell_nums)) %>% 
  arrange(site_id)

plot(cell_nums)
plot(sites2$geometry, add = TRUE)

# checking that ID returned from extract will match site num
stopifnot(sites2$site_id == 1:max(sites2$site_id)) 
cell_nums_sites <- extract(cell_nums, sites2[, ]) %>% 
  rename(site = ID) %>% 
  as_tibble()

site_ids <- cell_nums_sites$site[!is.na(cell_nums_sites$cell_num)]

# filter site level data

sw_site_bio1 <- pft5_bio2 %>% 
  calc_aherb(group_cols = c('run', 'years', 'RCP', 'graze', 'id', 
                            'site')) %>% 
  filter(.data$run == !!run,
         RCP == "Current",
         site %in% site_ids,
         PFT %in% unname(pft_lookup))

# * evaluate clim of stepwat sites ------------------------------------------
# see if the subset of stepwat sites, that overlap with the rap data from the
# fire data are representative of all sites

clim2 <- clim1 %>% 
  mutate(dataset_overlap = Site_id %in% site_ids)

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

# ** extract RAP data --------------------------------------------------------

rap_df2 <- rap_df1 %>% 
  select(year, cell_num, pfgAGB, afgAGB)

rap_sub1 <- cell_nums_sites %>% 
  filter(.data$site %in% site_ids) %>% 
  left_join(rap_df2, by = 'cell_num')

hist(rap_sub1$pfgAGB, xlim = c(0, 200), breaks = 200)

rap_sub2 <- rap_sub1 %>%  
  pivot_longer(cols = c('pfgAGB', 'afgAGB'),
               values_to = 'biomass',
               names_to = 'PFT') %>% 
  mutate(PFT = pft_lookup[PFT])

# pulling together interpolation level data -------------------------------

# identify stepwat grid-cells w/ good matching, and that are also cells
# included in the fire manuscript
cell_nums2 <- extend(cell_nums, r_qual)
mask0 <- r_qual < qual_cutoff & !is.na(cell_nums2)
mask <- classify(mask0, rcl = cbind(FALSE, NA))

r_sw2 <- mask(r_sw1, mask)
df_pherb <- r_sw2[[str_subset(names(r_sw2), 'Pherb')]]
df_aherb <- r_sw2
cell_nums_m <- mask(cell_nums2, mask)
cell_nums_mv <- as.numeric(cell_nums_m[!is.na(cell_nums_m)]) # vector

rap_df_m <- rap_df2 %>% 
  filter(cell_num %in% cell_nums_mv) %>% 
  pivot_longer(cols = c("pfgAGB", 'afgAGB'),
               names_to = 'PFT',
               values_to = 'biomass') %>% 
  mutate(PFT = pft_lookup[PFT]) %>% 
  select(-cell_num, -year)

df_sw1 <- as_tibble(r_sw2) 

info1 <- create_rast_info(r_sw2)

df_sw2 <- df_sw1 %>% 
  pivot_longer(cols = everything(),
               names_to = 'id',
               values_to = 'biomass') %>% 
  left_join(info1, by = 'id')

tmp1 <- df_sw2 %>% 
  select(biomass, PFT, graze) %>% 
  mutate(dataset = 'interpolated') 

sw_comb <- sw_site_bio1 %>% 
  select(biomass, PFT, graze) %>% 
  mutate(dataset = 'site level') %>% 
  bind_rows(tmp1)

# density figures ---------------------------------------------------------

pdf(paste0('figures/RAP/RAP-vs-sw_hists_', run, '.pdf'),
    width = 8)

# interpolated data comparison
ggplot() +
  geom_histogram(data = rap_df_m, 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = df_sw2, 
               aes(biomass, color = graze)) +
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
  theme(legend.position = 'top')

# stepwat site level comparison

ggplot() +
  geom_histogram(data = rap_sub2, 
                 aes(biomass, y = after_stat(density)),
                 bins = 100) + 
  geom_density(data = sw_site_bio1, 
               aes(biomass, color = graze)) +
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
  theme(legend.position = 'top')

# histograms of just stepwat data 
# showing these seperately because density curves may be hiding details

for (pft in pft_lookup) {
  g <- sw_comb %>% 
    filter(PFT == pft) %>% 
    ggplot(aes(x = biomass, y = after_stat(density))) +
    geom_histogram() +
    facet_grid(dataset~graze) +
    labs(title = pft,
         subtitle = 'Comparing site level and interpolated stepwat biomass')
  print(g)
}



dev.off()

