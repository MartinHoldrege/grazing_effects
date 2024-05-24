# Purpose: 
# compare climate data from stepwat with climate data from gridMet/MACA, that
# is part of the new RR data release from DS
# for the 200 stepwat simulations, particularly for determinging why 
# the spatial patterns in the correlation between montly temp and precip
# seem different between the two data sources. If the data don't match well for the 200 
# sites that means it isn't an interpolation problem

# Author: Martin Holdrege

# script started: May 24, 2024


# dependencies ------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(patchwork)
source('src/paths.R')
theme_set(theme_bw())

# read in data ------------------------------------------------------------

loc1 <- read_csv('data_raw/site_locations.csv') # stepwat sites
sw_clim1 <- read_csv("data_processed/site_means/clim_seasonality.csv")

# data from the RR data release, this is the mean correlation between 
# monthly precip and temperature
rr_cur1 <- rast(file.path(path_large, 'newRR3-analysis/examples/tmp_examples',
                          'Preds19__sc1_NA_ambient_1980-2020-clim__value-sim.tif'))

rr_fut1 <- rast(file.path(path_large, 'newRR3-analysis/examples/tmp_examples',
                          'Preds19__RCP45_2064-2099-clim__value-acrmod-med.tif'))


# summarize across GCMs ---------------------------------------------------

sw_clim_med <- sw_clim1 %>% 
  filter(RCP != 'RCP85',
         years != '2030-2060') %>% 
  group_by(site, years, RCP) %>% 
  summarise(PTcor_stepwat = median(PTcor))
# extract clim at points --------------------------------------------------

# convert site locations to sf object
loc2 <- st_as_sf(loc1[, c('site_id', 'X_WGS84', 'Y_WGS84')],
                 crs = 'EPSG:4326',
                 coords = c('X_WGS84', 'Y_WGS84')) %>% 
  st_transform(crs = crs(rr_cur1)) %>% 
  arrange(site_id)

# testing ~~~~~~~~~~~~~~
minmax(rr_fut1)
print(rr_fut1)

names(rr_fut1)
plot(rr_fut1) # won't plot
# end testing ~~~~~~~~~~

terra::activeCat(rr_fut1) <- 'CorTempPPT_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'CorTempPPT_mean'

rr_fut2 <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(PTcor_rr = CorTempPPT_mean_acrmoddistr)
rr_cur2 <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(PTcor_rr = CorTempPPT_mean)

comb <- sw_clim_med %>% 
  left_join(rr_fut2, by = c('site' = 'ID', 'RCP'))%>% 
  left_join(rr_cur2, by = c('site' = 'ID', 'RCP')) %>% 
  mutate(PTcor_rr = ifelse(is.na(PTcor_rr.x), PTcor_rr.y, PTcor_rr.x)) %>% 
  select(-PTcor_rr.y, -PTcor_rr.x)


# figures -----------------------------------------------------------------

cap1 = 'The STEPWAT2 T-P correlation was calculated using montly avg T and P,\nthe GridMet/MACA correlation was calculated for each year, then averaged'
g <- ggplot(comb, aes(PTcor_rr, PTcor_stepwat)) +
  geom_point() +
  facet_wrap(~RCP) +
  labs(x = 'GridMet/MACA',
       y = 'STEPWAT2',
       caption = cap1,
       subtitle = 'Correlation between monthly temperature and precip,\nfor 200 sites at which STEPWAT simulations are conducted')

g0 <- ggplot(comb, aes(PTcor_stepwat)) +
  geom_histogram() +
  facet_wrap(~RCP)+
  labs(x = 'PT correlation (STEPWAT2)')

g1 <- ggplot(comb, aes(PTcor_rr)) +
  geom_histogram() +
  facet_wrap(~RCP) +
  labs(x = 'PT correlation (MACA/GridMet)')

g3 <- g0/g1
g3

pdf('figures/climate/PTcor_comparison-between-sources.pdf')
g
g3
dev.off()
