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

# climate from the weather data base, summarized in the 00_query_weather_db.R script
db_clim1 <- read_csv('data_processed/site_means/dbWeather_200sites.csv')
# summarize across GCMs ---------------------------------------------------

sw_clim_med <- sw_clim1 %>% 
  filter(RCP != 'RCP85',
         years != '2030-2060') %>% 
  group_by(site, years, RCP) %>% 
  # PTcor2 is the mean of yearly correlations
  summarise(PTcor_stepwat = median(PTcor),
            MAT_stepwat = median(MAT),
            MAP_stepwat = median(MAP)) # median across GCMs
# extract clim at points --------------------------------------------------

# convert site locations to sf object
loc2 <- st_as_sf(loc1[, c('site_id', 'X_WGS84', 'Y_WGS84')],
                 crs = 'EPSG:4326',
                 coords = c('X_WGS84', 'Y_WGS84')) %>% 
  st_transform(crs = crs(rr_cur1)) %>% 
  arrange(site_id)


# * PT corr ---------------------------------------------------------------

terra::activeCat(rr_fut1) <- 'CorTempPPT_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'CorTempPPT_mean'

rr_fut2 <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(PTcor_rr = CorTempPPT_mean_acrmoddistr)
rr_cur2 <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(PTcor_rr = CorTempPPT_mean)


# * MAP -------------------------------------------------------------------

terra::activeCat(rr_fut1) <- 'PPT_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'PPT_mean'

rr_fut2a <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(MAP_rr = PPT_mean_acrmoddistr)
rr_cur2a <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(MAP_rr = PPT_mean)


stopifnot(! 'RCP' %in% names(db_clim1)) # db_clim1 is just for current climate
db_clim2 <- db_clim1 %>% 
  select(Site_id, 'MAT_C', 'MAP_mm', matches('CorrTP')) %>% 
  rename( 'MAT_db' = 'MAT_C', 'MAP_db' = 'MAP_mm', PTcor2_db = 'CorrTP2',
          PTcor_db = 'CorrTP') %>% 
  mutate('RCP' = 'Current')
comb <- sw_clim_med %>% 
  left_join(rr_fut2, by = c('site' = 'ID', 'RCP'))%>% 
  left_join(rr_cur2, by = c('site' = 'ID', 'RCP')) %>% 
  mutate(PTcor_rr = ifelse(is.na(PTcor_rr.x), PTcor_rr.y, PTcor_rr.x)) %>% 
  left_join(bind_rows(rr_fut2a, rr_cur2a), by = c('site' = 'ID', 'RCP')) %>% 
  select(-PTcor_rr.y, -PTcor_rr.x) %>% 
  left_join(db_clim2, by =  c('site' = 'Site_id', 'RCP'))


# figures -----------------------------------------------------------------


# *comparing stepwat weather to gridmet -----------------------------------


cap1 = 'The STEPWAT2 (type 1) and the GridMet/MACA P-T correlations (type 2)'
cap2 <- "'type 2' is the correlation between monthly values in a given yr, then avaraged\n type 1 is correlation of monthly means"

g <- ggplot(comb, aes(PTcor_rr, PTcor_stepwat)) +
  geom_point() +
  facet_wrap(~RCP) +
  labs(x = 'GridMet/MACA (type 2)',
       y = 'STEPWAT2 (type 1)',
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

g4 <- ggplot(comb, aes(PTcor_rr, PTcor2_db)) +
  geom_point() +
  labs(x = 'GridMet/MACA (type 2)',
       y = 'Weather database (from 200 sites, type 2)',
       caption = cap2,
       subtitle = 'Correlation (type 2) between monthly temperature and precip,\nfor 200 sites at which STEPWAT simulations are conducted') +
  geom_abline(slope = 1)
g4

g5 <- ggplot(comb, aes(MAP_rr, MAP_stepwat*10)) +
  geom_point() +
  labs(x = 'MAP from GridMet/MACA',
       y = 'MAP from STEPWAT')+
  geom_abline(slope = 1)
# *comparing weather db to stepwat output ----------------------------------
# the climate from the weather db should be closely comparable to what comes
# out of stepwat (which has gone through the wx generator)

base <- function() {
  list(
    geom_abline(slope = 1),
    geom_point(),
    labs(subtitle = 'comparing climate metrics from STEPWAT output\n vs directly from the weather database')
  )
}



pdf('figures/climate/PTcor_comparison-between-sources.pdf')
g
g3
g4
g5
ggplot(comb, aes(MAP_db, MAP_stepwat*10)) +
  base()

ggplot(comb, aes(MAT_db, MAT_stepwat)) +
  base()

ggplot(comb, aes(PTcor_db, PTcor_stepwat))+
  base() +
  labs(caption = cap2,
       x = 'PT cor from database (type 1)',
       y = 'PT cor from stepwat (type 1)')

ggplot(comb, aes(PTcor2_db, PTcor_stepwat))+
  base() +
  labs(caption = cap2,
       x = 'PT cor from database (type 2)',
       y = 'PT cor from stepwat (type 1)')


ggplot(comb, aes(PTcor_db, PTcor2_db))+
  base() + 
  labs(caption = cap2,
     x = 'PT cor from database (type 1)',
     y = 'PT cor from database (type 2)',
     subtitle = 'weather db data (comparing summary methods)')

dev.off()

