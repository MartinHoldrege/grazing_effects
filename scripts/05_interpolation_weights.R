# Purpose: Calculate the number of pixels each site is interpolated to
# in each ecoregion, and across the entire study area. The purpose is
# to come up with weights, these weights are then used for summary statistics
# plotting (boxplots etc)

# Author: Martin Holdrege

# Started: April 11, 2025

# params ------------------------------------------------------------------

source('src/params.R')
v <- opt$v_interp # interpolation version
vr <- opt$vr
vr_name <- opt$vr_name

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source("src/mapping_functions.R")
source('src/SEI_functions.R')
# load data ---------------------------------------------------------------

r_eco <- load_wafwa_ecoregions_raster(v = vr)

locs1 <- rast(paste0("data_processed/interpolation_data/interp_locations_200sites_",
                     v, ".tif"))

c3eco <- load_c3eco(v = vr)

# prepare rasters ---------------------------------------------------------

locs2 <- terra::mask(locs1, c3eco)

# calculate weights --------------------------------------------------------

c3eco_num <- as.numeric(c3eco)
names(c3eco_num) <- 'c3eco'
df_comb <- as.data.frame(c(c3eco_num, locs2))
stopifnot(all(!is.na(df_comb)))

wt_c3eco1 <- df_comb %>% 
  rename(site = site_id) %>% 
  left_join(cats(c3eco)[[1]], by = c('c3eco' = 'ID')) %>% 
  group_by(c3eco, site, region, c3) %>% 
  summarize(weight = n(),
            .groups = 'drop') 

# calculating for 'entire' region
wt_c3eco2 <- wt_c3eco1 %>% 
  group_by(site, c3) %>% 
  summarize(weight = sum(weight),
            .groups = 'drop') %>% 
  mutate(region = region_factor(return_levels = TRUE, v = vr)[1]) %>% 
  bind_rows(wt_c3eco1) %>% 
  df_factor() %>% 
  arrange(region, c3, site) %>% 
  select(-c3eco)

# calculating by region (across sei classes)
wt_eco <- wt_c3eco2 %>% 
  group_by(site, region) %>% 
  summarize(weight = sum(weight),
            .groups = 'drop') %>% 
  arrange(region, site)


stopifnot(1:200 %in% wt_eco$site) # make sure all sites show up


# save output -------------------------------------------------------------

write_csv(wt_eco, 
          paste0('data_processed/interpolation_data/interpolation_weights_', 
                 v, vr_name, '.csv'))

write_csv(wt_c3eco2, 
          paste0('data_processed/interpolation_data/interpolation_weights_c3_', 
                 v, '_', vr, '.csv'))
