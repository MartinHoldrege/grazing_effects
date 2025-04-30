# Purpose: Calculate the number of pixels each site is interpolated to
# in each ecoregion, and across the entire study area. The purpose is
# to come up with weights, these weights are then used for summary statistics
# plotting (boxplots etc)

# Author: Martin Holdrege

# Started: April 11, 2025

# params ------------------------------------------------------------------

v <- 'v4' # interpolation version

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source("src/mapping_functions.R")

# load data ---------------------------------------------------------------

r_eco <- load_wafwa_ecoregions_raster()

locs1 <- rast(paste0("data_processed/interpolation_data/interp_locations_200sites_",
                     v, ".tif"))

match1 <- rast(paste0("data_processed/interpolation_quality/",
                      "matching_quality", v, '.tif'))

# prepare rasters ---------------------------------------------------------

locs2 <- locs1
locs2[match1 > 1.5] <-NA

# calculate weights --------------------------------------------------------
# weights are the number of pixels a given site is interpolated to

area_wide0 <- as.data.frame(locs2)


area <- cellSize(locs2, unit = 'ha', transform = FALSE)
pixel_size <- unique(as.vector(values(area)))
stopifnot(length(pixel_size) == 1)

df1 <- tibble(ecoregion = as.vector(values(r_eco)) + 1, # for some reason converted to 0-2,
              site = as.vector(values(locs2)),
              area = as.vector(area))

# there are some pixels that fall just outside the boundary 
# of the ecoregions, for now these pixels are being discarded
# for this analyses. My guess is that they've arisen b/ of aggregation
# of pixels from 30 m to 1km when creating the study area map
# or some other projection change causing slight shifts at the margins
check <- sum(is.na(df1$ecoregion) & !is.na(df1$site))
stopifnot(check < 600)

lvls <- levels(r_eco)[[1]]
stopifnot(lvls$ID == sort(lvls$ID)) # below code only words if lvls is ordered
df2 <- df1 %>% 
  drop_na() %>% 
  mutate(region = lvls$ecoregion[ecoregion])

weights1 <- df2 %>% 
  group_by(region, site) %>% 
  summarize(weight = sum(area)/pixel_size,
            .groups = 'drop')

# run an area comparison
if (FALSE) {
  # here area is tied to entire pixels
  # note-still unclear how these two approaches differ under the hood
  weights1 %>% 
    group_by(region) %>% 
    summarize(area = sum(weight)*pixel_size)
  
  # here area calculated differently
  check <- read_csv(paste0("data_processed/area/ecoregion-area_v1.csv"))
  check
}

weights2 <- weights1 %>% 
  group_by(site) %>% 
  summarize(weight = sum(weight)) %>% 
  mutate(region = 'Entire study area') %>% 
  bind_rows(weights1)

stopifnot(1:200 %in% weights2$site) # make sure all sites show up


# save output -------------------------------------------------------------

write_csv(weights2, 
          paste0('data_processed/interpolation_data/interpolation_weights_', 
                 v, '.csv'))

