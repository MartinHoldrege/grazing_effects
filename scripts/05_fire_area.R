# Purpose: Calculate the expected burned area, and the amount of
# area in different categories of expected stand age (based on fire probability)

# Author: Martin Holdrege

# Date started: March 24, 2025

# params ------------------------------------------------------------------

run <- 'fire1_eind1_c4grass1_co20_2503'
test_run <- FALSE

v <- 'v1' # version for output files

v_interp <- 'v4' # interpolation version (input files)

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source("src/general_functions.R")
source("src/mapping_functions.R")
source("src/probability_functions.R")

# read in data ------------------------------------------------------------

# fire probability for each GCM

# data up-scaled for each GCM
fire_files <- list.files(
  file.path("data_processed/interpolated_rasters/fire/", v_interp),
  pattern = run, full.names = TRUE)

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
r_prob_gcm1 <- terra::rast(fire_files) # class SpatRast

if(test_run) {
  r_prob_gcm1 <- r_prob_gcm1[[1:2]]
}

into_gcm <- c("type", "RCP", "years",  "graze", "GCM")
info_gcm1 <- create_rast_info(r_prob_gcm1, into = into_gcm)

eco1 <- load_wafwa_ecoregions()

# vectors etc -------------------------------------------------------------

# for area by age class group
# these intervals of years
# to answer the question of (e.g) what is the expected area of the region
# that will have burned between 0 and 10 years ago
age_groups <- create_age_groups()


# the probability of falling into any one of the categories should be 1
# (b/ the categories include all possible years)
test <- map_dbl(age_groups, function(x) {
  prob_geometric_rule(p = 0.01, a = x['a'], b = x['b'])
})
stopifnot(all.equal(sum(test), 1))

# prepare rasters ---------------------------------------------------------

stopifnot(isTRUE(same.crs(r_prob_gcm1, vect(eco1))))

r_eco1 <- rasterize(
  vect(eco1),
  r_prob_gcm1[[1]],
  field = 'ecoregion',
  touches = TRUE
)

# expected burned area ----------------------------------------------------

area <- cellSize(r_prob_gcm1[[1]], unit = 'ha', transform = FALSE)
varnames(area) <- 'cellSize'
# amount of area expected to burn in each
# first convert probability from % to proportion
r_prob_gcm2 <- r_prob_gcm1/100 
r_area_exp <- r_prob_gcm2*area 

# expected burned area per ecoregion
ba_eco1 <- terra::extract(
  r_area_exp,
  vect(eco1),           # convert sf to SpatVector
  fun = sum,
  na.rm = TRUE,
  touches = TRUE
)

ba_eco1$ecoregion <- eco1$ecoregion

ba_eco2 <- ba_eco1 %>% 
  select(-ID) %>% 
  pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') %>% 
  left_join(info_gcm1, by = 'id') %>% 
  select(-id)

ba_eco3 <- ba_eco2 %>% 
  group_by(ecoregion, run, RCP, years, graze) %>% 
  summarize(area_median = median(area),
            area_low = calc_low(area),
            area_high = calc_high(area),
            .groups = "drop")

# area by by stand age ----------------------------------------------------

# probability of each pixel of falling into a given age group
prob_by_age_group_l <- map(age_groups, function(x) {
  prob_geometric_rule(r_prob_gcm2, a = x['a'], b = x['b'])
})

# excpecte area (for each pixel) by age group
area_by_age_group_l <- map(prob_by_age_group_l, \(x) x*area)

area_age_group1 <- map(area_by_age_group_l, 
                       \(x) terra::extract(x, vect(eco1), fun = sum, 
                                           na.rm = TRUE, touches = TRUE))

area_age_group2 <- map(area_age_group1, function(x) {
  x$ecoregion <- eco1$ecoregion[x$ID]
  x$ID <- NULL
  
  x %>% 
    pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') 
    
}) %>% 
  bind_rows(.id = 'age_group') %>% 
  left_join(info_gcm1, by = 'id')


# summarize across GCMs
area_age_group3 <- area_age_group2 %>% 
  group_by(age_group, ecoregion, run, RCP, years, graze) %>% 
  summarize(area_median = median(area),
            area_low = calc_low(area),
            area_high = calc_high(area),
            .groups = "drop")

# area by ecoregion -------------------------------------------------------
# for % of area calculations
mask <- r_prob_gcm1[[1]]
mask[!is.na(mask)] <- 1
area_masked <- mask(area, mask)

area_eco <- terra::extract(
  area_masked,
  vect(eco1),         
  fun = sum,
  na.rm = TRUE, 
  touches = TRUE
)
area_eco$ecoregion <- eco1$ecoregion
area_eco$ID <- NULL


# save output -------------------------------------------------------------

if(!test_run) {
  write_csv(ba_eco3, 
            paste0("data_processed/area/expected-burn-area_", v, "_", run, 
                   ".csv"))
  
  write_csv(ba_eco2, 
            paste0("data_processed/area/expected-burn-area_by-GCM_", v, "_", run, 
                   ".csv"))
  
  write_csv(area_eco, 
            paste0("data_processed/area/ecoregion-area_", v,".csv"))
  
  write_csv(area_age_group3, 
            paste0("data_processed/area/area-by-age-group_", v, "_", run, 
                   ".csv"))
}
