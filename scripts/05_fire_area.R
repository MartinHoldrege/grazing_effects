# Purpose: Calculate the expected burned area, and the amount of
# area in different categories of expected stand age (based on fire probability)

# Author: Martin Holdrege

# Date started: March 24, 2025


# params ------------------------------------------------------------------

run <- 'fire1_eind1_c4grass1_co20_2502'
test_run <- FALSE

v <- 'v1' # version for output files

# burn probabilities, and associated fire return intervals
# for area calculations
bp_bins <- c(Inf, 3, 2, 1, 0.5, 0) # neg probability in last bin so includes 0
fri_labels <- c("< 33", '33 - 50', '50 - 100', '100 - 200', '200 +')

# for area by age class group
# these intervals of years
# to answer the question of (e.g) what is the expected area of the region
# that will have burned between 0 and 10 years ago
age_groups <- list(
  c(a = 0, b = 10),
  c(a = 11, b = 25),
  c(a = 26, b = 75),
  c(a = 76, b = 150),
  c(a = 151, b = Inf)
)

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source("src/general_functions.R")
source("src/mapping_functions.R")
source("src/probability_functions.R")

# read in data ------------------------------------------------------------

path_r <- "data_processed/interpolated_rasters"
path <- file.path(path_r, 
                  paste0(run, '_fire-prob_future_median_across_GCMs.tif'))

r_prob1 <- rast(path)
if(test_run) {
  r_prob1 <- r_prob1[[1:2]]
}
into <- c("type", "RCP", "years", 
          "graze")

info_prob1 <- create_rast_info(r_prob1, into = into)

eco1 <- load_wafwa_ecoregions()



# vectors etc -------------------------------------------------------------

names(age_groups) <- map_chr(age_groups, function(x) {
  out <- if(is.finite(x['b'])) {
    paste0(x['a'], '-',  x['b'])
  } else {
    paste0(x['a'], ' +')
  }
  out
})

# the probability of falling into any one of the categories should be 1
# (b/ the categories include all possible years)
test <- map_dbl(age_groups, function(x) {
  prob_geometric_rule(p = 0.01, a = x['a'], b = x['b'])
})
stopifnot(all.equal(sum(test), 1))

# prepare rasters ---------------------------------------------------------

stopifnot(isTRUE(same.crs(r_prob1, vect(eco1))))

r_eco1 <- rasterize(
  vect(eco1),
  r_prob1[[1]],
  field = 'ecoregion'
)

# expected burned area ----------------------------------------------------

area <- cellSize(r_prob1[[1]], unit = 'ha', transform = FALSE)
varnames(area) <- 'cellSize'
# amount of area expected to burn in each
# first convert probability from % to proportion
r_prob2 <- r_prob1/100 
r_area_exp <- r_prob2*area 

# expected burned area per ecoregion
ba_eco1 <- terra::extract(
  r_area_exp,
  vect(eco1),           # convert sf to SpatVector
  fun = sum,
  na.rm = TRUE
)

ba_eco1$ecoregion <- eco1$ecoregion

ba_eco2 <- ba_eco1 %>% 
  select(-ID) %>% 
  pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') %>% 
  left_join(info_prob1, by = 'id')

# area by by stand age ----------------------------------------------------

# probability of each pixel of falling into a given age group
prob_by_age_group_l <- map(age_groups, function(x) {
  prob_geometric_rule(r_prob2, a = x['a'], b = x['b'])
})

# excpecte area (for each pixel) by age group
area_by_age_group_l <- map(prob_by_age_group_l, \(x) x*area)

area_age_group1 <- map(area_by_age_group_l, 
                       \(x) terra::extract(x, vect(eco1), fun = sum, 
                                           na.rm = TRUE))


area_age_group2 <- map(area_age_group1, function(x) {
  x$ecoregion <- eco1$ecoregion[x$ID]
  x$ID <- NULL
  
  x %>% 
    pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') 
    
}) %>% 
  bind_rows(.id = 'age_group') %>% 
  left_join(info_prob1, by = 'id')


# area by ecoregion -------------------------------------------------------
# for % of area calculations
mask <- r_prob1[[1]]
mask[!is.na(tmp)] <- 1
area_masked <- mask(area, mask)

area_eco <- terra::extract(
  area_masked,
  vect(eco1),         
  fun = sum,
  na.rm = TRUE
)
area_eco$ecoregion <- eco1$ecoregion
area_eco$ID <- NULL

# save output -------------------------------------------------------------

if(!test_run) {
  write_csv(ba_eco2, 
            paste0("data_processed/area/expected-burn-area_", v, "_", run, 
                   ".csv"))
  
  write_csv(area_eco, 
            paste0("data_processed/area/ecoregion-area_", v,".csv"))
  
  write_csv(area_age_group2, 
            paste0("data_processed/area/area-by-age-group_", v, "_", run, 
                   ".csv"))
}
