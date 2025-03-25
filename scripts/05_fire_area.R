# Purpose: Calculate the expected burned area, and the amount of
# area in different categories of fire return intervals (which we
# consider a proxy for the potential sagebrush 'stand' age in that location)

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


# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source("src/general_functions.R")
source("src/mapping_functions.R")

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
r_area_exp <- (r_prob1/100)*area 



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



# area by fire return interval bin ----------------------------------------

# Use terra::classify to bin probabilities

rcl <- cbind(bp_bins[-1], # from
             bp_bins[-length(bp_bins)],  # to
             1:length(fri_labels)) # becomes

r_binned <- terra::classify(r_prob1, rcl = rcl,
                            right = FALSE)

# calculates the area in each combination of 
# binned fire probability group and ecoregion
fri_area_l <- map(names(r_binned), function(lyr) {
  df <- zonal(area, c(r_binned[[lyr]], r_eco1), 
              fun = 'sum', wide = FALSE)
  df$fri_group <- df[[lyr]]
  df$id <- lyr
  df[[lyr]] <- NULL
  df
})

fri_area1 <- fri_area_l %>% 
  bind_rows() %>% 
  as_tibble()
 
# adding 0s ---------------------------------------------------------------

# adding missing (zero area) back in

cols2expand <- map(fri_area1[c("fri_group", "ecoregion", "id")], unique) 
  
expanded <- expand_grid(!!!cols2expand)
  
fri_area2 <- fri_area1 %>% 
  right_join(expanded,
             by = join_by(ecoregion, fri_group, id)) %>% 
  mutate(area = ifelse(is.na(area), 0, area),
         fri_label = fri_labels[fri_group]) %>% 
  left_join(info_prob1, by = 'id') %>% 
  select(-id)
  
# save output -------------------------------------------------------------

write_csv(ba_eco2, 
          paste0("data_processed/area/expected-burn-area_", v, "_", run, 
                 ".csv"))

write_csv(fri_area2, 
          paste0("data_processed/area/area-by-fri-class_", v, "_", run, 
                 ".csv"))
