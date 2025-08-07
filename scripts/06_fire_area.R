# Purpose: Calculate the expected burned area, and the amount of
# area in different categories of expected stand age (based on fire probability)
# additionally expected burned area by region and current SCD group (based
# on pixelwise and gcm-wise summaries)

# Author: Martin Holdrege

# Date started: March 24, 2025

# params ------------------------------------------------------------------

source('src/params.R')

test_run <- opt$test_run
vr_name <- opt$vr_name
run <- opt$run
v_interp <- opt$v_interp
v <- paste0('v2', opt$vr_name) # version for output files
runv <- paste0(opt$run, opt$v_interp)
ref_graze <- 'Moderate' # the reference grazing level for SEI calculations
# (i.e. current climate at ref graze means we are not adjusting
# the sei values with stepwat results)
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source("src/general_functions.R")
source("src/mapping_functions.R")
source("src/probability_functions.R")
source("src/SEI_functions.R")

# read in data ------------------------------------------------------------

# fire probability for each GCM

# data up-scaled for each GCM
fire_files <- list.files(
  file.path("data_processed/interpolated_rasters/fire/", v_interp),
  pattern = run, full.names = TRUE)

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
r_prob_gcm1 <- terra::rast(fire_files) # class SpatRast

into_gcm <- c("type", "RCP", "years",  "graze", "GCM")
info_gcm1 <- create_rast_info(r_prob_gcm1, into = into_gcm)

# fire probability summarized across GCMs
r_prob_smry1 <- rast(paste0("data_processed/interpolated_rasters/", 
                            run, "_fire-prob_future_summary_across_GCMs.tif"))
info_smry1 <- create_rast_info(r_prob_smry1, 
                               into = c("type", "RCP", "years",  "graze", 
                                        "summary"))

# computing wise it would be faster to have total_region = FALSE
# and then sum the summary stats for the individual
# ecoregions, instead of then seperately computing for the total_region
# (doing this for now, because it requires fewer lines of code)
eco1 <- load_wafwa_ecoregions(total_region = TRUE, v = vr)

r_eco1 <- load_wafwa_ecoregions_raster(v = vr) # raster of ecoregions

# * SEI ------------------------------------------------------------

r_sei1 <- rast(file.path("data_processed/interpolated_rasters/", v_interp,
                    paste0(runv, "_q-sei_scd-adj_summary.tif")))

band <- names(r_sei1) %>% 
  stringr::str_subset(paste0('_SEI_SEI_Current_Current_', ref_graze, '_median'))
stopifnot(length(band) == 1)
r_sei2 <- r_sei1[[band]]

# * if test run -----------------------------------------------------------

if(test_run) {
  r_prob_gcm1 <- downsample(r_prob_gcm1) 
  r_prob_smry1 <- downsample(r_prob_smry1)
  r_eco1 <- downsample(r_eco1)
  r_sei2 <- downsample(r_sei2)
}

# vectors etc -------------------------------------------------------------

# for area by age class group
# these intervals of years
# to answer the question of (e.g) what is the expected area of the region
# that will have burned between 0 and 10 years ago
age_groups <- create_age_groups()


# the probability of falling into any one of the categories should be 1
# (b/ the categories include all possible years)
test <- map_dbl(age_groups, function(x) {
  prob_geometric_rule(p = 0, a = x['a'], b = x['b'])
})
stopifnot(all.equal(sum(test), 1))

# prepare rasters ---------------------------------------------------------

stopifnot(isTRUE(same.crs(r_prob_gcm1, vect(eco1))))

# expected burned area ----------------------------------------------------

area <- cellSize(r_prob_gcm1[[1]], unit = 'ha', transform = FALSE)
varnames(area) <- 'cellSize'


# * gcm-wise --------------------------------------------------------------

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
  # using FALSE, otherwise summing indidividual ecoregions equals more area than total area
  touches = FALSE 
)

ba_eco1$ecoregion <- eco1$ecoregion

ba_eco2 <- ba_eco1 %>% 
  select(-ID) %>% 
  pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') %>% 
  left_join(info_gcm1, by = 'id') %>% 
  select(-id)


# * pixewlise -------------------------------------------------------------

r_prob_smry2 <- r_prob_smry1/100 # convert % to proportion

r_area_exp_smry <- r_prob_smry2*area 

# expected burned area per ecoregion
ba_eco_smry1 <- terra::extract(
  r_area_exp_smry ,
  vect(eco1),           # convert sf to SpatVector
  fun = sum,
  na.rm = TRUE,
  touches = FALSE 
)

ba_eco_smry1$ecoregion <- eco1$ecoregion

ba_eco_smry2 <- ba_eco_smry1 %>% 
  select(-ID) %>% 
  pivot_longer(cols = -ecoregion,
               names_to = 'id',
               values_to = 'area') %>% 
  left_join(info_smry1, by = 'id') %>% 
  select(-id)

ba_eco_smry3 <- ba_eco_smry2 %>% 
  pivot_wider(values_from = 'area',
              names_from = 'summary',
              names_prefix = 'area_') %>% 
  select(-type)

# * by current SEI class --------------------------------------------------

c3 <- sei2c3(r_sei2)
names(c3) <- 'c3'

c3eco <- c3*10 + as.numeric(r_eco1)
names(c3eco) <- 'c3eco'
c3eco_area1 <- zonal(area, c3eco, fun = 'sum')

# **gcm-wise --------------------------------------------------------------

c3eco_ba_gcm1 <- zonal(r_area_exp, c3eco, fun = 'sum')

c3eco_ba_gcm2 <- c3eco_ba_gcm1  %>% 
  pivot_longer(-c3eco,
               values_to = 'expected_ba',
               names_to = 'id') %>% 
  mutate(c3 = c3eco_to_c3(c3eco),
         region = c3eco_to_eco(c3eco, levels(r_eco1)[[1]]$ecoregion)) %>% 
  left_join(c3eco_area1, by = "c3eco") %>% 
  left_join(info_gcm1, by = 'id') %>% 
  select(-id, c3eco)

# ** pixelwise ------------------------------------------------------------

c3eco_ba_smry1 <- zonal(r_area_exp_smry, c3eco, fun = 'sum')

c3eco_ba_smry2 <- c3eco_ba_smry1  %>% 
  pivot_longer(-c3eco,
               values_to = 'expected_ba',
               names_to = 'id') %>% 
  mutate(c3 = c3eco_to_c3(c3eco),
         region = c3eco_to_eco(c3eco, levels(r_eco1)[[1]]$ecoregion)) %>% 
  left_join(c3eco_area1, by = "c3eco") %>% 
  left_join(info_smry1, by = 'id') %>% 
  select(-id, -c3eco)

# area by by stand age ----------------------------------------------------

# *gcm-wise ---------------------------------------------------------------

# probability of each pixel of falling into a given age group
prob_by_age_group_l <- map(age_groups, function(x) {
  prob_geometric_rule(r_prob_gcm2, a = x['a'], b = x['b'])
})

# excpected area (for each pixel) by age group
area_by_age_group_l <- map(prob_by_age_group_l, \(x) x*area)

area_age_group1 <- map(area_by_age_group_l, 
                       \(x) terra::extract(x, vect(eco1), fun = sum, 
                                           na.rm = TRUE, touches = FALSE))

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


# * pixelwise median fire prob --------------------------------
# using the abbreviation 'pw' to denote pixelwise
# these summaries are same as above but for median fire probability
# this way the median fire probability summed across groups
# sums to area of region

info_smry2 <- info_smry1 %>% 
  filter(summary == 'median')



# probability of each pixel of falling into a given age group
prob_by_age_group_l_pw <- map(age_groups, function(x) {
  prob_geometric_rule(r_prob_smry2[[info_smry2$id]], a = x['a'], b = x['b'])
})

# excpecte area (for each pixel) by age group
area_by_age_group_l_pw <- map(prob_by_age_group_l_pw, \(x) x*area)

area_age_group1_pw <- map(area_by_age_group_l_pw, 
                       \(x) terra::extract(x, vect(eco1), fun = sum, 
                                           na.rm = TRUE, touches = FALSE))

area_age_group2_pw <- map(area_age_group1_pw, function(x) {
  x$ecoregion <- eco1$ecoregion[x$ID]
  x$ID <- NULL
  
  x %>% 
    pivot_longer(cols = -ecoregion,
                 names_to = 'id',
                 values_to = 'area') 
  
}) %>% 
  bind_rows(.id = 'age_group') %>% 
  left_join(info_smry1, by = 'id')

# renaming (and widening if there were multiple summary stats)
# so format consistent with area_age_group3
area_age_group3_pw <- area_age_group2_pw %>% 
  pivot_wider(values_from = 'area',
              names_from = 'summary',
              names_prefix = 'area_')

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
  touches = FALSE
)
area_eco$ecoregion <- eco1$ecoregion
area_eco$ID <- NULL


# save output -------------------------------------------------------------

if(!test_run) {
  write_csv(ba_eco_smry3, 
            paste0("data_processed/area/expected-burn-area_smry_", v, "_", run, 
                   ".csv"))
  
  write_csv(c3eco_ba_smry2, 
            paste0("data_processed/area/expected-burn-area_by-c3-smry_", v, "_", run, 
                   ".csv"))
  
  write_csv(c3eco_ba_gcm2, 
            paste0("data_processed/area/expected-burn-area_by-c3-GCM_", v, "_", run, 
                   ".csv"))
  
  write_csv(ba_eco2, 
            paste0("data_processed/area/expected-burn-area_by-GCM_", v, "_", run, 
                   ".csv"))
  
  write_csv(area_eco, 
            paste0("data_processed/area/ecoregion-area_", v,".csv"))
  
  write_csv(area_age_group3, 
            paste0("data_processed/area/area-by-age-group_", v, "_", run, 
                   ".csv"))
  
  write_csv(area_age_group3_pw, 
            paste0("data_processed/area/area-by-age-group_pw_", v, "_", run, 
                   ".csv"))
}
