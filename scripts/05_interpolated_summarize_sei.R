# Martin Holdrege

# Script started April 18, 2025

# Purpose: Calculate Q and SEI for each GCM, then summarize across
# GCMs. Note this needs to be calculated for each GCM, b/ median
# SEI need arise from the median biomass of each component of SEI

# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")
source('src/SEI_functions.R')
# params ------------------------------------------------------------------

v <- 'v4' # interpolation version
# multiple runs may need to listed if different runs done for different
# grazing levels

run <- 'fire1_eind1_c4grass1_co20_2503'

pfts <- c('Sagebrush', 'Pherb', 'Aherb')

test_run <- FALSE

# read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM

regex <- paste0(run, '.*(', paste0(pfts, collapse = '|'), ")")
bio_files <- list.files(file.path("data_processed/interpolated_rasters/biomass", v),
                        pattern = regex,
                        full.names = TRUE)
length(bio_files)

# for testing use a subset for computational efficiency
# bio_files <- bio_files[sample(1:length(bio_files), 100)]

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(bio_files) # class SpatRast

r_eco1 <- load_wafwa_ecoregions_raster()

# params ------------------------------------------------------------------

# some terra operations can be done in parallel and need 
# to know num of cores
num.cores <- parallel::detectCores(logical = FALSE) 

# names of layers --------------------------------------------------

rast_info <- create_rast_info(rast1, id_noGCM = TRUE) %>% 
  filter_clim_extremes() %>% 
  mutate(layer_num = 1:nrow(.)) %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(graze, RCP, years, GCM, id) %>% 
  mutate(id_noPFT = str_replace(id, paste0("_(", paste0(pfts, collapse = '|'), ")"), 
                                ""))
# Calculate Q scores and SEI ----------------------------------------------

rast2 <- rast1[[rast_info$id]] # only keep clim scenarios of interest
rast2 <- mask(rast2, r_eco1)

info_l <- rast_info %>% 
  group_by(id_noPFT) %>% 
  group_split()

if(test_run) {
  info_l <- info_l[1:2]
}

test <- map_dbl(info_l, nrow)
stopifnot(test == 3)

r_sei_l1 <- map(info_l, function(df) {
  bio2qsei_raster(r = rast2[[df$id]],
                  eco_raster = r_eco1, 
                  type_string = 'biomass')
})
  
r_sei2 <- rast(r_sei_l1)

rast_info2 <- create_rast_info(r_sei2 , 
                               # here 'group' is the PFT or NA for SEI
                               into = c("group", "type", "RCP", "years", 
                                          "graze", "GCM"),
                               id_noGCM = TRUE)

# median by GCM -----------------------------------------------------------
# for each treatment/scenario combination calculate the median
# across GCMs

# list of rasters, each element is a SpatRaster including all GCMs
# for a given treatment/scenario combination
rast_gcm_l <- split(r_sei2, f = rast_info2$id_noGCM)

names(rast_gcm_l) <- unique(rast_info2$id_noGCM)

# note these median rasters now inlcude the median under current
# conditions (i.e. medians of one value). 
med2 <- rast(map(rast_gcm_l, app, fun = "median"))

# delta biomass cref ------------------------------------------------------
# change in q/sei (absolulte, not relative) relative to ambient climate (c)
# conditions, calculated within a grazing level
# naming: q-sei-rdiff-cref, biomass raw difference climate refence (i.e. difference
# from historical conditions)

rast_info_med <- create_rast_info(med2, 
                                  into = c("group", "type", "RCP", "years", 
                                           "graze")) %>% 
  mutate(id3 = paste(run2, group, type, sep = "_"))

info_med_l <- split(rast_info_med, f = rast_info_med$id3)

diff_cref1 <- map(info_med_l, function(df) {
  id_current <- df$id[df$RCP == 'Current']
  stopifnot(length(id) ==1)
  ids_future <- df$id[df$RCP != 'Current']
  
  delta <- med2[[ids_future]] - med2[[id_current]]
  new_names <- names(delta) %>% 
    str_replace("_Q", "_Q-rdiff-cref") %>% 
    str_replace('_SEI', "_SEI-rdiff-cref")
  names(delta) <- new_names
  delta
})

diff_cref2 <- rast(diff_cref1)
names(diff_cref2) <- map(diff_cref1, names) %>% 
  unlist()

# write rasters ---------------------------------------------------------------

# median across GCMs, for all future scenarios

names(med2) <- paste0(names(med2), '_median')

writeRaster(med2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0(run, "_q-sei_future_summary_across_GCMs.tif")),
            overwrite = TRUE)

# difference in sei/q (raw) relative to current conditions (within a grazing level)
names(diff_cref2) <- paste0(names(diff_cref2), '_median')

writeRaster(diff_cref2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0(run, "_q-sei-rdiff-cref_summary.tif")),
            overwrite = TRUE)

