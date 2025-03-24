# Martin Holdrege

# Script started May, 15 2024

# Purpose: Create maps of climate variables and change in climate across the west
# based on STEPWAT2 output for 200 sites that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script


# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")

# params ------------------------------------------------------------------

date <- '20250228' # date associated w/ the climate data (this date, in file name,
# would be updated if new climate data source, or sites, etc. were used. )

# some terra operations can be done in parallel and need 
# to know num of cores
num.cores <- parallel::detectCores(logical = FALSE) 
version <- 'v4' # interpolation version

# read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM
clim_files <- list.files(
  file.path("data_processed/interpolated_rasters/climate/", version),
  pattern = paste0('climate.*', date, '_', version, '.tif'),
  full.names = TRUE)
length(clim_files)

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(clim_files) # class SpatRast

# names of layers --------------------------------------------------
# descriptions of layers of raster stack
r_names <- names(rast1)

into <- c("type", "category", "RCP", "years", "GCM", 'date', 'version')

rast_info <- create_rast_info(rast1, into = into, run_regex = '^') %>% 
  # removing trailing GCM, and date
  mutate(id_noGCM = str_replace(id, '_[^_]+_\\d{8}_v\\d$', '')) %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(id) %>% 
  select(-date, -run, -category, -version)

# median by GCM -----------------------------------------------------------
# for each treatment/scenario combination calculate the median
# across GCMs

# split raster by scenario

rast2 <- rast1[[rast_info$id]] # make sure rows are ordered correctly

# list of rasters, each element is a SpatRaster including all GCMs
# for a given treatment/scenario combination
rast_gcm_l <- split(rast2, f = rast_info$id_noGCM)

names(rast_gcm_l) <- unique(rast_info$id_noGCM)
length(rast_gcm_l)

# check that the number of layers in each dataset is correct
stopifnot(map_dbl(rast_gcm_l, nlyr) == rast_info %>%
            group_by(id_noGCM) %>%
            summarise(n = n()) %>%
            pull(n)
)

# note these median rasters now inlcude the median under current
# conditions (i.e. medians of one value). 
med1 <- map(rast_gcm_l, app, fun = "median")

med2 <- rast(med1)

# delta cref ------------------------------------------------------
# change in climate variable relative to ambient climate (c)
# conditions, calculated within a grazing level
# naming: var-rdiff-cref, climate variable raw difference climate refence (i.e. difference
# from historical conditions)

rast_info_med <- create_rast_info(med2, 
                                  into = c("type", "climate", 'RCP', "years"),
                                  run_regex = '^')

info_med_l <- split(rast_info_med, f = rast_info_med$type)

diff_cref1 <- map(info_med_l, function(df) {
  id_current <- df$id[df$RCP == 'Current']
  stopifnot(length(id) ==1)
  ids_future <- df$id[df$RCP != 'Current']
  
  delta <- med2[[ids_future]] - med2[[id_current]]
  new_names <- names(delta) %>% 
    str_replace("_climate", "_rdiff")
  names(delta) <- new_names
  delta
})

diff_cref2 <- rast(diff_cref1)
names(diff_cref2) <- map(diff_cref1, names) %>% 
  unlist()

# save files ---------------------------------------------------------

# * rasters ---------------------------------------------------------------

# median across GCMs, for all future scenarios
writeRaster(med2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0("climate_median_across_GCMs_", version, ".tif")),
            overwrite = TRUE)

# difference in climate (raw) relative to current conditions 
writeRaster(diff_cref2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0("climate_rdiff-cref_median_", version, ".tif")),
            overwrite = TRUE)

