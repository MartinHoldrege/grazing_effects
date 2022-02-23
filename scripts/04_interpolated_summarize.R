# Martin Holdrege

# Script started Jan 7, 2022

# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script, which can takes ~ 4 hours to run on
# my laptop (when excluding RCP4.5)


# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")

# read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM, and for c4on and c4off
bio_files <- list.files("data_processed/interpolated_rasters/biomass/",
                        pattern = "^c4.*tif$",
                        full.names = TRUE)
length(bio_files)

# for testing use a subset for computational efficiency
# bio_files <- bio_files[sample(1:length(bio_files), 100)]

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(bio_files) # class SpatRast


# params ------------------------------------------------------------------

# some terra operations can be done in parallel and need 
# to know num of cores
num.cores <- parallel::detectCores(logical = FALSE) 

# names of layers --------------------------------------------------
# descriptions of layers of raster stack
r_names <- names(rast1)

rast_info <- tibble(id2 = r_names,
                   id = r_names) %>% 
  separate(col = id2,
           into = c("c4", "PFT", "type", "RCP", "years", "graze", "GCM"),
           sep = "_") %>% 
  df_factor() %>% 
  mutate(layer_num = 1:nrow(.),
         # remove GCM from the string
         id_noGCM = str_replace(id, "_[^_]*$", "")) %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(id)


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

# SpatRasterDataset (i.e. each dataset includes all GCMs for treatment/scenario)
# sds_gcm <- sds(rast_gcm_l) 
# class(sds_gcm)
# 
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

#  c3Pgrass/Pgrass -------------------------------------------------------

grass_info <- rast_info %>% 
  filter(PFT %in% c("C3Pgrass", "Pgrass"),
         # just looking at end of century and current for now
         years %in% c("2070-2100", "Current")) %>% 
  arrange(PFT, c4, graze, RCP, years) %>% # important so comparing correct lyrs
  # creates two identical df's except the PFT and id columns are different
  split(f = .$PFT, drop = TRUE)


grass_id <- map(grass_info, function(x) x$id) # id values
grass_id_noGCM <- map(grass_info, function(x) x$id_noGCM) # id values

# seperate rasters for Pgrass and C3Pgrass
grass_r1 <- map(grass_id, function(x) rast1[[x]])

# grass ratio 
gratio_r1 <- grass_r1$C3Pgrass/grass_r1$Pgrass

# check that number of layers has been preserved
stopifnot(nlyr(gratio_r1) == length(grass_id[[1]]))

# ** median ratio ----------------------------------------------------------

# info on the median lyrs
grass_info_med <- grass_info$C3Pgrass %>% 
  arrange(graze, RCP, years, c4) %>% 
  select(-id, -GCM, -layer_num) %>% 
  filter(!duplicated(.))

# split into list
gratio_r_l <- split(gratio_r1, grass_id_noGCM$C3Pgrass)

# naming so layers will be identifiable, C3Pgrass is in name
# but this actuall C3Pgrass/Pgrass ratio
names(gratio_r_l) <- unique(grass_id_noGCM$C3Pgrass)

# take median across GCMs
gratio_med1 <- map(gratio_r_l, app, fun = "median")

# this is operation takes a few minutes
gratio_med2 <- rast(gratio_med1) # put back into one SpatRaster
nlyr(gratio_med2)
names(gratio_med2)

# ** Pgrass ----------------------------------------------------------------

Pgrass_id_noGCM <- names(med2) %>% 
  str_subset("_Pgrass_") %>% #on Pgrass
  str_subset("_2030-2060_", negate = TRUE) # not mid century

# information on Pgrass, medians across GCMs (same as info for C3Pgrass)
Pgrass_info_med<- grass_info_med %>% 
  mutate(id_noGCM = str_replace(id_noGCM, "_C3Pgrass","_Pgrass")) %>% 
  select(-PFT)

Pgrass_target_lyrs <- Pgrass_info_med %>% 
  filter(RCP != "Current") %>% 
  pull(id_noGCM)

Pgrass_ref_lyrs <- create_ref_id(Pgrass_target_lyrs)

# scaled % change relative to current grazing of the same intensity
rast_d_Pgrass <- rast_diff(rast = med2,
                           target_layer = Pgrass_target_lyrs,
                           ref_layer = create_ref_id(Pgrass_target_lyrs))

names(rast_d_Pgrass)


# save files ---------------------------------------------------------

# * info files ------------------------------------------------------------

saveRDS(rast_info, 
        "data_processed/interpolated_rasters/raster_info.RDS")

saveRDS(grass_info_med,
        "data_processed/interpolated_rasters/grass_info_med.RDS")

# info on median pgrass
saveRDS(Pgrass_info_med,
        "data_processed/interpolated_rasters/Pgrass_info_med.RDS")

# * rasters ---------------------------------------------------------------

# median across GCMs, for all future scenarios
writeRaster(med2, "data_processed/interpolated_rasters/bio_future_median_across_GCMs.tif",
            overwrite = TRUE)


# c3/total grass ratio
writeRaster(gratio_med2, "data_processed/interpolated_rasters/C3Pgrass-Pgrass-ratio_by-scenario_median.tif",
            overwrite = TRUE)

# scaled percent change from current to future (wgraze), for Pgrass

writeRaster(rast_d_Pgrass, "data_processed/interpolated_rasters/Pgrass_bio-diff-wgraze_median.tif",
            overwrite = TRUE)
