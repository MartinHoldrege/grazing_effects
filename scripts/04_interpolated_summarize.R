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
  mutate(layer_num = 1:nrow(.),
         graze = graze2factor(graze),
         years = years2factor(years),
         RCP = rcp2factor(RCP),
         PFT = pft_all_factor(PFT),
         # remove GCM from the string
         id_noGCM = str_replace(id, "_[^_]*$", "")) %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(id)

# id noGCM values, exluding current conditions
id_noGCM_ftr <- rast_info %>% 
  filter(RCP != "Current") %>% 
  pull(id_noGCM) %>% 
  unique()
  

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

med1 <- map(rast_gcm_l[id_noGCM_ftr], 
            app, fun = "median")

med2 <- rast(med1)


# save misc files ---------------------------------------------------------

# not saving files while code testing
saveRDS(rast_info, 
        "data_processed/interpolated_rasters/raster_info.RDS")

# median across GCMs, for all future scenarios
writeRaster(med2, "data_processed/interpolated_rasters/bio_future_median_across_GCMs.tif",
            overwrite = TRUE)
