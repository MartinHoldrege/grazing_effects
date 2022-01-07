# Martin Holdrege

# Script started Jan 7, 2022

# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script, which can takes ~ an hour to run on
# my laptop


# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source("src/mapping_functions.R")

# read in data ------------------------------------------------------------

bio_files <- list.files("data_processed/interpolated_rasters/biomass/",
                        pattern = ".tif$",
                        full.names = TRUE)

bio_d_files <- list.files("data_processed/interpolated_rasters/bio_diff/",
                        pattern = ".tif$",
                        full.names = TRUE)

# this is kind of a 'raster stack', where each layers is a tif
rast1 <- terra::rast(c(bio_files, bio_d_files)) # class SpatRast



# names of layers --------------------------------------------------
# descriptions of layers of raster stack
r_names <- names(rast1)

names_df <- tibble(id2 = r_names,
                   id = r_names) %>% 
  separate(col = id2,
           into = c("PFT", "type", "RCP", "years", "graze"),
           sep = "_")

# biomass maps ------------------------------------------------------------


# biomass diff maps -------------------------------------------------------



# sample set of figures ---------------------------------------------------

## Set parameters and layout:
par( mar = c(1,1,2,1), mgp = c(3,0.3,0))
layout.matrix <- matrix(c(1,2,3,4,5,6),nrow = 2, ncol = 3, byrow = T)
layout(layout.matrix, widths = rep(1,3), heights = rep(1,2))


# testing
image_bio(rast1, subset = 1, title = "test")





