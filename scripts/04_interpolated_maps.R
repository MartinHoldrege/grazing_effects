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
source("src/general_functions.R")

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
           sep = "_") %>% 
  mutate(layer_num = 1:nrow(.),
         graze = graze2factor(graze),
         years = years2factor(years),
         RCP = rcp2factor(RCP),
         PFT = pft5_factor(PFT))

# names of layers for first set of of maps
# just want to compare current biomass to delta biomass in the future under
# RCP 8.5
names_plot_rcp8.5 <- names_df %>% 
  filter((RCP == "Current" & type == "biomass") | 
           (RCP == "RCP8.5" & type == "bio-diff")) %>% 
  arrange(PFT, graze, years)

names_plot_rcp4.5 <- names_df %>% 
  filter((RCP == "Current" & type == "biomass") | 
           (RCP == "RCP4.5" & type == "bio-diff")) %>% 
  arrange(PFT, graze, years)

# so can plot rcp 4.5 and 8.5 in separate (current conditions are plotted
# in both cases)
names_plot1 <- bind_rows(names_plot_rcp4.5, names_plot_rcp8.5)


# sample set of figures ---------------------------------------------------

## Set parameters and layout:
par( mar = c(1,1,2,1), mgp = c(3,0.3,0))
layout.matrix <- matrix(c(1,2,3,4,5,6),nrow = 2, ncol = 3, byrow = T)
layout(layout.matrix, widths = rep(1,3), heights = rep(1,2))


# testing
# either pass the name of the layer or the layer number
image_bio(rast1, subset = names(rast1)[1], title = "test")

image_bio_diff(rast1, subset = 150, title = "test2")


# diff from current, same grazing level -----------------------------------
# each page shows current biomass, and mid and end century delta biomass
# for a given RCP and grazing levels. Two grazing levels shown (top and bottom
# rows of maps). Delta biomass is calculated as scaled % difference from 
# current for that same grazing level.
# RCP 4.5 and 8.5 shown seperately

pdf("figures/biomass_maps/bio_current_vs_delta_by_graze_v1.pdf",
  width = 9, height = 6)


par( mar = c(1,1,2,1), mgp = c(3,0.3,0))
layout.matrix <- matrix(c(1,2,3,4,5,6),nrow = 2, ncol = 3, byrow = T)
layout(layout.matrix, widths = rep(1,3), heights = rep(1,2))


for (i in 1:nrow(names_plot1)) {
  row <- names_plot1[i, ]
  # current biomass
  if (row$type == "biomass" & row$RCP == "Current") {
    title <- paste("Current", row$PFT, "biomass,", row$graze, "grazing")
    
    image_bio(rast = rast1, subset = row$id, title = title)
    
  }
  
  # future delta biomass
  if (row$type == "bio-diff") {
    title <- substitute(paste(Delta, " ", PFT,  " (", RCP, " ", years,") ", 
                              graze, " grazing"),
                       list(PFT = as.character(row$PFT), 
                            years = as.character(row$years), 
                            graze = as.character(row$graze),
                            RCP = as.character(row$RCP)))

    image_bio_diff(rast = rast1, subset = row$id, title = title)
  }
 
}
dev.off()


# bio diff from light grazing -------------------------------------------------

# * current only ----------------------------------------------------------
# Scaled % change relative to current biomass under light grazing, for current
# conditions (e.g. % change going from light to moderate grazing, under
# current conditions)


# ** names of layers to compare -------------------------------------------

levs_pft <- levels(names_df$PFT)

names_current <- names_df %>% 
  filter(type == "biomass", RCP == "Current") %>% 
  arrange(PFT, graze) # arranging for figure creation

# information sets of layers that will be compared
lyrs_current <- map(levs_pft, function(pft) {
  
  out <- list()
  out$PFT <- pft
  df_pft <- names_current %>% 
    filter(PFT == pft)
  
  # name of reference layers
  out$ref_layer <- df_pft %>% 
    filter(graze == "Light") %>% 
    pull(id)
  
  stopifnot(length(out$ref_layer) == 1)
  
  target_rows <- df_pft %>% 
    filter(graze != "Light")
  
  out$graze <- target_rows$graze # grazing levels of non light grazing layers
  
  # ids of the non-light grazing layers
  out$target_layers <- target_rows$id
  
  out
})


# ** maps -----------------------------------------------------------------
# creates sets of 4 maps for each PFT
# Actual biomass (light grazing), and delta biomass for moderate,
# heavy and very heavy grazing

pdf("figures/biomass_maps/bio_current-light_and_delta-current_by_graze_v1.pdf",
    width = 9, height = 6)


# 'looping' over PFT's
map(lyrs_current, function(x) {
  par( mar = c(1,1,2,1), mgp = c(3,0.3,0))
  layout.matrix <- matrix(c(1,2,3,4,5,6),nrow = 2, ncol = 3, byrow = T)
  layout(layout.matrix, widths = rep(1,3), heights = rep(1,2))

  # current biomass, light grazing
  ctitle <- paste("Current", x$PFT, "biomass, light grazing")
  
  image_bio(rast = rast1, subset = x$ref_layer, title = ctitle)
    
  #  delta biomass under other grazing scenario ('looping' over grazing 
  # levels)
  map2(x$graze, x$target_layers, function(graze, target_layer) {
    title <- substitute(paste(Delta, " ", PFT,  ", current ", 
                               graze, " grazing"),
                         list(PFT = as.character(x$PFT), 
                              graze = as.character(graze)))
    
    # create raster of % scaled change
    r <- rast_diff(rast1, ref_layer = x$ref_layer, 
                   target_layer = target_layer)
    
    # r only has one layer
    image_bio_diff(rast = r, subset = 1, title = title)
  })

})

dev.off()
