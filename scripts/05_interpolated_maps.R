# Martin Holdrege

# Script started Jan 7, 2022

# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script which can takes ~ 4 hours to run on
# my laptop (when excluding RCP4.5) and rasters created in the
# the 04_inperpolated_summarize.R script, 

# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")

# read in data ------------------------------------------------------------

p <- "data_processed/interpolated_rasters"

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM, and for c4on and c4off
bio_files <- list.files(file.path(p, "biomass/"),
                        pattern = "^c4.*tif$",
                        full.names = TRUE)
length(bio_files)

# for now just using a subset for computational efficiency
#bio_files <- bio_files[sample(1:length(bio_files), 200)]

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(bio_files) # class SpatRast

# median biomass across GCMs
med1 <- terra::rast(file.path(p, "bio_future_median_across_GCMs.tif"))

# minimum grazing level at which biomass threshold is exceeded
min_graze_files <- list.files(file.path(p, "min_graze/"),
                              full.names = TRUE)

rast_min_gr1 <- rast(min_graze_files)

# * raster info -------------------------------------------------------------

# names of the raster layers, and treatment info, created in the 
# 04_interpolated_summarize.R script
rast_info <- readRDS(file.path(p, "raster_info.RDS"))


# prepare rasters ----------------------------------------------------------

# * median biomass ----------------------------------------------------------

# layers of current biomass
info_current <- rast_info %>% 
  filter(RCP == "Current") 

rast_current <- rast1[[info_current$id]]

# renaming without trailing 'current'
names(rast_current) <- info_current$id_noGCM 

# combing so have median across gcms in the future, and current biomass
med2 <- c(rast_current, med1) # med1 didn't include current


# * c3Pgrass/Pgrass -------------------------------------------------------
# this might belong in the 04_...summarize script. 

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

# grass ratio (this is throwing a warning--not sure why)
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


# maps--min graze -------------------------------------------------------

min_gr_info <- tibble(id = names(rast_min_gr1),
                    id2 = id) %>% 
  separate(col = id2,
           into = c("c4", "PFT", "type", "RCP", "years"),
           sep = "_") %>% 
  mutate(years = years2factor(years),
         RCP = rcp2factor(RCP),
         PFT = pft_all_factor(PFT)) %>% 
  arrange(PFT, RCP)

pdf("figures/min_graze_maps/min_graze_c4on.pdf",
    width = wfig6, height = hfig6)

par(mar = mar, mgp = mgp)
layout(layout.matrix4, widths = widths6, heights = heights6)

for (i in 1:nrow(min_gr_info)) {
  row <- min_gr_info[i, ]
  
  RCP <- if (row$RCP == "Current") {
    row$RCP
  } else {
    paste0(row$RCP, " (", row$years, ")")
  }
  title <- paste(row$PFT, RCP)
  image_min_gr(rast_min_gr1, subset = row$id,
               title = title)
}

dev.off()

# maps--Pgrass and C3Pgrass/Pgrass ------------------------------------------

# testing
# either pass the name of the layer or the layer number
image_bio(med2, subset = 1, title = "test",
          vec = c(0, 0.1, 1),
          show0legend = FALSE,
          legend_lab = "C3Pgrass/Pgrass")


pdf("figures/biomass_maps/C3-ratio_and_Pgrass_c4on-off.pdf",
  width = wfig6, height = hfig6)

# C3Pgrass/Pgrass maps (not showing % change just the ratio)

## Set parameters and layout:
par(mar = mar, mgp = mgp)
layout(layout.matrix6, widths = widths6, heights = heights6)

for(i in 1:nrow(grass_info_med)) {
  row <- grass_info_med[i, ]
  if(row$RCP == "Current") {
    title = paste0("Current, ", tolower(row$graze), " grazing")
  } else {
    title = paste0(row$c4, " (", row$RCP, " ",
                   row$years, "), ", tolower(row$graze), " grazing")
  }
  
  image_bio(gratio_med2, subset = row$id_noGCM, title = title,
            vec = c(0, 0.1, 1),
            show0legend = FALSE,
            legend_lab = "C3Pgrass/Pgrass",
            n_breaks = 31)
}

# Pgrass maps, showing current, and percent change from current
# for that grazing level for c4 on and c4 off

for(i in 1:nrow(Pgrass_info_med)) {
  row <- Pgrass_info_med[i, ]

  if(row$RCP == "Current") {
    title = paste0("Current Pgrass, ", tolower(row$graze), " grazing")
    
    image_bio(med2, subset = row$id_noGCM, title = title)
  } else {
    title <- substitute(paste(c4," ", Delta, " Pgrass (", RCP, " ", years,") ", 
                              graze, " grazing"),
                        list(c4 = row$c4,
                             PFT = as.character(row$PFT), 
                             years = as.character(row$years), 
                             graze = as.character(row$graze),
                             RCP = as.character(row$RCP)))
    
    image_bio_diff(rast_d_Pgrass, subset = row$id_noGCM, title = title)
  }
}

dev.off()


