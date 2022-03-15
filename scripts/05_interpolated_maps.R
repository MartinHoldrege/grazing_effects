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
gref <- "Moderate" # current grazing level to compare to (for filter below)

# median biomass across GCMs
med2 <- rast(file.path(p, "bio_future_median_across_GCMs.tif"))

# ratio of C3Pgrass/total Pgrass, median across GCMs
gratio_med2 <- rast(file.path(p, "C3Pgrass-Pgrass-ratio_by-scenario_median.tif"))

# scaled % change of PGrass relative to current conditions of the same intensity
# (i.e. 'wgraze'). Median across GCMs
rast_d_Pgrass <- rast(file.path(p, "Pgrass_bio-diff-wgraze_median.tif"))

# minimum grazing level at which biomass threshold is exceeded
min_graze_files <- list.files(file.path(p, "min_graze/"),
                              full.names = TRUE)

rast_min_gr1 <- rast(min_graze_files)

# withing GCM scaled % change in biomass
# these are medians across GCMs
# e.g. change from heavy to light or moderate for RCP8.5 mid-century
wgcm_files <- list.files(file.path(p, "bio_diff"),
                         pattern = "_bio-diff-wgcm-heavy_",
                         full.names = TRUE)

rast_wgcm1 <- rast(wgcm_files)

# within grazing level scaled % change in biomass
# (i.e. current to future, for given grazing level), medians across GCMs
wgraze_files <- list.files(file.path(p, "bio_diff"),
                           pattern = "bio-diff-wgraze.*median.tif$",
                           full.names = TRUE)

rast_wgraze1 <- rast(wgraze_files)

# scaled % change from current moderate or heavy grazing to future (RCP8.5-mid) 
# light graze
# median across GCMs

rast_diff_gref <- rast(file.path(p, "bio-diff-gref-cur-heavy_median.tif"))

# * raster info -------------------------------------------------------------
# These files have been created in the 
# 04_interpolated_summarize.R script

# names of the raster layers, and treatment info, 
rast_info <- readRDS(file.path(p, "raster_info.RDS"))

# names/info of median grass layers (for c3pgrass)
grass_info_med <- readRDS(file.path(p, "grass_info_med.RDS"))

# names info about perennial grass layers
Pgrass_info_med <- readRDS(file.path(p, "Pgrass_info_med.RDS"))

# Maps ----------------------------------------------------------

# maps: wgcm & w-graze & gref----------------------------------------------
# change in biomass within a grazing level
# 9 panels for each PFT. 
# top row current: moderat grazing bio, heavy grazing, and delta from heavy to ligth
# middle row, same thing but for one future scenario. 
# bottom row, scaled percent change in biomass from current to RCP8.5 mid century,
# for moderate grazing, and second figure the same but for heavy grazing.
# 3rd figure (bottom right corner): change from current heavy grazing to future
# moderate

into_vars <- c("c4", "PFT", "type", "RCP", "years", "graze")
wgcm_info1 <- create_rast_info(
  rast_wgcm1, into = into_vars) %>% 
  # change from heavy grazing to gref (ie. light or moderate)
  filter(graze == gref)

wgraze_info1 <- create_rast_info(rast_wgraze1) %>%  
  select(-GCM) %>% # this is just 'median'
  arrange(PFT, desc(graze)) %>% 
  group_by(PFT) %>% 
  # making a higher numbered order variable so these plotted last on the page
  mutate(order = 1:n() + 6) 

gref_info <- create_rast_info(rast_diff_gref, into = into_vars) %>% 
  group_by(PFT) %>% 
  # making a higher numbered order variable so these plotted last on the page
  mutate(order = 9) %>% 
  filter(graze == gref)

# combing information on median biomass and wgcm and wgraze biomass difference
wgcm_info2 <- rast_info %>% 
  dplyr::select(-GCM,-id, -layer_num) %>% 
  rename(id = id_noGCM) %>% 
  filter_rcp_c4() %>% 
  filter(graze %in% c(gref, "Heavy")) %>% 
  distinct() %>% 
  bind_rows(wgcm_info1) %>% 
  # putting in order want to use when plotting
  arrange(PFT, RCP, desc(type), desc(graze)) %>% 
  group_by(PFT) %>% 
  mutate(order = 1:n()) %>% 
  bind_rows(wgraze_info1) %>% 
  bind_rows(gref_info) %>% 
  arrange(PFT, order)
  
wgcm_info2

# this figure creation takes a few minutes.

# NOTE that currently, some of the change values (heavy to light),
# are >100%. Those are throwing warnings (in the image_bio_diff function).

pdf("figures/biomass_maps/bio_wgcm-bio-diff_c4on_v5.pdf",
    width = wfig6, height = hfig6*3/2)

par(mar = mar, mgp = mgp)
layout(layout.matrix9, widths = widths9, heights = heights9)

pft_ids <- ""
#for (i in 1:9) { # for testing
for (i in 1:nrow(wgcm_info2)) {
 
  row <- wgcm_info2[i, ]
  RCP <- if (row$RCP == "Current") {
    as.character(row$RCP)
  } else {
    paste0(row$RCP, " (", row$years, ")")
  }
  
  if(row$type == "biomass") {
    
    # extracting min/max for all rasters of that PFT, so panels 
    # have the same range
    pft_ids_new <- wgcm_info2 %>% 
      filter(row$PFT == PFT, row$type == type) %>% pull(id)
    
    # so vector isn't being extracted repeatedly for the same PFT
    if(any(pft_ids_new != pft_ids)) {
      # very inefficient, but minmax() doesn't work when the raster contains
      # zero's which are plotted differently
      vec <- values(med2[[pft_ids_new]]) %>% as.vector()
    }
    
    
    # biomass figures
    title <- paste0(row$PFT, " ", RCP, ", ", row$graze, " grazing")
    
    image_bio(med2, subset = row$id, title = title,
              vec = vec)
    
    pft_ids <- pft_ids_new
    # effect of grazing within a climate scenario
  } else if(row$type == "bio-diff-wgcm-heavy") {
    # bio diff (within gcm) figure
    title <- substitute(paste(Delta, PFT, ", ", RCP, ", ", "heavy to ", 
                               graze, " grazing"),
                        list(c4 = row$c4,
                             PFT = as.character(row$PFT), 
                             graze = tolower(as.character(row$graze)),
                             RCP = RCP))
    image_bio_diff(rast_wgcm1, subset = row$id, title = title)
    # climate effect withing grazing trmt
  } else if (row$type == "bio-diff-wgraze") {
    title <- substitute(paste(Delta, PFT, ", Current to ", RCP, ", ", 
                              graze, " grazing"),
                        list(c4 = row$c4,
                             PFT = as.character(row$PFT), 
                             graze = tolower(as.character(row$graze)),
                             RCP = RCP))
    image_bio_diff(rast_wgraze1, subset = row$id, title = title)
    # change in climate and grazing
  } else {
    title <- substitute(paste(Delta, PFT, ", Current heavy to ", RCP, ", ", 
                              graze),
                        list(c4 = row$c4,
                             PFT = as.character(row$PFT), 
                             graze = tolower(as.character(row$graze)),
                             RCP = as.character(row$RCP)))
    image_bio_diff(rast_diff_gref, subset = row$id, title = title)
  }
}

dev.off()

# maps--min graze -------------------------------------------------------
#showing two panels for each PFT, the min graze for current, and min
# graze, for one future scenario.

min_gr_info <- tibble(id = names(rast_min_gr1),
                    id2 = id) %>% 
  separate(col = id2,
           into = c("c4", "PFT", "type", "RCP", "years"),
           sep = "_") %>% 
  df_factor() %>% 
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


pdf("figures/biomass_maps/C3-ratio_and_Pgrass_c4on-off_v2.pdf",
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


