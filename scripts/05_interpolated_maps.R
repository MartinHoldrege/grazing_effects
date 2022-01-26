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
# bio_d_files <- list.files("data_processed/interpolated_rasters/bio_diff/",
#                         pattern = ".tif$",
#                         full.names = TRUE)


# for now just using a subset for computational efficiency
bio_files <- bio_files[sample(1:length(bio_files), 200)]

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(bio_files) # class SpatRast


# * AIM -------------------------------------------------------------------

# 2011-2015 AIM vegetation data that Rachel Renne has compiled
aim1 <- read_csv("data_raw/AIM/LMF_ALL.data.2011-2015.FINAL_v2.csv")


# params ------------------------------------------------------------------

# parameters used for layout
mar <- c(1,1,2,1)
mgp = c(3,0.3,0)

# parameters specific to 6 paneled figures
layout.matrix6 <- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, 
                         byrow = TRUE) # for 6 paneled figs
widths6 <- rep(1, 3)
heights6 = rep(1, 2)

wfig6 <- 9 # width of 9 figure panels (inches)
hfig6 <- 6 # height of 9 figure panels

# some terra operations can be done in parallel and need 
# to know num of cores
num.cores <- parallel::detectCores(logical = FALSE) 

# names of layers --------------------------------------------------
# descriptions of layers of raster stack
r_names <- names(rast1)

names_df <- tibble(id2 = r_names,
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


# median by GCM -----------------------------------------------------------
# for each treatment/scenario combination calculate the median
# across GCMs



# create SpatRasterDataset

rast2 <- rast1[[names_df$id]] # make sure rows are ordered correctly

# list of rasters, each element is a SpatRaster including all GCMs
# for a given treatment/scenario combination
rast_gcm_l <- split(rast2, f = names_df$id_noGCM)

names(rast_gcm_l) <- unique(names_df$id_noGCM)

# SpatRasterDataset (i.e. each dataset includes all GCMs for treatment/scenario)
sds_gcm <- sds(rast_gcm_l) 
class(sds_gcm)

# check that the number of layers in each dataset is correct
stopifnot(nlyr(sds_gcm) == names_df %>% 
            group_by(id_noGCM) %>% 
            summarise(n = n()) %>% 
            pull(n)
          )

# compute median across GCMs
system.time(
  sds_med1 <- terra::app(sds_gcm, f = "median", cores = num.cores)
)



# NEXT--take medians across GCM's in some kind of sensible way
# First make them into a spatrasterdataset (sds), then use app(x, median), 
# to take the median of each dataset
# m



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
par(mar = mar, mgp = mgp)
layout(layout.matrix6, widths = widths6, heights = heights6)


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

# pdf("figures/biomass_maps/bio_current_vs_delta_by_graze_v1.pdf",
#   width = wfig6, height = hfig6)


par(mar = mar, mgp = mgp)
layout(layout.matrix6, widths = widths6, heights = heights6)

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

# pdf("figures/biomass_maps/bio_current-light_and_delta-current_by_graze_v1.pdf",
#     width = wfig6, height = hfig6)


# 'looping' over PFT's
map(lyrs_current, function(x) {
  par(mar = mar, mgp = mgp)
  layout(layout.matrix6, widths = widths6, heights = heights6)

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


# AIM data maps -----------------------------------------------------------
# map(s) showing AIM C4Pgrass occurrence on top of interpolated STEPWAT2 
# C4Pgrass biomass data

# * process ---------------------------------------------------------------

# the WBG and WRG columns are c4 perennial bunch grasses and c4 perennial 
# rhizomatous grasses, respectively. 

aim2 <- aim1 %>% 
  select(PLOTKEY, latitude, longitude, WBG, WRG) %>% 
  # are C4 P grasses present?
  mutate(C4Pgrass = ifelse(WBG > 0 | WRG > 0, "present", "absent")) 

mean(aim2$C4Pgrass == "present") * 100 # % of sites w/ c4pgrass present


# creating a SpatVector (terra package)
# for some reason logical data column gets converted to character,
# C4Pgrass needed to be numeric
aim_sv1 <- vect(aim2, geom = c("longitude", "latitude"))


# * maps ------------------------------------------------------------------

# Next--separate maps showing the 200 sites, and c4 presence/absence
# pdf("figures/biomass_maps/C4Pgrass-bio_AIMC4-Pgrass-presence.pdf",
#     width = wfig6, height = hfig6)

pch1 <- 16
cex1 <- 0.6

par(mar = mar, mgp = mgp)
layout(layout.matrix6, widths = widths6, heights = heights6)

c4layer <- "C4Pgrass_biomass_Current_Current_Light" # Layer to plot
image_bio(rast1, subset = c4layer, title = "Current C4Pgrass, light grazing")

image_bio(rast1, subset = c4layer, title = "AIM plots where C4Pgrass absent (black points)")
plot(subset(aim_sv1, subset = aim_sv1$C4Pgrass == "absent"), 
     col = cols_c4present["absent"], 
     add = TRUE, pch = pch1, cex = cex1)

image_bio(rast1, subset = c4layer, title = "AIM plots where C4Pgrass present (blue points)")
plot(subset(aim_sv1, subset = aim_sv1$C4Pgrass == "present"), 
     col = cols_c4present["present"], 
     add = TRUE, pch = pch1, cex = cex1)

image_bio(rast1, subset = c4layer, title = "All AIM plots")
plot(aim_sv1, 
     col = cols_c4present[aim_sv1$C4Pgrass], 
     add = TRUE, pch = pch1, cex = cex1)

dev.off()
