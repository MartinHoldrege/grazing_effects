# Martin Holdrege

# Script started Jan 27, 2022

# Purpose: to plot AIM data showing observed presence/absence of C4 perennial
# grasses, and compare that to where stepwat2 is simulating c4pgrass when that
# data is upscaled. 

# Additional analyses will also be done to look at AIM sites that show C4grass
# but for which the upscaled data doesn't simulate c4 grass, and then look
# at what c4 grass species occur there. Also just tally for how mainy aim sites
# c4 grass presence/absence is mis-specified

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source("src/mapping_functions.R")

# Read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------
# upscaled stepwat2 output (created in the 03_interpolate.R script),
# for c4 grasses under current conditions, light grazing.

C4_r <- rast(file.path("data_processed/interpolated_rasters/biomass",
               "c4on_C4Pgrass_biomass_Current_Current_Light_Current.tif"))

# * AIM -------------------------------------------------------------------

# 2011-2015 AIM vegetation data that Rachel Renne has compiled
aim1 <- read_csv("data_raw/AIM/LMF_ALL.data.2011-2015.FINAL_v2.csv")


# * site locations --------------------------------------------------------
# location of the 200 sites for which simulations were run for 
# and whether c4pgrass is simulated to be present there under current 
# conditions

# file created in "03_interpolate.R" script
sites_c4 <- read_csv(
  "data_processed/site-num_C4Pgrass-presence_c4off.csv")

# process AIM data --------------------------------------------------------

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


# process sites -----------------------------------------------------------

# creating SpatVector
sites_c4_sv1 <- vect(sites_c4, geom = c("x", "y"))

# * maps ------------------------------------------------------------------

# map(s) showing AIM C4Pgrass occurrence on top of interpolated STEPWAT2 
# C4Pgrass biomass data, as well as the 200 simulation sites

# base of the map for several panels
base_map <- function(title) {
  image_bio(C4_r, subset = 1, title = title)
}

# due to some unknown issue this throws an error
# the first time it is run (when image resolution is high), so calling
# function here where it can't cause harm
# (see https://github.com/rspatial/terra/issues/30)
safely(base_map(title = "Current C4Pgrass, light grazing"))

pdf("figures/biomass_maps/C4Pgrass-bio_AIMC4-Pgrass-presence.pdf",
    width = wfig6, height = hfig6)

pch1 <- 16
cex1 <- 0.4

par(mar = mar, mgp = mgp)
layout(layout.matrix6, widths = widths6, heights = heights6)

base_map(title = "Current C4Pgrass, light grazing")

# AIM plots on top
base_map(title = "AIM plots where C4Pgrass absent (black points)")
plot(subset(aim_sv1, subset = aim_sv1$C4Pgrass == "absent"), 
     col = cols_c4present["absent"], 
     add = TRUE, pch = pch1, cex = cex1)

base_map(title = "AIM plots where C4Pgrass present (blue points)")
plot(subset(aim_sv1, subset = aim_sv1$C4Pgrass == "present"), 
     col = cols_c4present["present"], 
     add = TRUE, pch = pch1, cex = cex1)

base_map(title = "All AIM plots")
plot(aim_sv1, 
     col = cols_c4present[aim_sv1$C4Pgrass], 
     add = TRUE, pch = pch1, cex = cex1/2)

# showing 200 sites
base_map(title = "200 simulation sites (blue = C4Pgrass presence")
plot(sites_c4_sv1, 
     col = cols_c4present[sites_c4_sv1$C4Pgrass], 
     add = TRUE, pch = pch1, cex = 0.8)

dev.off()
