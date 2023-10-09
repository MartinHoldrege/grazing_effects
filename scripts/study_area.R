# Map of study area for Kyle, with one map also showing the 200 points

# Author: Martin Holdrege

# Script Started October 9, 2023


# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)

old_wd <- getwd()
setwd("../SEI")
source("src/figure_functions.r") # basemap functions
setwd(old_wd)

# read in data ------------------------------------------------------------

# template raster. This is area that we interpolated data across (for
# some analyses this complete area isn't used and instead additional
# mask is added)
path_template <- "data_processed/interpolation_data/cellnumbers.tif"
r1 <- rast(path_template)

# locations of the simulation sites
site_nums <- read_csv("data_raw/site_locations.csv", show_col_types = FALSE)



# prepare site location data ----------------------------------------------

sf_site1 <- st_as_sf(site_nums[, c('site_id', "X_WGS84", "Y_WGS84")], 
                     coords = c("X_WGS84", "Y_WGS84"),
                     crs = "EPSG:4326") %>% 
  st_transform(crs = crs(r1))


# create maps -------------------------------------------------------------


r2 <- classify(r1, c(0, 1e7, 1))

g <- plot_map(stars::st_as_stars(r2), 
                st_geom_state = states,
                add_coords = TRUE,
         show_legend = FALSE) +
  labs(x = NULL,
       y = NULL) +
  scale_fill_manual(values = '#99d8c9',
                      na.value = 'transparent')

g2 <- g +
  geom_sf(data = sf_site1, aes(color = 'Simulation locations')) +
  coord_sf(xlim = bbox1[c(1, 3)],
           ylim = bbox1[c(2, 4)]) +
  scale_color_manual(values = 'black', name = NULL) +
  theme(legend.position = c(0.22, 0.07),
        legend.background = element_rect(fill = 'white'))

jpeg("figures/study_area/study_area_v1.jpeg", units = 'in',
     width = 4, height = 4, res = 600)
g
dev.off()

jpeg("figures/study_area/study_area_points_v1.jpeg", units = 'in',
     width = 4, height = 4, res = 600)
g2
dev.off()

