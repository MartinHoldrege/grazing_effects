# Purpose: create map of study area showing ecoregions, these are the regions
# some the data analyses use for data summarization


# params ------------------------------------------------------------------

v <- 'v4' # interpolation version


# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source('src/mapping_functions.R')
source('src/fig_params.R')


# read in data ------------------------------------------------------------

match1 <- rast(paste0("data_processed/interpolation_quality/",
                      "matching_quality", v, '.tif'))

eco1 <-load_wafwa_ecoregions()

# raster of ecoregions ----------------------------------------------------

study_area <- match1
study_area[study_area > 1.5] <- NA
study_area[!is.na(study_area)] <- 1

r_eco <- terra::rasterize(eco1, study_area, field = "ecoregion")
r_eco[is.na(study_area)] <- NA

# map ---------------------------------------------------------------------

g <- plot_map2(r_eco) +
  scale_fill_manual(values = cols_ecoregion,
                    na.value = 'transparent',
                    na.translate = FALSE,
                    name = NULL) +
  theme(legend.position = 'bottom')

png(paste0("figures/ecoregions_", v, '.png'),
    width = 6, height = 7, units = 'in', res = 900)
g
dev.off()
