# climate and other stats useful for the study area description

# dependencies ------------------------------------------------------------

source('src/mapping_functions.R')
source('src/general_functions.R')
library(tidyverse)
library(terra)

# read in data ------------------------------------------------------------

# study area (for masking)
r_eco1 <- load_wafwa_ecoregions_raster() # the study area

# cell numbers corresponding to the climate raster
# file created in 01_interpolation_data.R
cellnums <- rast("data_processed/interpolation_data/cellnumbers.tif")

# Read in raw bioclim data for sagebrush extent and set cellnumbers as rownames
clim1 <- read_csv("data_processed/interpolation_data/clim_for_interp_daymet-v4_1991-2020.csv",
                col_types = cols(site_id = 'd'),
                show_col_types = FALSE)

# proportion summer precipitation data from the Holdrege et al 2024 Fire Ecology manuscript
# This is for estimating the mean and range of proportion summer ppt
# across the study area
psp1 <- terra::rast('../cheatgrass_fire/data_processed/daymet/daymet_climSummerAvg_1986-2019_1000m_sagebrush-biome-mask_v2.tif')[['prcpProp']]


# mask --------------------------------------------------------------------

cellnums2 <- crop(cellnums, r_eco1)
compareGeom(cellnums2, r_eco1)
cellnums2 <- mask(cellnums2, r_eco1)

# get MAP & MAT for study area
df_cellnums <- data.frame(cellnums2)
clim2 <- clim1 %>% 
  filter(cellnumber %in% df_cellnums$cellnumber)

# get psp for study area 
psp2 <- crop(psp1, r_eco1)
psp3 <- mask(psp2, r_eco1)
