# climate and other stats useful for the study area description

# dependencies ------------------------------------------------------------

source('src/mapping_functions.R')
source('src/general_functions.R')
library(tidyverse)
library(terra)

# read in data ------------------------------------------------------------

# study area (for masking)
size <- load_cell_size()

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

# Doherty et al. study area (for size comparison)
c3 <- rast('../SEI/data_publication1/rasters/SEIv11_2017_2020_30_Current_20220717.tif')

# mask --------------------------------------------------------------------

cellnums2 <- crop(cellnums, size)
compareGeom(cellnums2, size)
cellnums2 <- mask(cellnums2, size)

# get MAP & MAT for study area
df_cellnums <- data.frame(cellnums2)
clim2 <- clim1 %>% 
  filter(cellnumber %in% df_cellnums$cellnumber)

# get psp for study area 
psp2 <- crop(psp1, size)
psp3 <- mask(psp2, size)

clim3 <- c(clim2, psp3)

# size of study area ------------------------------------------------------

# total area
tot_area <- global(size, fun = 'sum', na.rm = TRUE)
tot_area
# sum
# area 90743900

global(psp3, 'notNA') # the psp3 data is from a more sparse that covers the whole
# region b/ >25% of grid-cells are missing 
global(size, 'notNA')

c3b <- ifel(c3 == 0, NA, c3)
area_scd <- cellSize(c3b, mask = TRUE,
                      unit = 'ha',
                      transform = FALSE)
tot_area_scd <- global(area_scd , fun = 'sum', na.rm = TRUE)
# sum
# area 102690249
100 - tot_area/tot_area_scd*100
#  11.63338

# climate -----------------------------------------------------------------

minmax(psp3)
# prcpProp
# min 0.03686094
# max 0.49768746
global(psp3, fun = 'mean', na.rm = TRUE)
# mean
# prcpProp 0.2099129
global(psp3, fun = median, na.rm = TRUE)
# global
# prcpProp 0.1946597

clim2 %>% 
  select(bio1, bio12) %>% 
  rename(MAT = 'bio1',
         MAP = 'bio12') %>% 
  summary()
# MAT                MAP       
# Min.   : 0.07887   Min.   :106.9  
# 1st Qu.: 6.83853   1st Qu.:282.7  
# Median : 8.02092   Median :356.6  
# Mean   : 8.22280   Mean   :357.7  
# 3rd Qu.: 9.74546   3rd Qu.:417.8  
# Max.   :15.11510   Max.   :958.7  