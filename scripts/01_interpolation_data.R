# Purpose: Prepare interpolation products
# Prepare climate variables (based on daymet) for interpolation, and
# determine interpolation mask by finding the union of the 


library(terra)
library(sf)
library(tidyverse)
source("src/mapping_functions.R")

# read in data ------------------------------------------------------------

# mean montly precip, and min and max temp, from daymet.
# calculated in the daymet_normals.js script, in the SEI repository
met1 <- rast("data_raw/daymet_monthly_normals_1981-2010.tif")
met1
crs(met1)

# SCD snap raster from Geoff Bedrosian (to get SCD mask from)
scd_snap1 <- terra::rast("../SEI/data_processed/snap_raster/SEIv11_2017_2020_30_20211228_Q5sc3_Albers_2bit.tif")

# template raster Renne and palmquist used
palm1 <- rast("data_raw/interpolation_data_RR/CoreARTR_combined_DayMet_cropped_trimmed.tif")

# file including site numbers (1:200) (for matching with simulation output)
site_nums <- read_csv("data_raw/site_locations.csv", show_col_types = FALSE)

# project/resample data ------------------------------------------------------------

palm2 <- project(palm1, crs_scd)
palm_mask1 <- palm2
palm_mask1[] <- 1
palm_mask1 <- mask(palm_mask1, palm2, updatevalue = 0)
palm_mask2 <- resample(palm_mask1, met1, method = 'max')
  
compareGeom(met1, palm_mask2, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)

# reclassify such that 1 means sagebrush habitat, 0 means masked out
scd_mask <- classify(scd_snap1, rcl = rbind(c(2, 1), c(3, 1)))

# resample to 1k resolution to match met data

# this provides the proportion of cells in the 1km2 that are 'scd' cells
scd_prop_cells <- resample(scd_mask, met1,method = 'average')


# conservatively using 'max' method so any 30x30 pixel considered sagebrush
# makes the larger 1x1km pixel be considered sagebrush (so won't be missing
# interpolated data anywhere)
scd_mask1k <- resample(scd_mask, met1, 
                       method = 'max')
compareGeom(met1, scd_mask1k, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)


# create mask -------------------------------------------------------------

union_mask1 <- scd_mask1k+palm_mask2
union_mask1 <- classify(union_mask1, rbind(c(2, 1)))
plot(union_mask1)

# reduce extent, because lot's of unneeded NAs around the periphery
union_mask2 <- trim(union_mask1, padding = 10)
plot(union_mask2)


# id raster ---------------------------------------------------------------
# Non NA grid cells contain the cell number

id1 <- union_mask2
id1[] <- 1
id1[] <- cells(id1)
id1 <- mask(id1, union_mask2, maskvalues = c(0, NA))
plot(id1)

# compute climate normals -------------------------------------------------
# bioclimatic variables based on avg monthly ppt and temp. (i.e. these
# are the worldclim variables)
met2 <- crop(met1, id1) %>% 
  mask(id1)

prcp <- raster::brick(met2['prcp_.*'])
tmin <- raster::brick(met2['tmin_.*'])
tmax <- raster::brick(met2['tmax_.*'])

bioclim <- dismo::biovars(prec = prcp,
                          tmin = tmin,
                          tmax = tmax
                          )
bioclim2 <- rast(bioclim) # convert back to spatraster


# * calculate PTcor --------------------------------------------------------
# correlation between monthly precip and temperature (not a bio clim variable)

tmean <- (tmin + tmax)/2

# P-T correlation (type 1)
ptcor <- app(c(terra::rast(tmean), terra::rast(prcp)), fun = function(x) {
  # first 12 layers are temp, remains ones are prcp
  cor(x[1:12], x[13:24])
})

names(ptcor) <- 'ptcor'


# combine layers ----------------------------------------------------------
names(scd_prop_cells) <- 'prop_scd'
scd_prop_cells <- crop(scd_prop_cells, bioclim2)
bioclim2 <- c(bioclim2, ptcor, scd_prop_cells)

stopifnot(all(cells(bioclim2) == cells(id1)))

df_bioclim1 <- as.data.frame(bioclim2, cells = TRUE) %>% 
  rename(cellnumber = cell)


# merge cell number and site number info ----------------------------------

sf_site1 <- st_as_sf(site_nums[, c('site_id', "X_WGS84", "Y_WGS84")], 
                     coords = c("X_WGS84", "Y_WGS84"),
                     crs = "EPSG:4326") %>% 
  st_transform(crs = crs_scd)

# determining which cell each site is located in
site2 <- as.data.frame(terra::extract(id1, vect(sf_site1), bind = TRUE)) %>% 
  rename(cellnumber = SEIv11_2017_2020_30_20211228_Q5sc3_Albers_2bit)

df_bioclim2 <- df_bioclim1 %>% 
  left_join(site2, by = 'cellnumber') %>% 
  filter(!is.na(bio1))

# save outputs ------------------------------------------------------------

write_csv(df_bioclim2, "data_processed/interpolation_data/clim_for_interpolation.csv")

writeRaster(id1, "data_processed/interpolation_data/cellnumbers.tif",
            overwrite = TRUE)

