# Purpose: Prepare interpolation products
# Prepare climate variables (based on daymet) for interpolation, and
# determine interpolation mask by finding the union of the 

library(terra)
library(sf)
library(tidyverse)
source("src/mapping_functions.R")
source("src/paths.R")

# params ------------------------------------------------------------------

sage_cutoffs <- c(1, 2, 5) # rcmap sage cover of 30 pixel (%), for it 
# to be 'counted' [this only applies to pixels that are within the SCD study area]

test_run <- FALSE # runs some code at low resolution for testing

# read in data ------------------------------------------------------------

# mean montly precip, and min and max temp, from daymet.
# calculated in the daymet_normals.js script, in the SEI repository
v <- 4 # daymet version
period <- '1991-2020' # period normals calculated for
met1 <- rast(paste0("data_raw/daymet_v", v, "_monthly_normals_",
                    period, ".tif"))
met1
crs(met1)

# file including site numbers (1:200) (for matching with simulation output)
site_nums <- read_csv("data_raw/site_locations.csv", show_col_types = FALSE)

# mean sage cover from RCMAP
# file created by SEI/scripts/export_RCMAP.js
# multiple tiles
paths_rcmap <- list.files(
  path = file.path(path_large, 
                   "grazing_effects/data_processed", 
                   "interpolation_data"),
  pattern = "RCMAPv6_sage-cover_mean_scd-extent_2001_2020_30m",
  full.names = TRUE
)

rcmap1 <- vrt(paths_rcmap)

if(test_run) {
  rcmap1 <- aggregate(rcmap1, fact = 10, fun = "mean")
}

# project/resample data ------------------------------------------------------------

# reclassify such that 1 means sagebrush habitat, 0 means masked out
scd_mask <- rcmap1
scd_mask <- classify(rcmap1, rcl = cbind(NA, 0))
scd_mask[!is.na(rcmap1)] <- 1
names(scd_mask) <- 'scd_all'
# plot(scd_mask)

rcmap_masks <- purrr::map(sage_cutoffs, function(x) {
  out <- scd_mask
  out[rcmap1 <= x] <- 0
  out
})
rcmap_masks <- rast(rcmap_masks)
names(rcmap_masks) <- paste0('sage_gt_', sage_cutoffs)

# resample to 1k resolution to match met data
masks_all <- c(scd_mask, rcmap_masks)

# this provides the proportion of cells in the 1km2 that belong to the
# given mask
met1 <- crop(met1, scd_mask)
prop_cells <- resample(masks_all, met1, method = 'average')
names(prop_cells) <- paste0("prop_", names(prop_cells))

# id raster ---------------------------------------------------------------
# Non NA grid cells contain the cell number

id1 <- prop_cells[['prop_scd_all']]
id1[] <- cells(id1)
id1 <- mask(id1, prop_cells[['prop_scd_all']], maskvalues = c(0, NA))
plot(id1)
names(id1) <- 'cellnumber'
# compute climate normals -------------------------------------------------
# bioclimatic variables based on avg monthly ppt and temp. (i.e. these
# are the worldclim variables)
met2 <- mask(met1, id1)

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

# scd_prop_cells <- crop(scd_prop_cells, bioclim2)
bioclim2 <- c(bioclim2, ptcor, prop_cells)

stopifnot(all(cells(bioclim2) == cells(id1)))

df_bioclim1 <- as.data.frame(bioclim2, cells = TRUE) %>% 
  rename(cellnumber = cell)

# merge cell number and site number info ----------------------------------

sf_site1 <- st_as_sf(site_nums[, c('site_id', "X_WGS84", "Y_WGS84")], 
                     coords = c("X_WGS84", "Y_WGS84"),
                     crs = "EPSG:4326") %>% 
  st_transform(crs = crs_scd)

# determining which cell each site is located in
site2 <- as.data.frame(terra::extract(id1, vect(sf_site1), bind = TRUE))

df_bioclim2 <- df_bioclim1 %>% 
  left_join(site2, by = 'cellnumber') %>% 
  filter(!is.na(bio1))

# save outputs ------------------------------------------------------------

if(!test_run) {
  write_csv(df_bioclim2, 
            paste0("data_processed/interpolation_data/clim_for_interp_", 
                   'daymet-v', v, '_', period, ".csv"))
  
  writeRaster(id1, "data_processed/interpolation_data/cellnumbers.tif",
              overwrite = TRUE)
}

