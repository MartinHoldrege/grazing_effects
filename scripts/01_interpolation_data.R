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
rerun <- TRUE # recreate intermediate objects (used so some slow code doesn't need to rerun)
include_palm <- TRUE # should the study area include the full extent of previous
# studies (e.g. Palmquist et al 2021), makes the study area larger
# this doesn't apply to the output that 


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
  path <- 'data_processed/temp_rds/rcmap1_testrun.rds'
  
  if(file.exists(path) & !rerun) {
    rcmap <- readRDS(path)
  } else {
    rcmap1 <- aggregate(rcmap1, fact = 20, fun = "mean")
    saveRDS(rcmap1, path)
  }
  
}

# template raster Renne and palmquist used
palm1 <- rast("data_raw/interpolation_data_RR/CoreARTR_combined_DayMet_cropped_trimmed.tif")

# project/resample data ------------------------------------------------------------

path1 <- 'data_processed/temp_rds/prop_cells_testrun.rds'
path2 <- 'data_processed/temp_rds/prop_cells.rds'

if(!rerun & test_run & file.exists(path1)) {
  prop_cells <- readRDS(path1)
} else if (!rerun & !test_run & file.exists(path2)) {
  prop_cells <- readRDS(path2)
} else {
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
  
  if(test_run) {
    saveRDS(prop_cells, path1)
  } else {
    saveRDS(prop_cells, path2)
  }
}


# id raster ---------------------------------------------------------------
# Non NA grid cells contain the cell number

# palmquist et al. mask
palm2 <- project(palm1, crs_scd)
palm3 <- trim(resample(palm2, met1, method = 'max'))

if(include_palm) {
  prop_cells <- extend(prop_cells, palm3)
  tmp <- c(prop_cells[['prop_scd_all']], 
           extend(palm3,prop_cells[['prop_scd_all']]))
  # summing layers together, so nonzero values are can be considered
  # the study area (i.e. extending the study area to be union of two layers)
  m <- app(tmp, fun = function(x) sum(x, na.rm = TRUE))
} else {
  m <- prop_cells[['prop_scd_all']]
}

id1 <- m
id1[] <- cells(m)
id1 <- mask(id1, m, maskvalues = c(0, NA))
names(id1) <- 'cellnumber'

# compute climate normals -------------------------------------------------
# bioclimatic variables based on avg monthly ppt and temp. (i.e. these
# are the worldclim variables)


path <- "data_processed/temp_rds/bioclim2.rds"

pass <- FALSE
if(file.exists(path) & !rerun) {
  bioclim2 <- readRDS(path)
  pass <- compareGeom(bioclim2, id1, ext = TRUE,
                      stopOnError = FALSE)
} 

# run this section if rds doesn't exist, or even if it does 
# and it's read in and the extend is wrong (e.g. if changed code above
# to change study area)
if(!pass) {
  met2 <- mask(crop(met1, id1), id1)
  
  prcp <- raster::brick(met2['prcp_.*'])
  tmin <- raster::brick(met2['tmin_.*'])
  tmax <- raster::brick(met2['tmax_.*'])
  bioclim <- dismo::biovars(prec = prcp,
                            tmin = tmin,
                            tmax = tmax
  )
  bioclim2 <- rast(bioclim) # convert back to spatraster
  saveRDS(bioclim2, path)
}

# combine layers ----------------------------------------------------------

# scd_prop_cells <- crop(scd_prop_cells, bioclim2)
bioclim2 <- c(bioclim2, prop_cells)

stopifnot(all(cells(bioclim2[[1]]) == cells(id1)))

df_bioclim1 <- as.data.frame(bioclim2, cells = TRUE) %>% 
  rename(cellnumber = cell)

# merge cell number and site number info ----------------------------------

sf_site1 <- st_as_sf(site_nums[, c('site_id', "X_WGS84", "Y_WGS84")], 
                     coords = c("X_WGS84", "Y_WGS84"),
                     crs = "EPSG:4326") %>% 
  st_transform(crs = crs_scd)

test <- sum(is.na(terra::extract(bioclim2[[1]], sf_site1)[[2]]))
if(test > 0) warning('output raster doesnt cover all stepwat sites')
# determining which cell each site is located in
site2 <- as.data.frame(terra::extract(id1, vect(sf_site1), bind = TRUE))

df_bioclim2 <- df_bioclim1 %>% 
  left_join(site2, by = 'cellnumber') %>% 
  filter(!is.na(bio1))

summary(df_bioclim1$prop_scd_all)


# only cell numbers where proportion scd cells >50% -----------------------

id_scd50 <- id1
tmp <- fill_raster(df_bioclim2[, c('cellnumber', 'prop_scd_all')], id1)
id_scd50[tmp < 0.5] <- NA

# save outputs ------------------------------------------------------------

if(!test_run) {
  write_csv(df_bioclim2, 
            paste0("data_processed/interpolation_data/clim_for_interp_", 
                   'daymet-v', v, '_', period, ".csv"))
  
  writeRaster(id1, "data_processed/interpolation_data/cellnumbers.tif",
              overwrite = TRUE)
  
  writeRaster(id_scd50, "data_processed/interpolation_data/cellnumbers_scd50.tif",
              overwrite = TRUE)
}

