# Purpose: Take gridmet climatologies, and calculate bioclim variables for
# multivariate matching

# Author: Martin Holdrege

# Date started: June 13, 2024


# dependencies ------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
source('src/paths.R')


# read in data ------------------------------------------------------------

# *gridmet ----------------------------------------------------------------

# gridmet data downloaded from
# http://thredds.northwestknowledge.net:8080/thredds/catalog/MET/climatologies/monthlyClimatologies/catalog.html
period <- '19912020'
vars <- c('tmmn', 'tmmx', "ppt")
paths <- paste0(path_large,
                    "/climate_data/gridmet-clim/gridmet_", vars,
                    "_monthlyClimatologies_", 
                    period, ".nc")
names(paths) <- vars

# b/ of issue w/ terra::rast crashing with these files,
# reading into memory with raster first
grid1 <- map(paths, function(x) {
  raster::brick(x) %>% 
    raster::readAll() %>% 
    rast()
})


# * study area polygon ----------------------------------------------------
 
# Jeffries and Finn, 2019 polygon (https://doi.org/10.5066/P950H8HS)
biome1 <- sf::st_read('data_raw/US_Sagebrush_Biome_2019/US_Sagebrush_Biome_2019.shp')


# prepare layers ----------------------------------------------------------

# * biome -----------------------------------------------------------------

# extending biome extent, so that get gridmet data from a wider climate envelope
# so have better matching potential with daymet (e.g in hot/dry cells in the south);

biome2 <- biome1$geometry
buffer_km <- 0 # put buffer in km here
buffer_distance <- buffer_km*1000 # buffer in m
biome_buf1 <- st_buffer(biome2, dist = buffer_distance)
plot(biome_buf1)
plot(biome2, add = TRUE)

biome_buf2 <- st_transform(biome_buf1, crs = crs(grid1[[1]]))
biome_buf3 <- vect(biome_buf2) # to spatvector

# * gridmet ---------------------------------------------------------------

grid2 <- map(grid1, \(x) crop(x, biome_buf3))

grid3 <- map(grid2, \(x) mask(x, biome_buf3)) # make cells outside polygon NA

# convert K to C
grid3$tmmn <- grid3$tmmn - 273.15
grid3$tmmx <- grid3$tmmx - 273.15


# bioclim vars ------------------------------------------------------------

# calculate bioclim variables
grid4 <- map(grid3, raster::brick) # dismo package doesn't use terra objects
bioclim <- dismo::biovars(prec =grid4$ppt,
                          tmin = grid4$tmmn,
                          tmax = grid4$tmmx
)
bioclim2 <- rast(bioclim) # back to terra object


# output ------------------------------------------------------------------

writeRaster(bioclim2, 
            paste0('data_processed/interpolation_data/',
                   'gridmet_bioclim-vars_', period, '_', buffer_km, 'buffer.tif'))
