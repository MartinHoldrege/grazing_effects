# Purpose--summarize the main SCD cover layers to the same
# resolution (1km) and extent as used for this study

# Author: Martin Holdrege 

# Started May 5, 2025


# parameters --------------------------------------------------------------


# dependencies ------------------------------------------------------------

library(terra)
library(purrr)
source('src/paths.R')
source("src/mapping_functions.R")

# load data ---------------------------------------------------------------

# output from the Doherty 2022 OFR, these are tiles outputed from GEE by Dave. T
file_suffixes <- c("560m-0000000000-0000000000.tif",
                   "560m-0000000000-0000065536.tif")

file_pfts <- list('Sagebrush' = 'sage', 'Pherb' = 'perennialG', 'Aherb' = 'annualG')

vrt_l <- map2(file_pfts, names(file_pfts), function(pft, name) {
  paths <- file.path(path_large, "SEI_rasters/WAFWA30mdata", 
                     paste0("SEIv11_2017_2020_30_20211228_", pft, file_suffixes))
  
  r <- vrt(paths)
  names(r) <- name
  r
})

cover30m <- rast(vrt_l) # this is the high resolution (smoothed to 560m) cover

template <- load_template()


# aggregate the data ------------------------------------------------------

cover30mb <- project(cover30m, crs(template))
cover1 <- resample(cover30mb, template, method = 'average', threads = TRUE)

cover2 <- cover1
cover2[is.na(template)] <- NA


# save output -------------------------------------------------------------

writeRaster(cover2, 'data_processed/scd/SEIv11_2017_2020_1000_20211228_cover560.tif')





