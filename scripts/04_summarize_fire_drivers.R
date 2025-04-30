# Purpose:
# calculate average of fire probability drivers (for each ecoregion) for each GCM
# downstream this will be used for helping understand the widespread
# fire probabilities across GCMs


# params ------------------------------------------------------------------

pfts <- c('Aherb', 'Pherb')
run <- 'fire1_eind1_c4grass1_co20_2503'
v_interp <- 'v4' # interpolation version (input files)

date <- '20250228' # date associated w/ the climate data (this date, in file name,
# would be updated if new climate data source, or sites, etc. were used. )

test_run <- FALSE
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')

# read in data ------------------------------------------------------------

# * climate ---------------------------------------------------------------

# data up-scaled for each GCM
clim_files <- list.files(
  file.path("data_processed/interpolated_rasters/climate", v_interp),
  pattern = paste0('climate.*', date, '_', v_interp, '.tif'),
  full.names = TRUE)
length(clim_files)

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
r_clim1 <- terra::rast(clim_files) # class SpatRast



# * biomass -----------------------------------------------------------------

bio_files <- list.files(
  file.path("data_processed/interpolated_rasters/biomass", v_interp),
  pattern = paste0(run, ".*", "(Pherb|Aherb)", '.*.tif'),
  full.names = TRUE)

r_bio1 <- rast(bio_files)

if (test_run) {
  r_clim1 <- r_clim1[[1:2]]
  r_bio1 <- r_bio1[[1:2]]
}

r_comb1 <- c(r_clim1, r_bio1)


# *ecoregions -------------------------------------------------------------

eco1 <- load_wafwa_ecoregions(total_region = TRUE)
eco2 <- vect(eco1)

# names of layers --------------------------------------------------

into_clim <- c("variable", "type", "RCP", "years", "GCM", 'date', 'version')

info_clim1 <- create_rast_info(r_clim1, into = into_clim, run_regex = '^') %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(id) %>% 
  select(-date, -run, -version)

info_bio1 <- create_rast_info(r_bio1) %>% 
  rename(variable = PFT) %>% 
  select(-run, -run2) # only 1 run, so not relevant

info_comb1 <- bind_rows(info_clim1, info_bio1)


# calculate means by ecoregions -------------------------------------------

# taking pixel means is ok, because we're using an equal area projection

# mean by ecoregion
eco_mean_l1 <- map(info_comb1$id, function(id) {
  terra::extract(r_comb1[[id]], vect(eco1), fun = mean,  na.rm = TRUE)
})

eco_mean1 <- map(eco_mean_l1, function(x) {
  x$ecoregion <- eco1$ecoregion[x$ID]
  x$ID <- NULL
  
  x %>% 
    pivot_longer(cols = -ecoregion,
                 names_to = 'id',
                 values_to = 'mean') 
  
}) %>% 
  bind_rows() %>% 
  left_join(info_comb1, by = 'id')

# save output -------------------------------------------------------------

write_csv(eco_mean1, paste0('data_processed/raster_means/', run, 
                            '_fire-driver-means_by-ecoregion.csv'))

