# Purpose:
# calculate five number summaries of Qs and SEI (remote sensing adjusted) values
# for creating boxplots and related figures for each ecoregion

# Author: Martin Holdrege
# Started: May 12, 2025

# params ------------------------------------------------------------------

source("src/params.R")
test_run <- FALSE
runv <- paste0(run, v_interp)

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')

# read in data ------------------------------------------------------------

 # file create in 05_interpolated_summarize_sei_scd-adj.R
r_qsei1 <- rast(file.path("data_processed/interpolated_rasters/", v_interp,
                     paste0(runv, "_q-sei_scd-adj_summary.tif")))

eco1 <- load_wafwa_ecoregions(wafwa_only = FALSE, total_region = TRUE)
eco2 <- vect(eco1['ecoregion'])

if(test_run) {
  size <- 100
  r_qsei1 <- spatSample(r_qsei1, size = size, method = 'regular',
                      as.raster = TRUE)
}

# names of layers --------------------------------------------------

into <- c("group", "type", "RCP", "years", "graze", 'summary')

info1 <- create_rast_info(r_qsei1, into = into) 


# calculate means by ecoregions -------------------------------------------

# taking  unweighted summary stats is ok, because we're using an equal area projection

# 5 number summary by ecoregion
eco_smry_l1 <- map(info1$id, function(id) {
  tmp <- terra::extract(r_qsei1[[id]], eco2, fivenum,  na.rm = TRUE)
  tmp <- as_tibble(tmp)
  colnames(tmp) <- c('ID', 'min', 'lower', 'middle', 'upper', 'max')
  tmp$'id' <- id
  tmp$region <- values(eco2)$ecoregion
  tmp$ID <- NULL
  tmp
})

eco_smry1 <- eco_smry_l1 %>% 
  bind_rows() %>% 
  left_join(info1, by = 'id')

# save output -------------------------------------------------------------

write_csv(eco_smry1, paste0('data_processed/raster_means/', runv, 
                            '_q-sei_scd-adj_summaries_by-ecoregion.csv'))
