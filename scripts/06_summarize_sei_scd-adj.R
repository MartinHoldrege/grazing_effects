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
source("src/SEI_functions.R")

# read in data ------------------------------------------------------------

 # file create in 05_interpolated_summarize_sei_scd-adj.R
r_qsei1 <- rast(file.path("data_processed/interpolated_rasters/", v_interp,
                     paste0(runv, "_q-sei_scd-adj_summary.tif")))

eco1 <- load_wafwa_ecoregions(wafwa_only = FALSE, total_region = TRUE)
eco2 <- vect(eco1['ecoregion'])

if(test_run) {
  r_qsei1 <- downsample(r_qsei1)
}

# names of layers --------------------------------------------------

into <- c("group", "type", "RCP", "years", "graze", 'summary')

info1 <- create_rast_info(r_qsei1, into = into) 

# calculate means by ecoregions -------------------------------------------

# taking  unweighted summary stats is ok, because we're using an equal area projection
sixnum <- function(x, na.rm = TRUE) {
  five <- fivenum(x, na.rm = na.rm)
  c(five, mean(x, na.rm = na.rm))
}
# 5 number summary by ecoregion
eco_smry_l1 <- map(info1$id, function(id) {
  tmp <- terra::extract(r_qsei1[[id]], eco2, sixnum,  na.rm = TRUE)
  tmp <- as_tibble(tmp)
  colnames(tmp) <- c('ID', 'min', 'lower', 'middle', 'upper', 'max', 'mean')
  tmp$'id' <- id
  tmp$region <- values(eco2)$ecoregion
  tmp$ID <- NULL
  tmp
})

eco_smry1 <- eco_smry_l1 %>% 
  bind_rows() %>% 
  left_join(info1, by = 'id')

# calculating % CSA & GOA ------------------------------------------------------

r_sei <- r_qsei1[[info1$id[info1$type == 'SEI']]]

# % core
sei_core1 <- terra::extract(r_sei, eco1, percent_csa,  na.rm = TRUE)
sei_core2 <- pivot_longer_extracted(sei_core1, as.character(eco1$ecoregion), 
                                    values_to = 'percent_csa') 
# % growth opportunity area
sei_goa1 <- terra::extract(r_sei, eco1, percent_goa,  na.rm = TRUE)
sei_goa2 <- pivot_longer_extracted(sei_goa1, as.character(eco1$ecoregion), 
                                    values_to = 'percent_goa') 

sei_pcent <- left_join(sei_core2, sei_goa2, by = c('region', 'id')) %>% 
  left_join(info1, by = 'id') %>% 
  select(-run2)

# save output -------------------------------------------------------------

prefix <- if(test_run) 'test' else runv

write_csv(eco_smry1, paste0('data_processed/raster_means/', prefix, 
                            '_q-sei_scd-adj_summaries_by-ecoregion.csv'))

write_csv(sei_pcent , paste0('data_processed/raster_means/', prefix, 
                            '_sei-class-pcent_scd-adj_summaries_by-ecoregion.csv'))
