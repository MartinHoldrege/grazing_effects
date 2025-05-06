# Martin Holdrege

# Script started Jan 7, 2022

# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script, which can takes ~ 4 hours to run on
# my laptop (when excluding RCP4.5)

# this script takes ~ 10 min to run


# params ------------------------------------------------------------------

source("src/params.R") # run and version
test_run <- FALSE
print(Sys.time())
# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")
source("src/SEI_functions.R")

# params ------------------------------------------------------------------

v <- v_interp # interpolation version

runv <- paste0(run, v)

pfts <- c('Sagebrush', 'Pherb', 'Aherb')

ref_graze <- 'Moderate' # all delta covers are relative to this grazing level

# read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM
regex <- paste0(runv, '.*(', paste(pfts, collapse = '|'), ")")

regex <- "fire1_eind1_c4grass1_co20_2503v4.*(Sagebrush|Pherb|Aherb)"
bio_files <- list.files(file.path("data_processed/interpolated_rasters/biomass", v),
                        pattern = regex,
                        full.names = TRUE)

rast1 <- terra::rast(bio_files) # class SpatRast

scd_cov1 <- load_scd_cover() # cover's used for the SCD
scd_q1 <- load_scd_q()

r_eco1 <- load_wafwa_ecoregions_raster(wafwa_only = TRUE)
# info --------------------------------------------------

rast_info <- create_rast_info(rast1, id_noGCM = TRUE) %>% 
  filter_clim_extremes() %>% 
  arrange(PFT, graze, RCP, years, GCM, id)

rast2 <- rast1[[rast_info$id]]

if(test_run) {
  size <- 100
  rast2 <- spatSample(rast2, size = size, method = 'regular',
                      as.raster = TRUE)
  scd_cov1 <- spatSample(scd_cov1, size = size, method = 'regular',
                        as.raster = TRUE)
  r_eco1 <- spatSample(r_eco1, size = size, method = 'regular',
                       as.raster = TRUE)
  scd_q1 <- spatSample(scd_q1, size = size, method = 'regular',
                       as.raster = TRUE)
}

compareGeom(scd_cov1, rast2)

# convert stepwat biomass to cover ----------------------------------------

get_pft <- function(df) {
  unique(as.character(df$PFT))
}

# 'stepwat' cover
info_pft_l <- unname(split(rast_info, as.character(rast_info$PFT)))
sw_cov1 <-  info_pft_l  %>% 
  map(.f = function(df) {
    pft <- get_pft(df)
    bio2cov_hold(rast2[[df$id]], pft = pft)
  }) %>% 
  rast()

# proportional change cover -----------------------------------------------

# proportional cover relative to the reference state
sw_cov_prop1 <- map(info_pft_l, .f = function(df, cover = sw_cov1) {
    pft <- get_pft(df)
    ref_id <- df %>% 
      filter(RCP == 'Current',
             graze == ref_graze) %>% 
      pull(id)
    stopifnot(length(ref_id) == 1)
    r_ref <- cover[[ref_id]]
    r <- cover[[df$id]]
    (r - r_ref)/r_ref
  }) %>% 
  rast()

# scd cover adjusted by prop change ---------------------------------------

cov_scd_adj1 <- map(info_pft_l, .f = function(df) {
  pft <- get_pft(df)
  prop <- sw_cov_prop1[[df$id]]
  out <- scd_cov1[[pft]] + scd_cov1*prop
  names(out) <- names(prop)
  out
}) %>% 
  rast()

# Q scores of 'adjusted' cover -------------------------------------------------

q_scd_adj1 <- map(info_pft_l, .f = function(df, cover = cov_scd_adj1) {
  r <- cover[[df$id]]/100 # convert from % to proportion
  cov2q_raster(r, eco_raster = r_eco1, pft = get_pft(df))
}) %>% 
  rast()

# Calculate SEI -----------------------------------------------------------

# check that ordering is correct (so can multiply Qs)
regex_pft <-  paste(pfts, collapse = '|')
test <- map(info_pft_l, function(df) {
  
  str_replace(df$id, regex_pft, "")
})

# if this fails the next steps will be wrong (need same correct
# order of scenarios, grazing etc)
stopifnot(length(unique(test)) == 1,
          length(info_pft_l) == 3)

sei_scd_adj1 <- q_scd_adj1[[info_pft_l[[1]]$id]]*
  q_scd_adj1[[info_pft_l[[2]]$id]]*
  q_scd_adj1[[info_pft_l[[3]]$id]]*
  scd_q1[['Q4']]*
  scd_q1[['Q5']]

names(sei_scd_adj1) <- str_replace(names(sei_scd_adj1), regex_pft, 'SEI') %>% 
  str_replace('biomass', 'scd-adj')
  

# Calculate summaries and differences -------------------------------------

print(Sys.time())
