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
# intermediate rasters saved, usually want true b/ of disc storage, except when testing/debugging
remove_tmp_files = FALSE
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

# summary stats calculated for these regions:
eco1 <- load_wafwa_ecoregions(total_region = TRUE, wafwa_only = FALSE)
r_eco4 <- load_wafwa_ecoregions_raster()
size <- load_cell_size()
# info --------------------------------------------------

rast_info <- create_rast_info(rast1, id_noGCM = TRUE) %>% 
  filter_clim_extremes() %>% 
  arrange(PFT, graze, RCP, years, GCM, id)

rast2 <- rast1[[rast_info$id]]

if(test_run) {
  rast2 <- downsample(rast2)
  scd_cov1 <- downsample(scd_cov1)
  r_eco1 <- downsample(r_eco1)
  scd_q1 <- downsample(scd_q1)
  r_eco4 <- downsample(r_eco4)
  # after downsampling the size is wrong, but this is just for testing
  size <- downsample(size) 
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
rm('sw_cov1')
gc()
# scd cover adjusted by prop change ---------------------------------------

cov_scd_adj1 <- map(info_pft_l, .f = function(df) {
  pft <- get_pft(df)
  prop <- sw_cov_prop1[[df$id]]
  out <- scd_cov1[[pft]] + scd_cov1[[pft]]*prop
  names(out) <- names(prop)
  out
}) %>% 
  rast()

# writing to disk so not stored in memory
dir.create('tmp/', showWarnings = FALSE)
rm('sw_cov_prop1', 'scd_cov1')
gc()
cov_scd_adj1 <- writeReadRast(cov_scd_adj1, 'cov_scd_adj1')


# Q scores of 'adjusted' cover -------------------------------------------------

q_scd_adj1 <- map(info_pft_l, .f = function(df, cover = cov_scd_adj1) {
  r <- cover[[df$id]]/100 # convert from % to proportion
  cov2q_raster(r, eco_raster = r_eco1, pft = get_pft(df))
}) %>% 
  rast()
q_scd_adj1 <- writeReadRast(q_scd_adj1, 'q_scd_adj1')

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

# leaving same number of name components as the Q layers, so can combine
names(sei_scd_adj1) <- str_replace(names(sei_scd_adj1), regex_pft, 'SEI') %>% 
  str_replace('biomass', 'SEI')

sei_scd_adj1 <- writeReadRast(sei_scd_adj1, 'sei_scd_adj1')



# * calculate gcm level summaries -----------------------------------------

info_sei <- create_rast_info(sei_scd_adj1, 
                             into = c("group", "type", "RCP", "years", 
                                      "graze", "GCM"), 
                             id_noGCM = TRUE)  




# mean SEI for each region and GCM
sei_mean1 <- terra::extract(sei_scd_adj1, eco1, mean,  na.rm = TRUE)
sei_core1 <- terra::extract(sei_scd_adj1, eco1, percent_csa,  na.rm = TRUE)
sei_goa1 <- terra::extract(sei_scd_adj1, eco1, percent_goa,  na.rm = TRUE)
regions <- as.character(eco1$ecoregion)
sei_mean2 <- pivot_longer_extracted(sei_mean1, regions, 
                                    values_to = 'SEI_mean')

sei_core2 <- pivot_longer_extracted(sei_core1, regions, 
                                    values_to = 'percent_csa')
sei_goa2 <- pivot_longer_extracted(sei_goa1, regions, 
                                   values_to = 'percent_goa')

sei_mean_core1 <- left_join(sei_mean2, sei_core2, by = c("id", "region")) %>% 
  left_join(info_sei, by = 'id') %>% 
  left_join(sei_goa2, by = c("id", "region")) %>% 
  select(-run2, -id_noGCM)


# * gcm-wise mean sei by c3 and ecoregion ---------------------------------
# mean SEI for each ecoregion and historical SEI class in that ecoregion

# 'SEI of the reference class (i.e. this is SEI as in Doherty, except for
# the artifacts of our approach that are do to coarse 1km scaling)
id_ref <- info_sei %>% 
  filter(RCP == 'Current', group == 'SEI', graze == ref_graze) %>% 
  pull(id)
stopifnot(length(id_ref) == 1)

c3_ref <- sei2c3(sei_scd_adj1[[id_ref]])
c3eco <- c3_ref*10 + as.numeric(r_eco4)
names(c3eco) <- 'c3eco'
c3eco_sei1 <- zonal(sei_scd_adj1, c3eco, fun = 'mean')
c3eco_area1 <- zonal(size, c3eco, fun = 'sum')

c3eco_sei2_gcm <- c3eco_sei1 %>% 
  pivot_longer(-c3eco,
               values_to = 'SEI_mean',
               names_to = 'id') %>% 
  mutate(c3 = c3eco_to_c3(c3eco),
         region = c3eco_to_eco(c3eco, levels(r_eco4)[[1]]$ecoregion)) %>% 
  left_join(c3eco_area1, by = "c3eco") %>% 
  left_join(info_sei, by = 'id') %>% 
  select(-c3eco)


# Calculate summaries and differences -------------------------------------


# *median/low/high --------------------------------------------------------
# pixelwise summaries across GCMs

cov_scd_adj_smry <- summarize_gcms_raster(r = cov_scd_adj1, info = rast_info,
                                          include_low_high = FALSE)
rm('cov_scd_adj1')

into_smry <- c("group", "type", "RCP", "years", "graze", "summary")
info_smry <- create_rast_info(cov_scd_adj_smry, into = into_smry)

q_scd_adj_smry <- summarize_gcms_raster(q_scd_adj1, info = rast_info,
                                          include_low_high = TRUE)
rm('q_scd_adj1'); gc()


sei_scd_adj_smry <- summarize_gcms_raster(sei_scd_adj1, info = info_sei,
                                          include_low_high = TRUE)

info_sei_smry <- create_rast_info(sei_scd_adj_smry, 
                 into = c("group", "type", "RCP", "years",  "graze", "summary"))

# * pixelwise summary mean sei by c3 and ecoregion ---------------------------------
# mean SEI for each ecoregion and historical SEI class in that ecoregion

c3eco_sei1_smry <- zonal(sei_scd_adj_smry, c3eco, fun = 'mean')

c3eco_sei2_smry <- c3eco_sei1_smry %>% 
  pivot_longer(-c3eco,
               values_to = 'SEI_mean',
               names_to = 'id') %>% 
  mutate(c3 = c3eco_to_c3(c3eco),
         region = c3eco_to_eco(c3eco, levels(r_eco4)[[1]]$ecoregion)) %>% 
  left_join(c3eco_area1, by = "c3eco") %>% 
  left_join(info_sei_smry, by = 'id') %>% 
  select(-c3eco)

# * differences (within grazing level) relative to current climate -----------
rm('sei_scd_adj1') # make more memory available

by <- c('run', 'group', 'type', 'graze') 
cov_scd_adj_diff <- calc_rast_cref(cov_scd_adj_smry, info_smry,
                                   by = by,
                                   type_from = '_biomass',
                                   type_to = '_cov-rdiff-cref')

q_scd_adj_diff <- calc_rast_cref(q_scd_adj_smry, info_smry,
                                   by = by,
                                   type_from = '_biomass',
                                   type_to = '_Q-rdiff-cref')

sei_scd_adj_diff <- calc_rast_cref(
  sei_scd_adj_smry, 
  info = info_sei_smry,
  by = by,
  type_from = '_SEI_SEI',
  type_to = '_SEI_SEI-rdiff-cref')


# Output -------------------------------------------------------------

# * prep files --------------------------------------------------------------
names(cov_scd_adj_smry) <- str_replace(names(cov_scd_adj_smry), '_biomass', 
                                       '_cover')

names(q_scd_adj_smry) <- str_replace(names(q_scd_adj_smry), '_biomass', 
                                     '_Q')

q_sei_smry <- c(q_scd_adj_smry, sei_scd_adj_smry)
q_sei_diff <- c(q_scd_adj_diff, sei_scd_adj_diff)


# * write files -----------------------------------------------------------

# summaries across GCMs
if(!test_run){
writeRaster(cov_scd_adj_smry, 
            file.path("data_processed/interpolated_rasters/", v_interp,
                      paste0(runv, "_cover_scd-adj_summary.tif")),
            overwrite = TRUE)

writeRaster(q_sei_smry, 
            file.path("data_processed/interpolated_rasters/", v_interp,
                      paste0(runv, "_q-sei_scd-adj_summary.tif")),
            overwrite = TRUE)

# differences relative to current climate (comparison within grazing level)
writeRaster(cov_scd_adj_diff, 
            file.path("data_processed/interpolated_rasters/", v_interp,
                      paste0(runv, "_cover-rdiff-cref_scd-adj_summary.tif")),
            overwrite = TRUE)

writeRaster(q_sei_diff, 
            file.path("data_processed/interpolated_rasters/", v_interp,
                      paste0(runv, "_q-sei-rdiff-cref_scd-adj_summary.tif")),
            overwrite = TRUE)

 write_csv(sei_mean_core1, paste0('data_processed/raster_means/', runv, 
                                  '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv'))
 
 write_csv(c3eco_sei2_smry, paste0('data_processed/raster_means/', runv, 
                                   '_sei-mean_scd-adj_smry-by-region-c3.csv'))
 
 write_csv(c3eco_sei2_gcm, paste0('data_processed/raster_means/', runv, 
                                   '_sei-mean_scd-adj_by-GCM-region-c3.csv'))

}


# remove temporary files --------------------------------------------------

if(remove_tmp_files) {
  paths <- list.files('tmp/', pattern = '.tif', full.names = TRUE)
  file.remove(paths)
}

