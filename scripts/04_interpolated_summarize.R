# Martin Holdrege

# Script started Jan 7, 2022

# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script, which can takes ~ 4 hours to run on
# my laptop (when excluding RCP4.5)

# this script takes ~ 10 min to run


# dependencies ------------------------------------------------------------

library(tidyverse) 
library(terra)
source("src/mapping_functions.R")
source("src/general_functions.R")

# params ------------------------------------------------------------------

v <- 'v4' # interpolation version
# multiple runs may need to listed if different runs done for different
# grazing levels

runs <- c('fire1_eind1_c4grass1_co20_2502')
run <- runs[[1]] # used for output file naming

# read in data ------------------------------------------------------------

# * rasters ---------------------------------------------------------------

# data up-scaled for each GCM
regex <- if(length(runs) == 1) runs else paste0('(', paste(runs, collapse = ')|('), ")")
bio_files <- list.files(file.path("data_processed/interpolated_rasters/biomass", v),
                        pattern = regex,
                        full.names = TRUE)
length(bio_files)

# for testing use a subset for computational efficiency
# bio_files <- bio_files[sample(1:length(bio_files), 100)]

# this is kind of a 'raster stack', where each layers is a tif
# it's not actually loaded into memory
rast1 <- terra::rast(bio_files) # class SpatRast

# params ------------------------------------------------------------------

# some terra operations can be done in parallel and need 
# to know num of cores
num.cores <- parallel::detectCores(logical = FALSE) 

# names of layers --------------------------------------------------
# descriptions of layers of raster stack
r_names <- names(rast1)

rast_info <- create_rast_info(rast1, id_noGCM = TRUE) %>% 
  mutate(layer_num = 1:nrow(.)) %>% 
  # the ordering is important for later creation of spatraster dataset
  arrange(id_noGCM, id)


# median by GCM -----------------------------------------------------------
# for each treatment/scenario combination calculate the median
# across GCMs

# split raster by scenario

rast2 <- rast1[[rast_info$id]] # make sure rows are ordered correctly

# list of rasters, each element is a SpatRaster including all GCMs
# for a given treatment/scenario combination
rast_gcm_l <- split(rast2, f = rast_info$id_noGCM)

names(rast_gcm_l) <- unique(rast_info$id_noGCM)
length(rast_gcm_l)

# SpatRasterDataset (i.e. each dataset includes all GCMs for treatment/scenario)
# sds_gcm <- sds(rast_gcm_l) 
# class(sds_gcm)
# 
# check that the number of layers in each dataset is correct

stopifnot(map_dbl(rast_gcm_l, nlyr) == rast_info %>%
            group_by(id_noGCM) %>%
            summarise(n = n()) %>%
            pull(n)
)

# note these median rasters now inlcude the median under current
# conditions (i.e. medians of one value). 
med2 <- rast(map(rast_gcm_l, app, fun = "median"))

# delta biomass cref ------------------------------------------------------
# change in biomass (absolulte, not relative) relative to ambient climate (c)
# conditions, calculated within a grazing level
# naming: bio-rdiff-cref, biomass raw difference climate refence (i.e. difference
# from historical conditions)

rast_info_med <- create_rast_info(med2, 
                                 into = c("PFT", "type", "RCP", "years", 
                                          "graze")) %>% 
  mutate(id3 = paste(run2, PFT, sep = "_"))

info_med_l <- split(rast_info_med, f = rast_info_med$id3)

diff_cref1 <- map(info_med_l, function(df) {
  id_current <- df$id[df$RCP == 'Current']
  stopifnot(length(id) ==1)
  ids_future <- df$id[df$RCP != 'Current']
  
  delta <- med2[[ids_future]] - med2[[id_current]]
  new_names <- names(delta) %>% 
    str_replace("biomass", "bio-rdiff-cref")
  names(delta) <- new_names
  delta
})

diff_cref2 <- rast(diff_cref1)
names(diff_cref2) <- map(diff_cref1, names) %>% 
  unlist()


# delta biomass gref  -----------------------------------------------------------
# delta biomass from current heavy to future light or moderate grazing
# specifically, future is RCP 8.5. Because the reference level
# is current, I'm calculating this using median rasters (which in this case
# should yield the same answer as calculate for each gcm then taking the median).

gref_info <- rast_info %>% 
  filter_clim_extremes() %>% 
  select(-GCM, -id, -layer_num) %>% 
  distinct() %>% 
  filter((RCP == 'Current' & graze == "Heavy")|
           (RCP == 'RCP45' & graze %in% c("Light", "Moderate")))

gref_info <- split(gref_info, 
                   paste(gref_info$RCP, gref_info$graze, sep = "_")) %>% 
  map(function(x) arrange(x, PFT))

stopifnot(
  gref_info$Current_Heavy$PFT == gref_info$RCP45_Moderat$PFT,
  gref_info$Current_Heavy$PFT ==  gref_info$RCP45_Light$PFT,
  # this only works if there is one raster per PFT (i.e. reference raster)
  unique(gref_info$Current_Heavy$PFT) == gref_info$Current_Heavy$PFT
  )

# current heavy and future light or moderate grazing
current_heavy <- subset(med2, subset = gref_info$Current_Heavy$id_noGCM)
future_light <- subset(med2, subset = gref_info$RCP45_Light$id_noGCM)
future_moderate <- subset(med2, subset = gref_info$RCP45_Moderate$id_noGCM)

# max under current conditions
max_heavy <- t(minmax(current_heavy))[, 2]

# % change to future light  relative to current heavy grazing
diff_gref2light <- (future_light - current_heavy)/max_heavy*100

# % change  to future moderate relative to current heavy
diff_gref2moderate <- (future_moderate - current_heavy)/max_heavy*100

stopifnot(names(max_heavy) == gref_info$Current_Heavy$id_noGCM)

# gref stands for fixed grazing reference
rast_diff_gref <- rast(list(diff_gref2light, diff_gref2moderate))

names(rast_diff_gref) <- names(rast_diff_gref) %>% 
  str_replace("_biomass_", '_bio-diff-gref-cur-heavy_')

#  c3Pgrass/Pgrass -------------------------------------------------------

# grass_info <- rast_info %>% 
#   filter(PFT %in% c("C3Pgrass", "Pgrass"),
#          # just looking at end of century and current for now
#          years %in% c("2070-2100", "Current")) %>% 
#   arrange(PFT, c4, graze, RCP, years) %>% # important so comparing correct lyrs
#   # creates two identical df's except the PFT and id columns are different
#   split(f = .$PFT, drop = TRUE)
# 
# 
# grass_id <- map(grass_info, function(x) x$id) # id values
# grass_id_noGCM <- map(grass_info, function(x) x$id_noGCM) # id values
# 
# # seperate rasters for Pgrass and C3Pgrass
# grass_r1 <- map(grass_id, function(x) rast1[[x]])
# 
# # grass ratio 
# gratio_r1 <- grass_r1$C3Pgrass/grass_r1$Pgrass
# 
# # check that number of layers has been preserved
# stopifnot(nlyr(gratio_r1) == length(grass_id[[1]]))

# ** median ratio ----------------------------------------------------------

# # info on the median lyrs
# grass_info_med <- grass_info$C3Pgrass %>% 
#   arrange(graze, RCP, years, c4) %>% 
#   select(-id, -GCM, -layer_num) %>% 
#   filter(!duplicated(.))
# 
# # split into list
# gratio_r_l <- split(gratio_r1, grass_id_noGCM$C3Pgrass)
# 
# # naming so layers will be identifiable, C3Pgrass is in name
# # but this actuall C3Pgrass/Pgrass ratio
# names(gratio_r_l) <- unique(grass_id_noGCM$C3Pgrass)
# 
# # take median across GCMs
# gratio_med1 <- map(gratio_r_l, app, fun = "median")
# 
# # this is operation takes a few minutes
# gratio_med2 <- rast(gratio_med1) # put back into one SpatRaster
# nlyr(gratio_med2)
# names(gratio_med2)

# ** Pgrass ----------------------------------------------------------------

# Pgrass_id_noGCM <- names(med2) %>% 
#   str_subset("_Pgrass_") %>% #on Pgrass
#   str_subset("_2030-2060_", negate = TRUE) # not mid century
# 
# # information on Pgrass, medians across GCMs (same as info for C3Pgrass)
# Pgrass_info_med<- grass_info_med %>% 
#   mutate(id_noGCM = str_replace(id_noGCM, "_C3Pgrass","_Pgrass")) %>% 
#   select(-PFT)
# 
# Pgrass_target_lyrs <- Pgrass_info_med %>% 
#   filter(RCP != "Current") %>% 
#   pull(id_noGCM)
# 
# Pgrass_ref_lyrs <- create_ref_id(Pgrass_target_lyrs)
# 
# # scaled % change relative to current grazing of the same intensity
# rast_d_Pgrass <- rast_diff(rast = med2,
#                            target_layer = Pgrass_target_lyrs,
#                            ref_layer = create_ref_id(Pgrass_target_lyrs))
# 
# names(rast_d_Pgrass)


# save files ---------------------------------------------------------

# * info files ------------------------------------------------------------

# saveRDS(rast_info, 
#         "data_processed/interpolated_rasters/raster_info.RDS")
# 
# saveRDS(grass_info_med,
#         "data_processed/interpolated_rasters/grass_info_med.RDS")
# 
# # info on median pgrass
# saveRDS(Pgrass_info_med,
#         "data_processed/interpolated_rasters/Pgrass_info_med.RDS")

# * rasters ---------------------------------------------------------------

# median across GCMs, for all future scenarios

writeRaster(med2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0(run, "_bio_future_median_across_GCMs.tif")),
            overwrite = TRUE)

# difference in biomass (raw) relative to current conditions (within a grazing level)
writeRaster(diff_cref2, 
            file.path("data_processed/interpolated_rasters", 
                      paste0(run, "_bio-rdiff-cref_median.tif")),
            overwrite = TRUE)


# c3/total grass ratio
# writeRaster(gratio_med2, "data_processed/interpolated_rasters/C3Pgrass-Pgrass-ratio_by-scenario_median.tif",
#             overwrite = TRUE)

# scaled percent change from current to future (wgraze), for Pgrass

# writeRaster(rast_d_Pgrass, "data_processed/interpolated_rasters/Pgrass_bio-diff-wgraze_median.tif",
#             overwrite = TRUE)

# # scaled % change from current light grazing to future (RCP8.5-mid) heavy graze
writeRaster(rast_diff_gref,
            file.path("data_processed/interpolated_rasters",
            paste(run, "_bio-diff-gref-cur-heavy_median.tif")),
            overwrite = TRUE)
