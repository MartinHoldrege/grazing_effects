# Purpose: 
# determine what SEI component (Q score) is most responsible for
# the climate or grazing driven pixel-wise change in SEI (both proportional
# change in Q for RGB maps, and which variable is the dominant driver of change)
# - Primary driver per GCM = PFT with largest |proportional change|

# --- params / setup -----------------------------------------------------------
source("src/params.R")                     # defines opt (run, v_interp, years, etc.)
# this script operates on the entire rasters so 'regions' don't apply here
v_interp <- opt$v_interp
run      <- opt$run
years    <- opt$years
test_run <- opt$test_run
yr_lab   <- opt$yr_lab                     # empty by default (see params.R)
runv     <- paste0(run, v_interp)

dir_tmp  <- file.path('tmp', stringr::str_replace(yr_lab, "_", ""))
dir_out  <- file.path("data_processed", "interpolated_rasters", "sei_attrib", 
                      v_interp)
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

graze_ref <- 'Moderate' # reference grazing level

# same minimum sei change cutoff used in Holdrege et al 2024 REM
diff_cutoff <- 0.01 # SEI change less than this not considered for attribution

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source("src/general_functions.R")
source("src/mapping_functions.R")
source('src/SEI_functions.R')
source('src/SEI_attrib_functions.R')

# read in data ------------------------------------------------------------

# Input produced by 05_interpolated_summarize_sei_scd-adj.R
path_q <- file.path(dir_tmp, "q_scd_adj1.tif")
stopifnot(file.exists(path_q))

r_q <- rast(path_q)                        # layers: Q for Sagebrush/Pherb/Aherb, by GCM & RCP
info <- create_rast_info(r_q, id_noGCM = TRUE) %>%
  arrange(PFT, graze, RCP, years, GCM, id) %>% 
  select(-run2, -type)


path_sei <- file.path(dir_tmp, 'sei_scd_adj1.tif')
r_sei1 <- rast(path_sei)

sei_diff_clim_smry1 <- rast(file.path("data_processed/interpolated_rasters", v_interp,
                       paste0(runv, yr_lab, "_q-sei-rdiff-cref_scd-adj_summary.tif")))

# optional downsample for test runs
if (isTRUE(test_run)) {
  r_q <- downsample(r_q)
  r_sei1 <- downsample(r_sei1)
  sei_diff_clim_smry1 <- downsample(sei_diff_clim_smry1)
}


# driver of climate sei change --------------------------------------------

# * --- filter to current vs selected future period ------------------------------

# Q
info_cur <- info %>%
  filter(RCP == "Current") %>%
  select(-RCP, -years, -GCM, -id_noGCM)

info_fut <- info %>%
  filter(RCP != "Current") 

info_pair <- left_join(info_fut, info_cur,
                        by = c('run', 'PFT', 'graze'),
                        suffix = c('_fut', '_cur'))

# SEI
info_sei <- create_rast_info(r_sei1, into = c('variable', 'type', 'RCP', 'years', 'graze',
                                              'GCM'),
                             id_noGCM = TRUE) %>% 
  select(-variable, -type, -run2)

info_pair_sei_clim <- info_sei %>% 
  filter(RCP == 'Current') %>% 
  select(-RCP, -years, -GCM, -id_noGCM) %>% 
  right_join(filter(info_sei, RCP != 'Current'),
             by = join_by(run, graze),
             suffix = c('_cur', '_fut'))


# * reference vs other grazing levels -------------------------------------

# Q
info_pair_gref <- info %>%
  filter(graze == graze_ref) %>%
  select(-graze, -id_noGCM) %>% 
  right_join(filter(info, graze != graze_ref),suffix = c('_ref', ''),
            by = join_by(run, PFT, RCP, years, GCM)) 

# SEI
info_pair_sei_gref <- info_sei %>% 
  filter(graze == graze_ref) %>% 
  select(-graze, -id_noGCM) %>% 
  right_join(filter(info_sei, graze != graze_ref),suffix = c('_ref', ''),
             by = join_by(run,  RCP, years, GCM))


# * sei change ------------------------------------------------------------
# only interested in 'q' changes that are in the same direction as the sei
# change, and interested pix

sei_delta_clim1 <- r_sei1[[info_pair_sei_clim$id_fut]] - r_sei1[[info_pair_sei_clim$id_cur]]

sei_delta_gref1 <- r_sei1[[info_pair_sei_gref$id]] - r_sei1[[info_pair_sei_gref$id_ref]]


# ** median delta gref  ---------------------------------------------------

# not used here but created for other scripts
sei_delta_gref_smry1 <- summarize_gcms_raster(sei_delta_gref1,
                                              info = info_pair_sei_gref,
                                              include_low_high = FALSE)
# * proportional Q |change| ---------------------------------------------------


# ** climate reference ----------------------------------------------------

q_prop1 <- prop_change(fut = r_q[[info_pair$id_fut]],
                       cur = r_q[[info_pair$id_cur]])

info_pair2 <- info_pair_sei_clim %>% 
  select(-id_cur, -id_noGCM) %>% 
  rename(id_fut_sei = id_fut) %>% 
  right_join(info_pair, 
             by = join_by(run, graze, RCP, years, GCM))

# only considering proportional change in Q when 
# it is in the same direction as the change in SEI (i.e.
# it contributed to that directional sei change)
# then relativize (range 0-1)
q_prop_rel1 <- relativize_q_prop(
  qprop = q_prop1[[info_pair2$id_fut]],
  sei_diff = sei_delta_clim1[[info_pair2$id_fut_sei]],
  info = info_pair2,
  id = 'id_fut')

# calculating the mean instead of median as in Holdrege et al. b/ 
# then each sums to 1 (for RGB map making)
q_prop_rel_smry <- info_pair2 %>% 
  rename(id = id_fut) %>% 
  raster_gcm_smry(r = q_prop_rel1,fun = 'mean')

# sei change mask (within climate) --basing on median change
info_sei_diff_clim_smry <- create_rast_info(sei_diff_clim_smry1, 
                                            into = c('variable', 'type', 'RCP', 
                                                     'years', 'graze','summary'))

id <- info_sei_diff_clim_smry %>% 
  filter(variable == 'SEI', summary == 'median') %>% 
  pull(id)
sei_diff_clim_smry2 <- sei_diff_clim_smry1[[id]]

sei_diff_clim_mask <- ifel(abs(sei_diff_clim_smry2) >= diff_cutoff, 1, 0)

# ** grazing reference ----------------------------------------------------

q_prop_gref1 <- prop_change(fut = r_q[[info_pair_gref$id]],
                       cur = r_q[[info_pair_gref$id_ref]])

info_pair_gref2 <- info_pair_sei_gref %>% 
  select(-id_ref, -id_noGCM) %>% 
  rename(id_sei = id) %>% 
  right_join(info_pair_gref, 
             by = join_by(run, graze, RCP, years, GCM))


# relativized proportional change in q (corrected for direction of sei change)
q_prop_rel_gref1 <- relativize_q_prop(
  qprop = q_prop_gref1[[info_pair_gref2$id]],
  sei_diff = sei_delta_gref1[[info_pair_gref2$id_sei]],
  info = info_pair_gref2,
  id = 'id')

# calculating the mean instead of median as in Holdrege et al. b/ 
# then each sums to 1 (for RGB map making)
q_prop_rel_gref_smry <- raster_gcm_smry(info_pair_gref, r = q_prop_rel_gref1,
                                        fun = 'mean')

# median SEI change in response grazing
sei_diff_gref_smry <- raster_gcm_smry(info = info_pair_sei_gref, r = sei_delta_gref1, 
                        fun = 'median')


sei_diff_gref_mask <- ifel(abs(sei_diff_gref_smry) >= diff_cutoff, 1, 0)


# dominant driver of change -----------------------------------------------


# *climate effect ---------------------------------------------------------

info_fut_l <- prepare_info_dom(info_fut)  

# **mode across GCMs -------------------------------------------------------

# calculate dominant driver and then take the mode
q_dom3 <- dom_driver_mode(info_fut_l, q_prop_rel1)

# then mask locations with less than 0.01 change in SEI
tmp <- filter(info_sei_diff_clim_smry, summary == 'median',
       variable == 'SEI') %>% 
  select(RCP, years, graze, id) %>% 
  rename(id_sei = id)

info_fut_dom_smry <- prepare_info_dom_mode(info_fut_l) %>% 
  bind_rows() %>% 
  select(run, RCP, years, graze, id_noGCMPFT) %>% 
  distinct() %>% 
  left_join(tmp,
    # only works assuming run's are the same (run suffix's are different
    # so can't merge by them)
    by = c('RCP', 'years', 'graze'))

# '4' means no dominant driver, b/ negligible sei change
q_dom4 <- ifel(sei_diff_clim_mask == 1,
               q_dom3[[info_fut_dom_smry$id_noGCMPFT]], 4)

# helper to stamp the same labels onto every layer
pfts <- unique(as.character(info_fut_l[[1]]$PFT))
stopifnot(length(pfts) ==3)
pfts2 <- c(pfts, 'None')
set_pft_levels <- function(x, pfts) {
  lev <- data.frame(value = seq_along(pfts), PFT = c(pfts))
  for (i in seq_len(nlyr(x))) {
    levels(x)[[i]] <- lev
  }
  x
}

q_dom5 <- set_pft_levels(q_dom4, pfts2)
names(q_dom5) <- names(q_dom4)

# *grazing effect ---------------------------------------------------------

info_gref_l <- prepare_info_dom(info_pair_gref2)

# **mode across GCMs -------------------------------------------------------

q_dom_gref3 <- dom_driver_mode(info_gref_l, q_prop_rel_gref1)

# then mask locations with less than 0.01 change in SEI
tmp <- info_pair_sei_gref %>% 
  select(RCP, years, graze, id_noGCM) %>% 
  rename(id_sei = id_noGCM) %>% 
  distinct()

info_gref_dom_smry <- prepare_info_dom_mode(info_gref_l) %>% 
  bind_rows() %>% 
  select(run, RCP, years, graze, id_noGCMPFT) %>% 
  distinct() %>% 
  left_join(tmp,
            # only works assuming run's are the same (run suffix's are different
            # so can't merge by them)
            by = c('RCP', 'years', 'graze'))

# '4' means no dominant driver, b/ negligible sei change
q_dom_gref4 <- ifel(sei_diff_gref_mask[[info_gref_dom_smry$id_sei]] == 1,
               q_dom_gref3[[info_gref_dom_smry$id_noGCMPFT]], 4)

q_dom_gref5 <- set_pft_levels(q_dom_gref4, pfts2)
names(q_dom_gref5) <- names(q_dom_gref4)

# --- write outputs ------------------------------------------------------------
if (!test_run) {
# Categorical raster: primary driver mode across GCMs (values 1=Sagebrush, 2=Pherb, 3=Aherb, 4 = None)
# this is w/ respect to historical conditions under the same grazing level
# (i.e. within grazing level comparison)
writeRaster(
  names_replace(q_dom5, 'SEI_SEI-rdiff-cref', 'SEI_dom-driver'),
  # cref means 'climate reference'
  filename = file.path(dir_out, paste0(runv, '_', years, "_dom-driver-sei_cref-mode.tif")),
  overwrite = TRUE
)

#gref --grazing reference (i.e. across grazing trmt comparison,
  #relative to moderate grazing)
writeRaster(
  names_replace(q_dom_gref5, 'SEI_SEI', 'SEI_dom-driver'),
  filename = file.path(dir_out, paste0(runv, '_', years, "_dom-driver-sei_gref-mode.tif")),
  overwrite = TRUE
)

# (2) Continuous rasters: mean proportional change across GCMs for each PFT
#     (These three bands can be used as RGB input for visualization)
# climate effect
writeRaster(
  names_replace(q_prop_rel_smry, '_biomass_', '_qprop_'),
  filename = file.path(dir_out, paste0(runv, '_', years, "_qprop-cref_mean.tif")),
  overwrite = TRUE
)
# grazing effect
writeRaster(
  names_replace(q_prop_rel_gref_smry, '_biomass_', '_qprop_'),
  filename = file.path(dir_out, paste0(runv, '_', years, "_qprop-gref_mean.tif")),
  overwrite = TRUE
)

# change in SEI (grazing reference)
writeRaster(
  names_replace(sei_delta_gref_smry1, '_SEI_SEI_', '_SEI_SEI-rdiff-gref_'),
  filename = file.path(dir_out, paste0(runv, '_', years, "_delta-SEI-gref_smry.tif")),
  overwrite = TRUE
)
}
