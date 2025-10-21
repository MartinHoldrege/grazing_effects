# Purpose: Create summary figures and maps
# showing drivers of SEI change (RGB maps, as well as dominant driver maps)
# this is attribution of changes in SEI due to changes in climate
# and due to change in grazing

# Author: Martin Holdrege

# Date Started: Sept 16, 2025


# params ------------------------------------------------------------------

source("src/params.R")                     # defines opt (run, v_interp, years, etc.)
# this script operates on the entire rasters so 'regions' don't apply here
v_interp <- opt$v_interp
runv      <- opt$runv
years    <- opt$years
yr_lab <- opt$yr_lab
test_run <- opt$test_run

dir_dat  <- file.path("data_processed", "interpolated_rasters", "sei_attrib", 
                      v_interp)
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source('src/general_functions.R')
source('src/mapping_functions.R')
source('src/fig_params.R')
source('src/SEI_attrib_functions.R')

# read in data ------------------------------------------------------------


# *dominant driver of SEI change ------------------------------------------

# cref = climate reference (i.e. effect of climate w/ in a grazing level)\
# create in "scripts/06_sei_attrib.R"
d_cref1 <- rast(file.path(dir_dat, 
                            paste0(runv, '_', years, "_dom-driver-sei_cref-mode.tif")))

#gref --grazing reference (i.e. across grazing trmt comparison,
#relative to moderate grazing)
# create in "scripts/06_sei_attrib.R"
d_gref1 <- rast(file.path(dir_dat, 
               paste0(runv, '_', years, "_dom-driver-sei_gref-mode.tif")))


# change in SEI -----------------------------------------------------------

# cref
# produced by 05_interpolated_summarize_sei_scd-adj.R
sei_cref1 <- rast(
  file.path("data_processed/interpolated_rasters", v_interp,
            paste0(runv, yr_lab, "_q-sei-rdiff-cref_scd-adj_summary.tif")))

# gref
# create in "scripts/06_sei_attrib.R"
sei_gref1 <- rast(file.path(dir_dat, 
                            paste0(runv, '_', years, "_delta-SEI-gref_smry.tif")))

# * proportional Q change -------------------------------------------------
# NOT used at the moment
#  mean proportional change across GCMs for each PFT
#     (These three bands can be used as RGB input for visualization)
# climate effect
# values relativized to between 0 and 1

# qprop_cref1 <- rast(file.path(dir_dat,
#                               paste0(runv, '_', years, "_qprop-cref_mean.tif")))
# grazing effect
# qprop_gref1 <- rast(file.path(dir_dat,
#                               paste0(runv, '_', years, "_qprop-gref_mean.tif")))


# *downsample -------------------------------------------------------------


if(test_run) {
  d_cref1 <- downsample(d_cref1)
  d_gref1 <- downsample(d_gref1)
  sei_gref1 <- downsample(sei_gref1)
  sei_cref1 <- downsample(sei_cref1)
  # qprop_cref1 <- downsample(qprop_cref1)
  # qprop_gref1 <- downsample(qprop_gref1)
}
comb_gref1 <- c(sei_gref1, d_gref1)
comb_cref1 <- c(sei_cref1, d_cref1)

# apply mask --------------------------------------------------------------

# the dominant driver layers are masked, using that mask and applying
# it to the q proportional change layers
# qprop_cref2 <- mask_using_dom_layer(r = qprop_cref1, dom = d_cref1)
# qprop_gref2 <- mask_using_dom_layer(r = qprop_gref1, dom = d_gref1)

# prepare info ------------------------------------------------------------

into <- c("variable", "type", "RCP", "years", "graze")

info_gref1 <- bind_rows(
  create_rast_info(sei_gref1, into = c(into, 'summary')),
  create_rast_info(d_gref1, into = into)
) %>% 
  filter(is.na(summary) | summary == 'median')

rcps <- unique(info_gref1$RCP)
rcps <- rcps[rcps != 'Current']

info_cref1 <- create_rast_info(sei_cref1, into = c(into, 'summary')) %>% 
  filter(variable == 'SEI') %>% 
  bind_rows(create_rast_info(d_cref1, into = c(into, 'summary'))) %>% 
  filter(summary == 'median') 

# * gref: dSEI, driver 12 panel -------------------------------------------------

# 12 panel plot where top two rows show change in SEI, relative
# to moderate grazing (columns are grazing levels) for historical and future conditions (rows)
# next two rows are the same structure but showing the primary driver 
# of sei change

for(rcp in rcps) {
  info_gref2 <- info_gref1 %>% 
    filter(RCP %in% c('Current', rcp),
           years %in% c('Current', !!years)
    )
  
  p2 <- plot_dsei_drvr_gref(info = info_gref2,
                            r = comb_gref1)
  
  filename <- paste0('figures/sei/attrib/map_dSEI_drvr_gref_12panel_',
                     rcp, '_', years, '_', runv, '.png')
  
  ggsave(filename = filename,
         plot = p2,
         width = 7,
         height = 8,
         dpi = 800)
}


# * cref: dSEI, driver 8 panel -----------------------------------------------

# 8 panel plot, showing (left column) the change in sei
# and (right column) the driver of the change
# relative to historical climate ('climate reference 'cref')
# rows are grazing levels

info_cref2 <- info_cref1 %>% 
  filter(RCP == rcp,
         years == !!years
  )
info <- info_cref2
r <- comb_cref1

for(rcp in rcps) {
  info_cref2 <- info_cref1 %>% 
    filter(RCP == rcp,
           years == !!years
    )
  
  p2 <- plot_dsei_drvr_cref(info = info_cref2,
                            r = comb_cref1)
  
  filename <- paste0('figures/sei/attrib/map_dSEI_drvr_cref_8panel_',
                     rcp, '_', years, '_', runv, '.png')
  
  ggsave(filename = filename,
         plot = p2,
         width = 5.1,
         height = 8,
         dpi = 800)
}

