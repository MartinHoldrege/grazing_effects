# Martin Holdrege

# Script started 1/6/2022

# Purpose: to use the rMultivariateMatching package to do nearest neighbor
# matching of the 200 simulated sites onto gridded weather data in the west,
# and use those matches to interpolate the STEPWAT2 biomass
# output across the west.

# Note: all the code below strongly borrows from code provided by Rachel Renne
# and uses her R package. 


# params ------------------------------------------------------------------

rerun <- FALSE # re-create rasters that have already been interpolated?
test_run <- FALSE
date <- "20250228" # for appending to select file names
version <- 'v4'
# v1--1st version of criteria (i.e. based on calculating mathcing criteria across the region)
# v2--second version of criteria (i.e. basing matching criteria just on 200 sites)
# v3--criteria used in palmquist et al 2021 and renne et al 2024
# v4-- criteria same as v3, but new study area mask used (more restrictive)
run_climate <- FALSE # whether to upscale the climate data (doesn't need to be
run_climate_daymet <- FALSE # create a climate interpolation, not interpolating
run_fire <- FALSE
# interpolate data showing which vars caused changes in fire?
run_fire_attribution <- TRUE 
run_util <- FALSE # utilization

# for filtering output, put NULL if don't want to filter that variable
PFT2run =c('Sagebrush', 'C3Pgrass', 'C4Pgrass', 'Cheatgrass', 'Pforb', 'Shrub', 
           'Pherb', 'Aherb', 'Aforb','Pgrass')

PFT2run_util <- c('Sagebrush', 'Pherb', 'Aherb') # for utilization

run2run = 'fire1_eind1_c4grass1_co20_2503' # NULL # 
years2run =  NULL #'Current'

v2paste <- if(version == 'v1') "" else version # for pasting to interpolated tiffs

subset_in_target <- FALSE
exclude_poor_matches <- TRUE

# dependencies ------------------------------------------------------------

library(tidyverse)
library(raster)
library(rMultivariateMatching)
library(doParallel)
library(dtplyr)
source("src/general_functions.R")
source("src/mapping_functions.R")
source("src/interpolation_functions.R")

# read in data ------------------------------------------------------------
# list dataframes of stepwat2 output created in "scripts/02_summarize_bio.R" 
bio <- readRDS('data_processed/site_means/summarize_bio.RDS') 

# template raster, only 1km pixels that contain at least 50% scd pixels
# file created in 01_interpolation_data.R

# actual area for interpolation
path_template <- "data_processed/interpolation_data/cellnumbers_scd50.tif"
template <- raster(path_template)

# cell numbers corresponding to the climate data (broader b/ some
# sites outside of the new mask, but cellnumbers still match that in the template)
template_clim <- raster("data_processed/interpolation_data/cellnumbers.tif")

if(version < 'v4') {
  template <- template_clim
}

# Read in raw bioclim data for sagebrush extent and set cellnumbers as rownames
# (tc stands for 'targetcells')
tc0 <- read_csv("data_processed/interpolation_data/clim_for_interp_daymet-v4_1991-2020.csv",
                col_types = cols(site_id = 'd'),
                show_col_types = FALSE)

if(run_fire_attribution) {
  # created in 02b_fire_attribution.R
  attrib <- readRDS("data_processed/site_means/fire_dominant_drivers.RDS")
}

# Prep data -------------------------------------------------------------
# prep input data files for interpolation

# * Target cell data ------------------------------------------------------
# prepare the bioclim data

tc0 <- as.data.frame(tc0) # colnmaes for tibbles are deprecated

# need x and y coordinates for the rmultivarmatch
coords <- as.data.frame(raster::coordinates(template_clim))
coords$cellnumber <- values(template_clim) 
tc1 <- tc0 %>% 
  left_join(coords, by = 'cellnumber')

rownames(tc1) <- tc1$cellnumber # needed for multivarmatch fun

# limit to just the 6 variables of interest (keep lat/long)
# these are the same bioclim variables R.R. used
bioclim_vars <- c("bio1", "bio4", "bio9", "bio12", 
                  "bio15", "bio18")

tc2a <- tc1[, c("cellnumber", 'site_id', 'x', 'y', bioclim_vars)]


# criteria for matchingvars function (here using 10% of range of data)
# this is for scaling the variables
if (version == 'v1') {
  criteria <- map_dbl(tc2[, bioclim_vars], function(x) {
    (max(x) - min(x))*0.1
  })
} else if(version == 'v2') {
  criteria_v2 <- map_dbl(tc2[!is.na(tc2$site_id), bioclim_vars], function(x) {
    (max(x) - min(x))*0.1
  })
  criteria <- criteria_v2
  # for v2, switching to base matching variables on the grid-cells with sites
} else if (version == 'v3' | version == 'v4') {
  # criteria used in Renne et al 2024 
  criteria_v3 <- c(bio1 = 1.55, 
                   bio4 = 61.1, 
                   bio9 = 3.53, 
                   bio12 = 84.9, 
                   bio15 = 9.21, 
                   bio18 = 33.8
  )
  criteria <- criteria_v3
}
criteria

write_csv(data.frame(variable = names(criteria), criteria = criteria), 
          paste0("data_processed/interpolation_data/criteria-for-interp_", version, '.csv'))

# * subset cell data ------------------------------------------------------
# location and climate data for the 200 sites where simulations were actually
# done

# sc stands for 'subset cell'
sc1 <-  tc2a %>% 
  filter(!is.na(site_id))


# filter down the tc2a dataframe to only include the 
# cells in the template
cellnumbers <- as.numeric(values(template))
cellnumbers <- cellnumbers[!is.na(cellnumbers)]
tc2 <- tc2a %>% 
  filter(cellnumber %in% cellnumbers)

# check (should be 200 sites)
if(nrow(sc1) != 200 | any(duplicated(sc1$site_id))) {
  stop("Problem with dataframe")
}

# stepwat2 outputs --------------------------------------------------------

# files created in 02_summarize_bio.R script
# this is data for the 5 plant functional types (PFTs) of interest

# * biomass (g/m2) ----------------------------------------------------------
# actual biomass (ie not difference)
# converting to wide format (w stands for 'wide'), 
# for using with interpolationPoints function below

# this is every PFT, for all GCMs--i.e. all derived biomass change 
# variables could be calculated from these values after up-scaling,
# therefore %change of biomass no longer up-scaled

pft5_bio_w1 <- bio$pft5_bio1 %>% 
  filter_scenarios(PFT = PFT2run, run = run2run, years = years2run) %>% 
  ungroup() %>% 
  mutate(
    run2 = paste0(run, v2paste),
    id = paste(run2, PFT, "biomass", id, GCM, sep = "_")) %>% 
  dplyr::select(site, id, biomass) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "biomass")

# joining in cell numbers
pft5_bio_w2 <- join_subsetcells(step_dat = pft5_bio_w1, sc_dat = sc1,
                                subset_in_target = subset_in_target)


# * utilization -----------------------------------------------------------
# just interpolating the across GCM summaries
if(run_util) {
  tmp <- bio$pft5_bio2 %>% 
    filter(PFT %in% PFT2run_util) %>% 
    dplyr::select(run, years, RCP, graze, id, site, PFT, matches('utilization')) %>% 
    pivot_longer(cols = matches('utilization')) %>% 
    mutate(summary = str_extract(name, "(?<=_).+"),
           summary = ifelse(is.na(summary), "median", summary),
           run2 = paste0(run, v2paste),
           id = paste(run2, PFT, "util", id, summary, sep = "_"))
  
  s <- c('median', 'low', 'high')
  stopifnot(tmp$summary %in% s,
            s %in% tmp$summary)
  
  util_w1 <- tmp %>% 
    dplyr::select(site, id, value) %>% 
    pivot_wider(id_cols = "site",
                names_from = "id",
                values_from = "value") %>% 
    join_subsetcells(step_dat = ., sc_dat = sc1,
                     subset_in_target = subset_in_target)
    
}

# * wildfire --------------------------------------------------------------

if(run_fire) {
fire_w1 <- bio$fire0 %>% 
  filter_scenarios(run = run2run, years = years2run) %>% 
  ungroup() %>% 
  # don't want to bother upscale fire prob on runs where fire not simulated
  filter(!str_detect(run, 'fire0')) %>% 
  mutate(id = paste(run, "fire-prob", id, GCM, sep = "_")) %>% 
  dplyr::select(site, id, fire_prob) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "fire_prob") %>% 
  join_subsetcells(sc_dat = sc1, subset_in_target = subset_in_target)

}


# * fire attribution ------------------------------------------------------

if (run_fire_attribution) {
  # climate effect
  one_change_w1 <- attrib$one_change_smry1 %>% 
    filter_scenarios(run = run2run, years = years2run) %>% 
    mutate(id = paste0(run, v2paste, "_","fire-delta-", 
                       pred_var_cur, "_", id, "_median")) %>% 
    dplyr::select(site, id, delta_1var) %>% 
    pivot_wider(id_cols = "site",
                names_from = "id",
                values_from = "delta_1var")%>% 
    join_subsetcells(sc_dat = sc1, subset_in_target = subset_in_target)
  
  dominant_driver_w1 <- attrib$dominant_driver1 %>% 
    filter_scenarios(run = run2run, years = years2run) %>% 
    mutate(id = paste0(run, v2paste, "_","fire-dom-driver", "_", id, "_mode"),
           dominant_driver = as.numeric(dominant_driver)) %>% 
    dplyr::select(site, id, dominant_driver) %>% 
    pivot_wider(id_cols = "site",
                names_from = "id",
                values_from = "dominant_driver") %>% 
    join_subsetcells(sc_dat = sc1, subset_in_target = subset_in_target)
  
  # combine
  fire_attrib_w <- dominant_driver_w1 %>% 
    left_join(one_change_w1, by = 'cellnumber')
  
  # grazing effect
  one_change_gref_w1 <- attrib$one_change_gref_smry1 %>% 
    filter_scenarios(run = run2run, years = years2run) %>% 
    mutate(id = paste0(run, v2paste, "_","fire-delta-gref-", 
                       pred_var_ref, "_", id, "_median")) %>% 
    dplyr::select(site, id, delta_1var) %>% 
    pivot_wider(id_cols = "site",
                names_from = "id",
                values_from = "delta_1var")%>% 
    join_subsetcells(sc_dat = sc1, subset_in_target = subset_in_target)
  
  dominant_driver_gref_w1 <- attrib$dominant_driver_gref1 %>% 
    filter_scenarios(run = run2run, years = years2run) %>% 
    mutate(id = paste0(run, v2paste, "_","fire-dom-driver-gref", "_", id, "_mode"),
           dominant_driver = as.numeric(dominant_driver)) %>% 
    dplyr::select(site, id, dominant_driver) %>% 
    pivot_wider(id_cols = "site",
                names_from = "id",
                values_from = "dominant_driver") %>% 
    join_subsetcells(sc_dat = sc1, subset_in_target = subset_in_target)
  
  # combine
  fire_attrib_gref_w <- dominant_driver_gref_w1 %>% 
    left_join(one_change_gref_w1, by = 'cellnumber')
  
}

# * climate ---------------------------------------------------------------
# climate data doesn't have a 'run' name attached b/ runs have same climate
# but date attached in case in future runs are based on a different climate
# dataset
clim_all_w1 <- bio$clim_all2 %>% 
  filter(.data$years %in% years2run | is.null(years2run)) %>% 
  pivot_longer(cols = c("MAP", "MAT", 'psp')) %>% 
  ungroup() %>% 
  mutate(id = paste(name, "climate", RCP, years, GCM, date, version, sep = "_")) %>% 
  dplyr::select(site, id, value) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "value")

# joining in cell numbers
clim_all_w2 <- join_subsetcells(step_dat = clim_all_w1, sc_dat = sc1,
                                subset_in_target = subset_in_target)

# code removed for following sections--see 2022 commits to get code
# *crossing threshold ----------------------------------------------------

# * wgcm bio diff ---------------------------------------------------------

# * bio diff --------------------------------------------------------------

# identify matches --------------------------------------------------------
# ID matches from subset cells for all Target cells (i.e. calculates distance)

# check that subset_in_target is true, for fun call below
if(!all(rownames(sc1 %in% rownames(tc2)))) {
  warning("Not all subset cells in target cells")
}

# note--consider saving the figure of matching quality
# this only takes ~ min to run
# not sure why "In max(raster::values(x), na.rm = TRUE..."
# warning is being thrown, but output looks ok
#no non-missing arguments to max; returning -In

# STOP: continue here--troubleshoot the
# the problem that subset_in_target = FALSE flag
# arg should be used, but that's causing additional problems
sc2 <- sc1[, c("site_id", 'x', 'y', bioclim_vars)]
rownames(sc2) <- sc2$site_id

# creating a hash so multivarmatch doesn't need to be re-run unecessarily
# these are the arguments passed to multivarmatch below
args <- list(
  matchingvars = drop_na(dplyr::select(tc2, -site_id)),
  # subset cells include cells that are not
  subsetcells = sc2,
  matchingvars_id = "cellnumber",
  subsetcells_id = "site_id",
  criteria = criteria,
  raster_template = template,
  subset_in_target = subset_in_target,
  addpoints = FALSE
)
hash <- digest::digest(args)

path_hash <- "data_processed/interpolation_data/multivarmatch_hash.txt"
path_match <- "data_processed/interpolation_data/multivarmatch.csv"
rerun_match <- TRUE

if(!file.exists(path_hash)) {
  write_lines(hash, path_hash)
} else {
  old_hash <- read_lines(path_hash)
  
  if(old_hash == hash & file.exists(path_match)) {
    rerun_match <- FALSE
    match1 <- read_csv(path_match, show_col_types = FALSE) %>% 
      as.data.frame()
  } else {
    write_lines(hash, path_hash)
  }
}

if(rerun_match) {
  match1 <- rlang::exec(multivarmatch, !!!args)
  write_csv(match1, path_match)
}

# *plotting interpolation quality -----------------------------------------

tmp_site_id <- match1 %>% 
  mutate(
    # site from which data pulled
    site_id = as.numeric(subset_cell),
    # cellnumber of the cell interpolated to
    cellnumber = as.numeric(target_cell)) %>% 
  dplyr::select(site_id, cellnumber)

# for each pixels this raster shows which stepwat site the data was drawn
# from for interpolation
interp_ordered <- fill_raster(tmp_site_id, template = template,
                              use_template_values = TRUE)

writeRaster(interp_ordered, 
            paste0("data_processed/interpolation_data/interp_locations_200sites_", version, ".tif"),
            overwrite = TRUE)

set.seed(1234)

rand <- sample(1:200, 200)
pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, 'RdBu')))(200)
jpeg(paste0("figures/interpolation_quality/interpolation_areas_", version, ".jpeg"),
     res = 600, units = 'in', width = 8, height = 8)
  plot(interp_ordered, main = 'interpolation areas (random colors)',
       col = pal[rand])
dev.off()


# proportion cells with decent matching quality
mean(match1$matching_quality < 1.5) # ~87%

match_quality <-  match1 %>% 
  mutate(cell_num = as.numeric(target_cell)) %>% 
  dplyr::select(cell_num, matching_quality) %>% 
  fill_raster(template, use_template_values = TRUE)

jpeg(paste0("figures/interpolation_quality/matching_quality_", version, ".jpeg"),
     res = 600, units = 'in', width = 8, height = 8)
plot(match_quality, main = 'matching quality', breaks = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 10),
     col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))[3:10])
dev.off()

writeRaster(match_quality, 
            paste0("data_processed/interpolation_quality/matching_quality", 
                   version, '.tif'),
            overwrite = TRUE)
# interpolation ---------------------------------------------------------
# use matches from multivarmatch to interpolate the STEPWAT2 output
# across all the grid cells

# code removed here--see 2022 commits (for interpolating non-biomass dataframes)

# * biomass ---------------------------------------------------------------

# on my computer to run all PFTs/scenarios/GCMs would take 21.5 hours!, 
# when run in parallel should take ~6.5 hours, 3.26 hrs if RCP4.5 is not included. 
# It actually took 3 hrs, 40 min, as it is now. 

# I'm using ideas from here:
# https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
# to run this in parallel

# num cores, this includes logical cores (threads)
num.cores <- parallel::detectCores() 
registerDoParallel(num.cores)

# if file was corrupted, and really small then re-create
min_size <- file.size(path_template)*0.2 # 0.5 multiplier is arbitrary, and maybe not strict enough

# which columns haven't already been upscaled
path_bio <- file.path("data_processed/interpolated_rasters/biomass", version)

# separately running interpolation on different sets of columns
run_interpolatePoints(df = pft5_bio_w2,
                      match = match1,
                      template = template,
                      path = path_bio,
                      min_size = min_size,
                      rerun = rerun,
                      exclude_poor_matches = exclude_poor_matches)

# * utilization ------------------------------------------------------------------
if(run_util) {
  
  path_util <- file.path("data_processed/interpolated_rasters/util", version)
  
  # separately running interpolation on different sets of columns
  run_interpolatePoints(df = util_w1,
                        match = match1,
                        template = template,
                        path = path_util,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)
  
}


# * fire ------------------------------------------------------------------
if(run_fire) {

  path_fire <- file.path("data_processed/interpolated_rasters/fire", version)

  # separately running interpolation on different sets of columns
  run_interpolatePoints(df = fire_w1,
                        match = match1,
                        template = template,
                        path = path_fire,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)

}

# * fire attribution ------------------------------------------------------

if(run_fire_attribution) {

  path_attrib <- file.path("data_processed/interpolated_rasters/fire_attrib", version)

  run_interpolatePoints(df = fire_attrib_w,
                        match = match1,
                        template = template,
                        path = path_attrib,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)
  
  run_interpolatePoints(df = fire_attrib_gref_w,
                        match = match1,
                        template = template,
                        path = path_attrib,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)

}


# * climate ---------------------------------------------------------------

path_clim <- file.path("data_processed/interpolated_rasters/climate", version)

if (run_climate) {

  run_interpolatePoints(df = clim_all_w2,
                        match = match1,
                        template = template,
                        path = path_clim,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)
  
}

if(run_climate_daymet){
  daymet1 <- tc1 %>% 
    filter(!is.na(site_id)) %>% 
    dplyr::select(cellnumber, ptcor, bio1, bio12) %>% 
    rename(PTcor = ptcor, MAT = bio1, MAP = bio12)

  names(daymet1)[-1] <- paste0(names(daymet1)[-1], "_daymet-climate_", date, "_", version)
  
  run_interpolatePoints(df = daymet1,
                        match = match1,
                        template = template,
                        path = path_clim,
                        min_size = min_size,
                        rerun = rerun,
                        exclude_poor_matches = exclude_poor_matches)
}

# When you're done, clean up the cluster
stopImplicitCluster()
