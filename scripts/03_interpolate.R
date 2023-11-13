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
test_run <- FALSE # TRUE # 
date <- "20230712" # for appending to select file names
run_climate <- FALSE # whether to upscale the climate data (doesn't need to be
# rerun unless climate variables are changed/updated)

# dependencies ------------------------------------------------------------

source("scripts/02_summarize_bio.R") # creates dataframes of stepwat2 output
library(raster)
library(rMultivariateMatching)
library(doParallel)

# read in data ------------------------------------------------------------

# template raster
path_template <- "data_processed/interpolation_data/cellnumbers.tif"
template <- raster(path_template)


# Read in raw bioclim data for sagebrush extent and set cellnumbers as rownames
# (tc stands for 'targetcells')
tc1 <- read_csv("data_processed/interpolation_data/clim_for_interpolation.csv",
                col_types = cols(site_id = 'd'),
                show_col_types = FALSE)

# Prep data -------------------------------------------------------------
# prep input data files for interpolation

# * Target cell data ------------------------------------------------------
# prepare the bioclim data

tc1 <- as.data.frame(tc1) %>%  # col names for tibbles are deprecated
  # creating dummy x and y coordinates (using projected data here
  # so lat/lon isn't approprate anymore), but location identifying x,y
  # is needed by multivarmatch()
  mutate(x = cellnumber,
         y = cellnumber)

# need x and y coordinates for the rmultivarmatch
coords <- raster::coordinates(template)
tc1$x <- coords[, 'x'][tc1$cellnumber]
tc1$y <- coords[, 'y'][tc1$cellnumber]

rownames(tc1) <- tc1$cellnumber # needed for multivarmatch fun

# limit to just the 6 variables of interest (keep lat/long)
# these are the same bioclim variables R.R. used
bioclim_vars <- c("bio1", "bio4", "bio9", "bio12", 
                  "bio15", "bio18")

tc2 <- tc1[, c("cellnumber", 'site_id', 'x', 'y', bioclim_vars)]



# criteria for matchingvars function (here using 10% of range of data)
# this is for scaling the variables
criteria <- map_dbl(tc2[, bioclim_vars], function(x) {
  (max(x) - min(x))*0.1
})

criteria


# * subset cell data ------------------------------------------------------
# location and climate data for the 200 sites where simulations were actually
# done


# sc stands for 'subset cell'
sc1 <-  tc2 %>% 
  filter(!is.na(site_id))

# check (should be 200 sites)
if(nrow(sc1) != 200 | any(duplicated(sc1$site_id))) {
  warning("Problem with dataframe")
}


# determine w/ sites have c4 grasses --------------------------------------
# creating a dataframe of where c4 grasses are present under current conditions
# --this is to be used in a later script

# not in use at the moment:

# sites_c4_present <- pft5_bio2 %>% 
#   filter(RCP == "Current", graze == "Light", PFT == "C4Pgrass",
#          c4 == "c4on") %>% 
#   mutate(C4Pgrass = ifelse(biomass > 0, "present", "absent")) %>% 
#   dplyr::select(site, C4Pgrass) %>% 
#   inner_join(site_nums2, by = c("site" = "site_id"))
# 
# stopifnot(nrow(sites_c4_present) == 200)
# 
# write_csv(sites_c4_present, 
#           "data_processed/site-num_C4Pgrass-presence_c4off.csv")

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
pft5_bio_w1 <- pft5_bio1 %>% 
  ungroup() %>% 
  # c4off and on are the same for RCP current, so don't want redundant 
  # up-scaling
  #filter(!(c4 == "c4off" & RCP == "Current")) %>% 
  mutate(id = paste(run, PFT, "biomass", id, GCM, sep = "_")) %>% 
  dplyr::select(site, id, biomass) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "biomass")

# joining in cell numbers
pft5_bio_w2 <- join_subsetcells(step_dat = pft5_bio_w1, sc_dat = sc1)


# * wildfire --------------------------------------------------------------

fire_w1 <- fire0 %>% 
  ungroup() %>% 
  # don't want to bother upscale fire prob on runs where fire not simulated
  filter(!str_detect(run, 'fire0')) %>% 
  mutate(id = paste(run, "fire-prob", id, GCM, sep = "_")) %>% 
  dplyr::select(site, id, fire_prob) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "fire_prob") %>% 
  join_subsetcells(sc_dat = sc1)


# * climate ---------------------------------------------------------------
# climate data doesn't have a 'run' name attached b/ runs have same climate
# but date attached in case in future runs are based on a different climate
# dataset
clim_all_w1 <- clim_all2 %>% 
  pivot_longer(cols = c("MAP", "MAT")) %>% 
  ungroup() %>% 
  mutate(id = paste(name, "climate", RCP, years, GCM, '20230919', sep = "_")) %>% 
  dplyr::select(site, id, value) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "value")

# joining in cell numbers
clim_all_w2 <- join_subsetcells(step_dat = clim_all_w1, sc_dat = sc1)

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
match1 <- multivarmatch(
  matchingvars = dplyr::select(tc2, -site_id),
  # I'm not sure the site ID column is needed/doing anything here
  subsetcells = sc1[, c("x", "y", "site_id")],
  matchingvars_id = "cellnumber",
  subsetcells_id = "site_id",
  criteria = criteria,
  raster_template = template,
  subset_in_target = TRUE
)



# *plotting interpolation quality -----------------------------------------

interp <- template

interp_ordered <- template

tmp_site_id <- match1 %>% 
  mutate(subset_cell = as.numeric(subset_cell)) %>% 
  left_join(sc1, by = c("subset_cell" = "cellnumber")) %>% 
  pull(site_id)

# for each pixels this raster shows which stepwat site the data was drawn
# from for interpolation
interp_ordered[as.numeric(match1$target_cell)] <- tmp_site_id
writeRaster(interp_ordered, "data_processed/interpolation_data/interp_locations_200sites.tif",
            overwrite = TRUE)

set.seed(1234)
# map of where site_ids are interpolated to
interp[as.numeric(match1$target_cell)] <- tmp_site_id %>% 
  factor() %>% 
  forcats::fct_shuffle()

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, 'RdBu')))(200)
jpeg(paste0("figures/interpolation_quality/interpolation_areas_", date, ".jpeg"),
     res = 600, units = 'in', width = 8, height = 8)
  plot(interp, main = 'interpolation areas (site numbers jumbled)',
       col = pal)
dev.off()


# proportion cells with decent matching quality
mean(match1$matching_quality < 1.5) # ~93%

match_quality <- template
match_quality[as.numeric(match1$target_cell)] <- match1$matching_quality

jpeg(paste0("figures/interpolation_quality/matching_quality_", date, ".jpeg"),
     res = 600, units = 'in', width = 8, height = 8)
plot(match_quality, main = 'matching quality', breaks = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 10),
     col = rev(RColorBrewer::brewer.pal(10, 'RdBu'))[3:10])
dev.off()


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


# if file was corrupted, and really small then re-create
min_size <- file.size(path_template)*0.2 # 0.5 multiplier is arbitrary, and maybe not strict enough

# which columns haven't already been upscaled
todo2 <- which_todo(df = pft5_bio_w2,
                   path = "data_processed/interpolated_rasters/biomass",
                   pattern = ".tif$",
                    min_size = min_size,
                   rerun = rerun)

todo2

# num cores, this includes logical cores (threads)
num.cores <- parallel::detectCores() 

registerDoParallel(num.cores)

# if all have been upscaled then todo just includes 'cellnumbers'
if(length(todo2) > 1) {
  

if(test_run) {
  todo2 <- todo2[1:3]
}
pft5_bio_w3 <- pft5_bio_w2[, todo2]



# making list of columns to use each time through dopar loop
vecs_l <- col_nums_parallel(pft5_bio_w3, by = 50)

print('biomass start')
print(Sys.time())

# seperately running interpolation on different sets of columns
foreach (x = vecs_l) %dopar% {
  rMultivariateMatching::interpolatePoints(
    matches = match1,
    output_results = pft5_bio_w3[, c(1, x)], 
    exclude_poor_matches = FALSE,
    subset_cell_names = "subset_cell",
    quality_name = "matching_quality",
    matching_distance = 1.5,
    raster_template = template,
    plotraster = FALSE,
    saveraster = TRUE,
    filepath = "./data_processed/interpolated_rasters/biomass",
    overwrite = TRUE
  )
}

print('biomass done')
print(Sys.time())
}
# * fire ------------------------------------------------------------------

# which columns haven't already been upscaled
todo3 <- which_todo(df = fire_w1,
                    path = "data_processed/interpolated_rasters/fire",
                    pattern = ".tif$",
                    min_size = min_size,
                    rerun = rerun)

# if all have been upscaled then todo just includes 'cellnumbers'
if(length(todo3) > 1) {
if(test_run) {
  todo3 <- todo3[1:3]
}
fire_w2 <- fire_w1[, todo3]


# making list of columns to use each time through dopar loop
vecs_l2 <- col_nums_parallel(fire_w2, by = 20)

print('fire start')
print(Sys.time())

# seperately running interpolation on different sets of columns
foreach (x = vecs_l2) %dopar% {
  rMultivariateMatching::interpolatePoints(
    matches = match1,
    output_results = fire_w2[, c(1, x)], 
    exclude_poor_matches = FALSE,
    subset_cell_names = "subset_cell",
    quality_name = "matching_quality",
    matching_distance = 1.5,
    raster_template = template,
    plotraster = FALSE,
    saveraster = TRUE,
    filepath = "./data_processed/interpolated_rasters/fire",
    overwrite = TRUE
  )
}
print('fire end')
print(Sys.time())

# When you're done, clean up the cluster
stopImplicitCluster()
}


# * climate ---------------------------------------------------------------

if (run_climate) {
  rMultivariateMatching::interpolatePoints(
    matches = match1,
    output_results = clim_all_w2, 
    exclude_poor_matches = FALSE,
    subset_cell_names = "subset_cell",
    quality_name = "matching_quality",
    matching_distance = 1.5,
    raster_template = template,
    plotraster = FALSE,
    saveraster = TRUE,
    filepath = "./data_processed/interpolated_rasters/climate",
    overwrite = TRUE
  )
}
