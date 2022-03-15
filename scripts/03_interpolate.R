# Martin Holdrege

# Script started 1/6/2022

# Purpose: to use the rMultivariateMatching package to do nearest neighbor
# matching of the 200 simulated sites onto gridded weather data in the west,
# and use those matches to interpolate the STEPWAT2 biomass
# output across the west.

# Note: all the code below strongly borrows from code provided by Rachel Renne
# and uses her R package. 

# dependencies ------------------------------------------------------------

source("scripts/02_summarize_bio.R") # creates dataframes of stepwat2 output
library(raster)
library(rMultivariateMatching)
library(doParallel)

# read in data ------------------------------------------------------------

# folder where I put data that R.R. provided for interpolation of the 200 sites
path1 <- "data_raw/interpolation_data_RR"

# template raster
template <- raster(file.path(path1, "CoreARTR_combined_DayMet_cropped_trimmed.tif"))


# Read in 200 points (this data from RR has bioclim variables already)
pts <- read.csv(file.path(path1, "Simulated_sites_final200.csv"))

# file including site numbers (1:200) (for matching with simulation output)
site_nums <- read_csv("data_raw/site_locations.csv", show_col_types = FALSE)

# Read in raw bioclim data for sagebrush extent and set cellnumbers as rownames
# (tc stands for 'targetcells')
tc1 <- readRDS(file.path(path1, "bioclim_coreARTR_DayMet_trimmed.csv"))


# Prep data -------------------------------------------------------------
# prep input data files for interpolation

# * Target cell data ------------------------------------------------------
# prepare the bioclim data

rownames(tc1) <- tc1$cellnumbers

# limit to just the 6 variables of interest (keep lat/long)
# these are the same bioclim variables R.R. used
bioclim_vars <- c("bioclim_01", "bioclim_04", "bioclim_09", "bioclim_12", 
                  "bioclim_15", "bioclim_18")

tc2 <- tc1[, c("cellnumbers", "x", "y", bioclim_vars)]

# rounding (so x, y used will match)
digits <-  4
tc2 <- tc2 %>% 
  mutate(x = round(x, digits = digits),
         y = round(y, digits = digits))


# criteria for matchingvars function (here using 10% of range of data)
# this is for scaling the variables
criteria <- map_dbl(tc2[, bioclim_vars], function(x) {
  (max(x) - min(x))*0.1
})

criteria


# * subset cell data ------------------------------------------------------
# location and climate data for the 200 sites where simulations were actually
# done


# round 

pts2 <- pts %>% 
  rename(cellnumbers = X) %>% 
  mutate(x = round(x, digits = digits),
         y = round(y, digits = digits))

site_nums2 <- site_nums %>% 
  mutate(x = round(X_WGS84, digits = digits),
         y = round(Y_WGS84, digits = digits)) %>% 
  dplyr::select(-X_WGS84, -Y_WGS84)

# sc stands for 'subset cell'
sc1 <-  inner_join(pts2, site_nums2, by = c("x", "y"))
rownames(sc1) <- sc1$cellnumbers # rownames needed for multivarmatch function

# check that the join worked (should be 200 sites)
if(nrow(sc1) != 200 | any(duplicated(sc1$site_id))) {
  warning("Problem with join")
}


# determine w/ sites have c4 grasses --------------------------------------
# creating a dataframe of where c4 grasses are present under current conditions
# --this is to be used in a later script

sites_c4_present <- pft5_bio2 %>% 
  filter(RCP == "Current", graze == "Light", PFT == "C4Pgrass",
         c4 == "c4on") %>% 
  mutate(C4Pgrass = ifelse(biomass > 0, "present", "absent")) %>% 
  dplyr::select(site, C4Pgrass) %>% 
  inner_join(site_nums2, by = c("site" = "site_id"))

stopifnot(nrow(sites_c4_present) == 200)

write_csv(sites_c4_present, 
          "data_processed/site-num_C4Pgrass-presence_c4off.csv")

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
  filter(!(c4 == "c4off" & RCP == "Current")) %>% 
  # for now excluding RCP 4.5, to save on computation time
  filter(RCP != "RCP4.5") %>%  
  mutate(id = paste(c4, PFT, "biomass", id, GCM, sep = "_")) %>% 
  dplyr::select(site, id, biomass) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "biomass")

# joining in cell numbers
pft5_bio_w2 <- join_subsetcells(step_dat = pft5_bio_w1, sc_dat = sc1)


# *crossing threshold ----------------------------------------------------
# What is mildest (minimum) grazing treatment that causes biomass to go 
# below threshold, within a climate scenario

thresh_min_graze_w <- thresh_min_graze1 %>% 
  filter_rcp_c4() %>% 
  mutate(id = paste(c4, PFT, "min-graze", RCP, years, sep = "_")) %>% 
  dplyr::select(site, id, min_graze) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "min_graze") %>% 
  join_subsetcells(sc_dat = sc1) # joining in cell numbers 
  

# * wgcm bio diff ---------------------------------------------------------

# change in biomass from light to heavy grazing, within a climate scenario
pft5_d_wgcm_w <- pft5_d_wgcm %>% 
  filter_rcp_c4() %>% 
  filter(graze == "Heavy") %>% 
  mutate(id = paste(c4, PFT, "bio-diff-wgcm", id, sep = "_")) %>% 
  dplyr::select(site, id, bio_diff) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "bio_diff") %>% 
  join_subsetcells(sc_dat = sc1) 

# change in biomass from heavy to moderate or light grazing, within a scenario
pft5_d_wgcm_heavy_w <- pft5_d_wgcm_heavy %>% 
  filter_rcp_c4() %>% 
  filter(graze %in% c("Light", "Moderate")) %>% 
  mutate(id = paste(c4, PFT, "bio-diff-wgcm-heavy", id, sep = "_")) %>% 
  dplyr::select(site, id, bio_diff) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "bio_diff") %>% 
  join_subsetcells(sc_dat = sc1) 


# * bio diff --------------------------------------------------------------
# change in biomass from current to RCP 8.5, under a given grazing scenario/
# this is the % change scaled to current max biomass. 

pft5_bio_d2_w <- pft5_bio_d2 %>% 
  filter_rcp_c4(current = FALSE) %>% 
  filter(graze %in% c("Light", "Heavy")) %>% 
  # here 'wgraze' means bio-diff within a grazing treatment.
  # this is the median across GCMS
  mutate(id = paste(c4, PFT, "bio-diff-wgraze", id, 'median', sep = "_")) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "bio_diff") %>% 
  join_subsetcells(sc_dat = sc1) 

# identify matches --------------------------------------------------------
# ID matches from subset cells for all Target cells (i.e. calculates distance)

# check that subset_in_target is true, for fun call below
if(!all(rownames(sc1 %in% rownames(tc2)))) {
  warning("Not all subset cells in target cells")
}

# note--consider saving the figure of matching quality
# this only takes ~ min to run
match1 <- multivarmatch(
  matchingvars = tc2,
  # I'm not sure the site ID column is needed/doing anything here
  subsetcells = sc1[, c("x", "y", "site_id")],
  matchingvars_id = "cellnumbers",
  subsetcells_id = "site_id",
  criteria = criteria,
  raster_template = template,
  subset_in_target = TRUE
)

head(match1)

# proportion cells with decent matching quality
mean(match1$matching_quality < 1.5) # ~94%

# interpolation ---------------------------------------------------------
# use matches from multivarmatch to interpolate the STEPWAT2 output
# across all the grid cells

# list of dataframes to upscale
output_results_l <- list(
  pft5_bio_d2_w = pft5_bio_d2_w, # within graze bio-diff
  pft5_d_wgcm_w = pft5_d_wgcm_w, # wgcm bio-diff (relative to light grazing)
  pft5_d_wgcm_heavy_w = pft5_d_wgcm_heavy_w,# wgcm bio-diff (rel to heavy graze)
  thresh_min_graze_w = thresh_min_graze_w #crossing threshold
)

# locations where the upscaled rasters will be put
output_paths <- c(
  pft5_bio_d2_w = "./data_processed/interpolated_rasters/bio_diff",
  pft5_d_wgcm_w = "./data_processed/interpolated_rasters/bio_diff",
  pft5_d_wgcm_heavy_w = "./data_processed/interpolated_rasters/bio_diff",
  thresh_min_graze_w = "./data_processed/interpolated_rasters/min_graze"
)

stopifnot(names(output_paths) == names(output_results_l))

# select dataframes to upscale (i.e. code is slow so 
# if other dfs have been previously, they don't need to be re-run)
names_to_upscale <- "pft5_d_wgcm_heavy_w"
# names_to_upscale <- names(output_paths) # un-comment if want to upscale all

# interpolate
map2(output_results_l[names_to_upscale], 
     output_paths[names_to_upscale], 
     function(df, path) {
      rMultivariateMatching::interpolatePoints(
        matches = match1,
        output_results = df, 
        exclude_poor_matches = TRUE,
        subset_cell_names = "subset_cell",
        quality_name = "matching_quality",
        matching_distance = 1.5,
        raster_template = template,
        plotraster = FALSE,
        saveraster = TRUE,
        filepath = path,
        overwrite = TRUE
      )
})

# * biomass ---------------------------------------------------------------

# on my computer to run all PFTs/scenarios/GCMs would take 21.5 hours!, 
# when run in parallel should take ~6.5 hours, 3.26 hrs if RCP4.5 is not included. 
# It actually took 3 hrs, 40 min, as it is now. 

# I'm using ideas from here:
# https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
# to run this in parallel

# currently not re-running this code here (change if inputs are updated)
if (FALSE){
  
# num cores, this includes logical cores (threads)
num.cores <- parallel::detectCores() 

registerDoParallel(num.cores)

# making list of columns to use each time through dopar loop
n <- ncol(pft5_bio_w2)
by = 50 # each time through loop 50 variables will be up-scaled
vecs <- seq(from = 2, to = n, by = by)
vecs_l <- map(vecs, function(x) {
  to <- x + by - 1
  
  # for the last set of columns
  if(to > n) {
    to <- n
  }
  x:to # column numbers
})

# check that all columns accounted for
stopifnot(unlist(vecs_l) == 2:n)


print(Sys.time())

# seperately running interpolation on different sets of columns
foreach (x = vecs_l) %dopar% {
  rMultivariateMatching::interpolatePoints(
    matches = match1,
    output_results = pft5_bio_w2[, c(1, x)], 
    exclude_poor_matches = TRUE,
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
print(Sys.time())

# When you're done, clean up the cluster
stopImplicitCluster()

}