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
# rounding (so x, y used everwhere will match)
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
pts2
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

# stepwat2 outputs --------------------------------------------------------

# files created in 02_summarize_bio.R script
# this is data for the 5 plant functional types (PFTs) of interest

# * biomass (g/m2) ----------------------------------------------------------
# actual biomass (ie not difference)
# converting to wide format (w stands for 'wide'), 
# for using with interpolationPoints function below
pft5_bio_w1 <- pft5_bio2 %>% 
  mutate(id = paste(PFT, "biomass", id, sep = "_")) %>% 
  dplyr::select(site, id, biomass) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "biomass")

pft5_bio_w2 <- sc1[, c("cellnumbers", "site_id")] %>% 
  rename(site = site_id) %>% 
  inner_join(pft5_bio_w1, by = "site") %>% 
  dplyr::select(-site)

# rownames needed for interpolatePoints
rownames(pft5_bio_w2) <- pft5_bio_w2$cellnumbers

stopifnot(nrow(pft5_bio_w2) == 200) # check for join issues

# * biomass difference ----------------------------------------------------
# % scaled biomass difference

# converting to wide format
pft5_bio_d_w1 <- pft5_bio_d2 %>% 
  mutate(id = paste(PFT, "bio-diff", id, sep = "_")) %>% 
  dplyr::select(site,  id, bio_diff) %>% 
  pivot_wider(id_cols = "site",
              names_from = "id",
              values_from = "bio_diff")

pft5_bio_d_w2 <- sc1[, c("cellnumbers", "site_id")] %>% 
  rename(site = site_id) %>% 
  inner_join(pft5_bio_d_w1, by = "site") %>% 
  dplyr::select(-site)

rownames(pft5_bio_d_w2) <- pft5_bio_d_w2$cellnumbers

# identify matches --------------------------------------------------------
# ID matches from subset cells for all Target cells (i.e. calculates distance)

# check that subset_in_target is true, for fun call below
if(!all(rownames(sc1 %in% rownames(tc2)))) {
  warning("Not all subset cells in target cells")
}

# note--consider saving the figure of matching quality
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


# * biomass ---------------------------------------------------------------

interpolatePoints(
  matches = match1,
  output_results = pft5_bio_w2, 
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


# * biomass diff ----------------------------------------------------------
# scaled % biomass difference

interpolatePoints(
  matches = match1,
  output_results = pft5_bio_d_w2, 
  exclude_poor_matches = TRUE,
  subset_cell_names = "subset_cell",
  quality_name = "matching_quality",
  matching_distance = 1.5,
  raster_template = template,
  plotraster = FALSE,
  saveraster = TRUE,
  filepath = "./data_processed/interpolated_rasters/bio_diff",
  overwrite = TRUE
)
