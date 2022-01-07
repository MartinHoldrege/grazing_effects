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

# Read in 200 points
pts <- read.csv(file.path(path1, "Simulated_sites_final200.csv"), row.names= 1)

# file including site numbers (1:200)
site_nums <- read_csv("data_raw/site_locations.csv")

# Read in raw bioclim data for sagebrush extent and set cellnumbers as rownames
# (tc stands for 'targetcells')
tc1 <- readRDS(file.path(path1, "bioclim_coreARTR_DayMet_trimmed.csv"))


# Prep inputs -------------------------------------------------------------
# input files for interpolation

# * Target cell data ------------------------------------------------------
# prepare the bioclim data

rownames(tc1) <- tc1$cellnumbers

# limit to just the 6 variables of interest (keep lat/long)
# these are the same bioclim variables R.R. used
bioclim_vars <- c("bioclim_01", "bioclim_04", "bioclim_09", "bioclim_12", 
                  "bioclim_15", "bioclim_18")

tc2 <- tc1[, c("cellnumbers", "x", "y", bioclim_vars)]

# criteria for matchingvars function (here using 10% of range of data)
# this is for scaling the variables
criteria <- map_dbl(tc2[, bioclim_vars], function(x) {
  (max(x) - min(x))*0.1
})


# * subset cell data ------------------------------------------------------
# location and climate data for the 200 sites where simulations were actually
# done

sc1 <-  inner_join(pts, site_nums, b = c("x" = "X_WGS84", "y" = "Y_WGS84"))
plot(pts$y, site_nums$Y_WGS84)

#Continue here: to join must first round, likely floating point issues.
sum(round(pts$x, 5) %in% round(site_nums$X_WGS84, 5))
# data on

