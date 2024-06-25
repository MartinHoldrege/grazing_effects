# Purpose: Select sites for doing stepwat simulations, with the goal of 
# representing the whole study area well (1 km grid cells containing >50% SCD 30 m cells)

# Author: Martin Holdrege

# Started: June 17, 2024

# parameters --------------------------------------------------------------

scd_prop <- 0.5 # proportion of grid cell that contains scd cells, to be part of study area
period <- '1991-2020' # period from which climate normals were calculated

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(rMultivariateMatching)
source('src/mapping_functions.R')

# load data ---------------------------------------------------------------

template <- rast("data_processed/interpolation_data/cellnumbers.tif")

# daymet data
day1 <- day1 <- rast(paste0("data_raw/daymet_v4_monthly_normals_", period, ".tif"))

# * grid-met --------------------------------------------------------------

buffer <- 0
# file created in 00_gridmet_normals.R
grid1 <- rast(paste0("data_processed/interpolation_data/gridmet_bioclim-vars_19912020_",
                     buffer, "buffer.tif"))


# testing


# Load targetcells data for Target Cells
data(targetcells)

# Create data frame of potential matching variables for Target Cells
allvars <- makeInputdata(targetcells)

# Restrict data to matching variables of interest
matchingvars <- allvars[,c("cellnumbers","x","y","bioclim_01","bioclim_04",
                           "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]

# Create raster_template
raster_template <- targetcells[[1]]

# Create vector of matching criteria
criteria <- c(0.7,42,3.3,66,5.4,18.4)

# Verify stopping criteria for 200 points
# Note: n_starts should be >= 10, it is 1 here to reduce run time.
results1 <- kpoints(matchingvars,criteria = criteria, klist = 200,
                    n_starts = 1, min_area = 50, iter = 50,
                    raster_template = raster_template,
                    verify_stop = TRUE, savebest = FALSE)

