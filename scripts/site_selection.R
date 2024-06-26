# Purpose: Select sites for doing stepwat simulations, with the goal of 
# representing the whole study area well (1 km grid cells containing >50% SCD 30 m cells)

# Author: Martin Holdrege

# Started: June 17, 2024

# parameters --------------------------------------------------------------

scd_prop <- 0.5 # proportion of grid cell that contains scd cells, to be part of study area
period <- '1991-2020' # period from which climate normals were calculated

# criteria used in Renne et al 2024 
criteria <- c(bio1 = 1.55, 
                 bio4 = 61.1, 
                 bio9 = 3.53, 
                 bio12 = 84.9, 
                 bio15 = 9.21, 
                 bio18 = 33.8
)
crit_descript <- 'Renne et al matching criteria'
v <- 'v1'
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
theme_set(theme_bw())
library(rMultivariateMatching)
source('src/mapping_functions.R')
source("src/matching_functions.R")

# load data ---------------------------------------------------------------

template <- rast("data_processed/interpolation_data/cellnumbers.tif")

# daymet data
r_day1 <- rast(paste0("data_raw/daymet_v4_monthly_normals_", period, ".tif"))
day1 <- read_csv(paste0("data_processed/interpolation_data/clim_for_interp_daymet-v4_",
                        period, ".csv"),
                 col_types = paste0(rep('n', 23), collapse = ""))

# * grid-met --------------------------------------------------------------

buffer <- 0
# file created in 00_gridmet_normals.R
r_grid1 <- rast(paste0("data_processed/interpolation_data/gridmet_bioclim-vars_19912020_",
                     buffer, "buffer.tif"))

# preparing dataframes ----------------------------------------------------

grid1 <- as.data.frame(r_grid1)
grid1$cellnumber <- row.names(grid1)
day2 <- day1 %>% 
  filter(prop_scd > scd_prop) # filter to mask

# normalize variables
grid_norm <- apply_criteria(grid1, criteria = criteria)
day_norm <- apply_criteria(day2, criteria = criteria)
# kpoints -----------------------------------------------------------------
# finding optimum number of sites 

# columns used for evenly spacing initial 'centers' for kmeans
# MAT, MAP, ppt warmest quarter
cols_for_centers <- c("bio1", "bio12", "bio18")

# test <- kpoints2(target = day_norm, 
#          potential_sites = grid_norm,
#          criteria_vars = names(criteria),
#          digits = 4,
#          k = 50,
#          verbose = TRUE)

kseq <- seq(25, 800, by = 25)

set.seed(123)
kresults <- map(kseq, \(k) {
  print(k)

  kpoints2(target = day_norm, 
    potential_sites = grid_norm,
    criteria_vars = names(criteria),
    cols_for_centers = cols_for_centers,
    digits = 4,
    k = k,
    verbose = TRUE)
  })

# saving b/ previous step is slow
saveRDS(kresults, 
        paste0('data_processed/site_selection/kpoints_', v, '.rds'))

kresults <- readRDS(paste0('data_processed/site_selection/kpoints_', v, '.rds'))

# * existing coverage -----------------------------------------------------

cuttoff <- 1

cov200 <- prop_coverage(day_norm,
              sites = day_norm[!is.na(day_norm$site_id), ],
              criteria_vars = names(criteria),
              cuttoff = cuttoff)

cov200 <- tibble(coverage = cov200,
                 k = 200)

# * kpoints figs ----------------------------------------------------------
cap1 <- paste0('sites selected from gridmet (', buffer, ' buffer).\n',
               'Study area = 1km pixels with >', scd_prop*100, '% SCD pixels.\n',
               crit_descript)
cov1 <- map(kresults, \(l) {
  tibble(coverage = l$coverage,
         converged = l$ifault == 0,
         k = nrow(l$sites))
}) %>% 
  bind_rows()

ggplot(cov1, aes(k, coverage)) +
  geom_point(aes(color = converged)) +
  geom_point(data = cov200, aes(color = "200 published sites")) +
  labs(y = 'Proportion of daymet cells covered by k sites',
       caption = cap1) +
  scale_color_manual(values = c('blue', 'gray', 'black'),
                     name = 'kmeans converged')
  


