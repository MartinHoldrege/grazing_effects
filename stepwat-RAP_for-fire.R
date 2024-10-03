# Purpose: Compare stepwat and RAP biomass, at the 200 simulation sites 
# so that some sort of quantile matching relationship can determined
# so that stepwat biomass is translated to biomass levels used for 
# fitting the Holdrege et al. 2024 fire equation

# Author: Martin Holdrege

# Oct 3, 2024


# params ------------------------------------------------------------------



# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(patchwork)
# read in data ------------------------------------------------------------

# data from the Holdrege et al. 2024 publication
rap_df1 <- read_csv("../cheatgrass_fire/data_processed/data_publication/Wildfire_Climate_Biomass_Data.csv")

cell_nums <- rast("../cheatgrass_fire/data_processed/data_publication/cell_nums.tif")

# created in grazing_effets/scripts/02_summarize_bio.R
pft5_bio2 <- readRDS(file.path(dir_graze, 'data_processed/site_means/pft5_bio2.RDS'))

sites1 <- read_csv("data_raw/site_locations.csv")
clim1 <- read_csv("data_processed/site_means/dbWeather_200sites.csv")
# select gridcells of rap data  ------------------------------------------

sites2 <- sf::st_as_sf(sites1, coords = c("X_WGS84", "Y_WGS84"),
                       crs = 'EPSG:4326') %>% 
  st_transform(crs = crs(cell_nums)) %>% 
  arrange(site_id)

plot(cell_nums)
plot(sites2$geometry, add = TRUE)

# checking that ID returned from extract will match site num
stopifnot(sites2$site_id == 1:max(sites2$site_id)) 
cell_nums_sites <- extract(cell_nums, sites2[, ]) %>% 
  rename(site = ID) %>% 
  as_tibble()

site_ids <- cell_nums_sites$site[!is.na(cell_nums_sites$cell_num)]


# evaluate clim of stepwat sites ------------------------------------------
# see if the subset of stepwat sites, that overlap with the rap data from the
# fire data are representative of all sites

clim2 <- clim1 %>% 
  mutate(dataset_overlap = Site_id %in% site_ids)

g <- ggplot(clim2, aes(color = dataset_overlap)) +
  theme_bw() +
  scale_color_manual(name = 'STEPWAT site falls on \nHoldrege et al 2024 1km cell',
                     values = c('darkgrey', 'blue' ))

g1 <- g + geom_point(aes(MAT_C, MAP_mm))

g2 <- g + geom_point(aes(CorrTP2, MAP_mm)) +
  labs(x = 'T-P correlation')
g3 <- g + geom_point(aes(CorrTP2, MAT_C)) +
  labs(x = 'T-P correlation')
p <- g1 + g2 + g3 + plot_layout(guides = 'collect',
                           ncol = 2) 
p&theme(legend.position = 'top')

# * extract RAP data --------------------------------------------------------

rap_df2 <- rap_df1 %>% 
  select(year, cell_num, pfgAGB, afgAGB)

rap_sub1 <- cell_nums_sites %>% 
  filter(.data$site %in% site_ids) %>% 
  left_join(rap_df2, by = 'cell_num')

hist(rap_sub1$pfgAGB, xlim = c(0, 200), breaks = 200)
hist(rap_df2$pfgAGB, xlim = c(0, 200))


