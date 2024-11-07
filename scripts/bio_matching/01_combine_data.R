# Purpose: Compare stepwat and RAP biomass, at the 200 simulation sites 
# so that some sort of quantile matching relationship can determined
# so that stepwat biomass is translated to biomass levels used for 
# fitting the Holdrege et al. 2024 fire equation

# Author: Martin Holdrege

# Oct 3, 2024


# params ------------------------------------------------------------------

version_interp <- 'v3' # interpolation version (see 03_interpolate.R)
run <- "fire1_eind1_c4grass1_co20_2311"
runs <- c("fire1_eind1_c4grass1_co20_2311", "fire0_eind1_c4grass1_co20")

# matching quality cutoff, this decides
# the cells from which both stepwat and rap data will be pulled from
qual_cutoff <- 0.5 
RCP <- 'Current'
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(patchwork)
source("src/general_functions.R")
source('src/fig_params.R')
source('src/fig_functions.R')
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

# data from the Holdrege et al. 2024 publication
rap_df1 <- read_csv("../cheatgrass_fire/data_processed/data_publication/Wildfire_Climate_Biomass_Data.csv")

cell_nums <- rast("../cheatgrass_fire/data_processed/data_publication/cell_nums.tif")

# created in grazing_effects/scripts/02_summarize_bio.R
pft5_bio2 <- readRDS('data_processed/site_means/summarize_bio.RDS')$pft5_bio2

sites1 <- read_csv("data_raw/site_locations.csv")
clim1 <- read_csv("data_processed/site_means/dbWeather_200sites.csv")

# partial dependence plot values from fire equation
# file created in cheatgrass_fire/scripts/06_figures_pdp_vip_quant.R
df_pdp2 <- readRDS("../cheatgrass_fire/data_processed/df_pdp2.rds");

# * interpolated stepwat --------------------------------------------------
PFTs <- c('Aherb', 'Pherb')

suffix <- paste0(version_interp, ".*", RCP, ",*.tif")
pattern <- paste0("(", paste0(runs, collapse = '|'), ")", suffix)

paths <- list.files("data_processed/interpolated_rasters/biomass/",
                    pattern = pattern,
                    full.names = TRUE)
basename(paths)
r_sw1 <- rast(paths)

r_qual <- rast(paste0("data_processed/interpolation_quality/matching_quality",
                 version_interp, ".tif"))

# vectors etc -------------------------------------------------------------

pft_lookup <- c('pfgAGB' = 'Pherb','afgAGB' = 'Aherb')

# prepare pdp predictions -------------------------------------------------

df_pdp3 <- df_pdp2 %>% 
  filter(inter_var == 'Mean prediction',
         variable %in% c("pfgAGB", 'afgAGB')) %>% 
  mutate(PFT = pft_lookup[variable]) %>% 
  select(-variable)

# stepwat site level analysis ---------------------------------------------

# * select gridcells of rap data  ------------------------------------------

sites2 <- sf::st_as_sf(sites1, coords = c("X_WGS84", "Y_WGS84"),
                       crs = 'EPSG:4326') %>% 
  st_transform(crs = crs(cell_nums)) %>% 
  arrange(site_id)


# checking that ID returned from extract will match site num
stopifnot(sites2$site_id == 1:max(sites2$site_id)) 
cell_nums_sites <- extract(cell_nums, sites2[, ]) %>% 
  rename(site = ID) %>% 
  as_tibble()

site_ids <- cell_nums_sites$site[!is.na(cell_nums_sites$cell_num)]

# filter site level data

sw_site_bio1 <- pft5_bio2 %>% 
  calc_aherb(group_cols = c('run', 'years', 'RCP', 'graze', 'id', 
                            'site')) %>% 
  filter(.data$run %in% !!runs,
         RCP == "Current",
         site %in% site_ids,
         PFT %in% unname(pft_lookup)) %>% 
  mutate(fire = str_extract(run, 'fire\\d'))

# * evaluate clim of stepwat sites ------------------------------------------
# see if the subset of stepwat sites, that overlap with the rap data from the
# fire data are representative of all sites

clim2 <- clim1 %>% 
  mutate(dataset_overlap = Site_id %in% site_ids)



# ** extract RAP data --------------------------------------------------------

rap_df2 <- rap_df1 %>% 
  select(year, cell_num, pfgAGB, afgAGB)

rap_sub1 <- cell_nums_sites %>% 
  filter(.data$site %in% site_ids) %>% 
  left_join(rap_df2, by = 'cell_num')

rap_sub2 <- rap_sub1 %>%  
  pivot_longer(cols = c('pfgAGB', 'afgAGB'),
               values_to = 'biomass',
               names_to = 'PFT') %>% 
  mutate(PFT = pft_lookup[PFT])

# pulling together interpolation level data -------------------------------

# identify stepwat grid-cells w/ good matching, and that are also cells
# included in the fire manuscript
cell_nums2 <- extend(cell_nums, r_qual)
mask0 <- r_qual < qual_cutoff & !is.na(cell_nums2)
mask <- classify(mask0, rcl = cbind(FALSE, NA))

r_sw2 <- mask(r_sw1, mask)
df_pherb <- r_sw2[[str_subset(names(r_sw2), 'Pherb')]]
df_aherb <- r_sw2
cell_nums_m <- mask(cell_nums2, mask)
cell_nums_mv <- as.numeric(cell_nums_m[!is.na(cell_nums_m)]) # vector

rap_df_m <- rap_df2 %>% 
  filter(cell_num %in% cell_nums_mv) %>% 
  pivot_longer(cols = c("pfgAGB", 'afgAGB'),
               names_to = 'PFT',
               values_to = 'biomass') %>% 
  mutate(PFT = pft_lookup[PFT]) %>% 
  select(-cell_num, -year)

df_sw1 <- as_tibble(r_sw2) 

info1 <- create_rast_info(r_sw2)

df_sw2 <- df_sw1 %>% 
  pivot_longer(cols = everything(),
               names_to = 'id',
               values_to = 'biomass') %>% 
  left_join(info1, by = 'id') %>% 
  mutate(fire = str_extract(run, 'fire\\d'))

tmp1 <- df_sw2 %>% 
  select(run, fire, biomass, PFT, graze) %>% 
  mutate(dataset = 'interpolated') 

sw_comb <- sw_site_bio1 %>% 
  select(run, fire, biomass, PFT, graze) %>% 
  mutate(dataset = 'site level') %>% 
  bind_rows(tmp1)

