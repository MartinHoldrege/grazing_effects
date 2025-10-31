# Purpose:
# calculate average of fire probability drivers (for each ecoregion) for each GCM
# downstream this will be used for helping understand the widespread
# fire probabilities across GCMs. Within each ecoregion also seperately
# calculating means by SEI class ('c3')

# params ------------------------------------------------------------------

source('src/params.R')
v_interp <- opt$v_interp
run <- opt$run
vr <- opt$vr
vr_name <- opt$vr_name
pfts <- c('Sagebrush', 'Aherb', 'Pherb')

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')
source('src/SEI_functions.R')

# read in data ------------------------------------------------------------

# weights for sites (how many cells each site interpolated to by ecoregion
# and sei class)
# output from "scripts/05_interpolation_weights.R"
wt_c3eco <- read_csv(
  paste0('data_processed/interpolation_data/interpolation_weights_c3_', 
         v_interp, '_', vr, '.csv'))

# list dataframes of stepwat2 output created in "scripts/02_summarize_bio.R" 
bio <- readRDS('data_processed/site_means/summarize_bio.RDS') 

# combine dataframes ------------------------------------------------------

clim_long <- bio$clim_all2 %>% 
  pivot_longer(cols = c('MAP', 'MAT', 'psp'),
               names_to = 'variable')

bio_long <- bio$pft5_bio1 %>%
  ungroup() %>% 
  select(-indivs, -utilization, -id) %>% 
  filter(PFT %in% pfts,
         run %in% !! run) %>% 
  rename(value = biomass,
         variable = PFT) %>% 
  mutate(variable = as.character(variable))

comb1 <- bind_rows(bio_long, clim_long)

# weighted means ----------------------------------------------------------
# weighted by area

# mean across SEI class and ecoregions
comb2 <- comb1 %>%
  select(-run) %>% 
  right_join(wt_c3eco,   relationship = "many-to-many",
             by = 'site') 

#summarize across sei class
comb2_eco <- comb2 %>% 
  group_by(site, years, RCP, graze, GCM, region, variable) %>% 
  summarise(
    # not really doing anything, taking the mean of equal values
    # at a site (i.e. if a site was interpolated to areas with differen SEI classes)
    value = mean(value),
    weight = sum(weight),
    .groups = 'drop') 
means_c3eco1 <- comb2 %>% 
  group_by(years, RCP, graze, GCM, region, variable, c3) %>% 
  # calculated weighted percentiles
  summarize_weighted(varname = 'value')

# mean across regions
means_eco1 <- comb2_eco %>% 
  group_by(years, RCP, graze, GCM, region, variable) %>% 
  summarize_weighted(varname = 'value')


# change relative to current ----------------------------------------------


delta1 <- comb2_eco %>% 
  mutate(region = region_factor(region)) %>% 
  # for current conditions just want the site value 
  filter(RCP == 'Current', region == levels(region)[1]) %>% 
  select(-weight, -years, -RCP, -GCM, -region) %>% # only care about the weight under future conditions
  right_join(filter(comb2_eco, RCP != 'Current'),
             by = join_by(variable, site, graze),
             suffix = c('_cur', '_fut'),
             relationship = "many-to-many") %>% 
  mutate(delta = value_fut - value_cur)

smry_delta_eco1 <- delta1 %>% 
  group_by(years, RCP, graze, GCM, region, variable) %>% 
  summarize_weighted(varname = 'delta')
  

# save output -------------------------------------------------------------

write_csv(means_eco1, paste0('data_processed/raster_means/', run, vr_name,
                             '_fire-driver-means_by-ecoregion.csv'))

write_csv(smry_delta_eco1, paste0('data_processed/raster_means/', run, vr_name,
                             '_delta-fire-driver_by-ecoregion.csv'))

write_csv(means_c3eco1, paste0('data_processed/raster_means/', run, '_', vr,
                               '_fire-driver-means_by-ecoregion-c3.csv'))

