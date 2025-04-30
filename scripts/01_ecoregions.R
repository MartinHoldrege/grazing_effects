# Purpose: Load level 3 ecoregions, and dissolve and extract the
# regions for our analysis (to create 4 regions that cover our region,
# keeping the regions used in the SCD but splitting the intermountain
# west region into two), used in load_wafwa_ecoregions() function

# Author: Martin Holdrege

# Date Started: April 29, 2025


# parameters --------------------------------------------------------------

# regions as listed in the SCD (Doherty et al. 2022)
wafwa_names <- list(
  'Southern Great Basin' = c("Central Basin and Range", 
                             "Mojave Basin and Range", 
                             "Sonoran Basin and Range"),
  'Great Plains' = c("High Plains", "Southwestern Tablelands", 
                     "Northwestern Glaciated Plains", 
                     "Northwestern Great Plains"),
  'Intermountain West' = c("Cascades", "Sierra Nevada", "Eastern Cascades Slopes and Foothills", "Columbia Plateau", "Blue Mountains", "Snake River Plain", "Northern Rockies", "Idaho Batholith", "Middle Rockies", "Wyoming Basin", "Wasatch and Uinta Mountains", "Colorado Plateaus", "Southern Rockies", "Arizona/New Mexico Plateau", "Arizona/New Mexico Mountains", "Northern Basin and Range")

)


# dependencies ------------------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
source('src/mapping_functions.R')
# read in data ------------------------------------------------------------

# downloaded from https://www.epa.gov/eco-research/ecoregions-north-america
lev3a <- read_sf("data_raw/ecoregion_level_3/us_eco_l3.shp")


# filter data -------------------------------------------------------------

# check
L3NAMES <- unlist(unname(wafwa_names))
stopifnot(L3NAMES %in% lev3a$US_L3NAME)


# 'eastern' intermountain
east_im <- c('6.2.10', '10.1.4', '6.2.13', '10.1.6', '6.2.14', 
             '10.1.7', '13.1.1')

# * View regions ------------------------------------------------------------

im <- lev3a %>% 
  filter(US_L3NAME %in% wafwa_names$`Intermountain West`) %>% 
  group_by(US_L3CODE, US_L3NAME, NA_L3CODE) %>% 
  summarise(geometry = sf::st_combine(geometry), .groups = 'drop')


ggplot(im) +
  geom_sf(fill = "lightgray", color = "black") +
  geom_sf_text(aes(label = NA_L3CODE), size = 3) +
  theme_minimal()


im %>% 
  filter(!NA_L3CODE %in% east_im) %>% 
  ggplot() +
  geom_sf(fill = "lightgray", color = "black") +
  geom_sf_text(aes(label = NA_L3CODE), size = 3) +
  theme_minimal()


# *adding wafwa region names ---------------------------------------------------

in_which_element <- function(x, l) {
  for(i in 1:length(l)) {
    if(x %in% l[[i]]) return(names(l)[i])
  }
  return(NA_character_)
}

lev3b <- lev3a %>% 
  filter(US_L3NAME %in% L3NAMES)

lev3c <- lev3b %>% 
  mutate(wafwa_NAME = map_chr(US_L3NAME, in_which_element, l = wafwa_names),
         # splitting intermountain west
         region = case_when(
           wafwa_NAME == "Intermountain West" & NA_L3CODE %in% east_im ~ 'Eastern Intermountain',
           wafwa_NAME == "Intermountain West" ~ 'Western Intermountain',
           TRUE ~ wafwa_NAME
         )
  ) 


# dissolving --------------------------------------------------------------

eco1 <- lev3c %>% 
  group_by(wafwa_NAME, region) %>% 
  summarise(geometry = st_union(geometry), .groups = 'drop')

eco2 <- st_transform(eco1, crs_scd)

ggplot(eco2) +
  geom_sf(fill = "lightgray", color = "black") +
  geom_sf_text(aes(label = region), size = 3) +
  theme_minimal()
# write output ------------------------------------------------------------

st_write(eco2, 'data_processed/ecoregions/four_regions_v1.gpkg',
         overwrite = TRUE, append = FALSE)


