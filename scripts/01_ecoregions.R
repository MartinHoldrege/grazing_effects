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
  'Intermountain West' = c("Cascades", "Sierra Nevada", 
                           "Eastern Cascades Slopes and Foothills", 
                           "Columbia Plateau", "Blue Mountains", 
                           "Snake River Plain", "Northern Rockies", 
                           "Idaho Batholith", "Middle Rockies", 
                           "Wyoming Basin", "Wasatch and Uinta Mountains", 
                           "Colorado Plateaus", "Southern Rockies", 
                           "Arizona/New Mexico Plateau",
                           "Arizona/New Mexico Mountains", 
                           "Northern Basin and Range")

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



# v1.1 regions ------------------------------------------------------------

# regions version 1.1
r1 <- list(# Western intermountain
  'NW Deserts' = c('10.1.3', '10.1.2', '10.1.8'),
  'NW Mountains' = c("6.2.9", "6.2.8", "6.2.15", "6.2.3", "6.2.12", "6.2.7"),
  # Plains
  'Great Plains' = c('9.3.3', '9.3.1', '9.4.1'),
  # GB
  'Southern Great Basin' = c('10.1.5', '10.2.1'),
  # Eastern Intermountain
  'Middle Rockies' = '6.2.10',
  'WY Basin' = '10.1.4',
  'Plateaus' = c('10.1.7', '10.1.6'),
  'Southern Mountains' = c('6.2.14', '6.2.13', '13.1.1'))


df_region2 <- tibble::enframe(r1, name = "region2", value = "NA_L3CODE") %>%
  tidyr::unnest(NA_L3CODE)    

stopifnot(all(df_region2$NA_L3CODE %in% lev3c$NA_L3CODE))

df_region2b <- lev3c %>% 
  right_join(df_region2) %>% 
  group_by(wafwa_NAME, region, region2) %>% 
  summarise(geometry = st_union(geometry), .groups = 'drop') %>% 
  mutate(region2 = factor(region2, levels = names(r1)))

# same naming convention as original 
df_region2c <- df_region2b %>% 
  select(-region) %>% 
  rename(region = region2)

# this code relies on downstream code so won't work as is
# just for development/testing
if(FALSE) {

  r_eco1 <- load_wafwa_ecoregions_raster()
  r_template <- load_template()
  eco4 <- load_wafwa_ecoregions_raster()
  
  # Create a named vector for mapping levels to integers
  lvls <- levels(df_region2b$region2)
  df_region2b$region2_num <- as.integer(df_region2b$region2)  # map to 1, 2, ...
  
  # Rasterize using the numeric codes
  vct_r2b <- terra::vect(df_region2b['region2_num'])
  
  r2_eco <- terra::rasterize(
    vct_r2b,
    r_template,
    field = 'region2_num',
    touches = FALSE
  )
  
  # Set levels (ordered factor)
  levels(r2_eco) <- data.frame(
    value = seq_along(lvls),
    region2 = lvls
  )
  
  r2_eco <- mask(r2_eco, r_template)
  plot(r2_eco, col = RColorBrewer::brewer.pal('Set1', n = length(r1)))
  plot(df_region2b$geometry, add = TRUE, border = "gray", lwd = 2)
  
  
  r_eco <- terra::rasterize(
    terra::vect(lev3c),
    r_template,
    field = 'NA_L3CODE',
    touches = FALSE
  )
  size <- load_cell_size()
  tota_area <- global(size, fun = 'sum', na.rm = TRUE)
  
  tmp <- lev3c[c('NA_L2CODE', 'NA_L3CODE', 'NA_L3NAME', 'region')]
  tmp$geometry <- NULL
  tmp <- distinct(tmp)
  eco_size <- zonal(size, r_eco, fun = 'sum', na.rm = TRUE) %>% 
    filter(!is.na(area)) %>% 
    mutate(area_perc = round(area/as.numeric(tot_area)*100, 1)) %>% 
    filter(!duplicated(.)) %>% 
    left_join(tmp) %>% 
    left_join(df_region2) %>% 
    arrange(region, 'NA_L2CODE', desc(area_perc))
  
  eco_size
  
  # regions that only contain 1% of the total study area
  # (definitely need to be combined)
  too_small <- eco_size %>% 
    as_tibble() %>% 
    filter(area_perc < 1)
  
  r_eco2 <- mask(r_eco, r_template)
  
  lev3d <- lev3c %>% 
    filter(NA_L3CODE %in% eco_size$NA_L3CODE)
  
  plot(eco4)
  plot(lev3d$geometry, add = TRUE, border = "gray", lwd = 2)
  
  
  lev3c_centroids <- lev3d %>%
    st_cast("POLYGON") %>%                          # split multipolygons into separate rows
    mutate(area_ha = as.numeric(st_area(.)) / 1e4) %>%
    group_by(US_L3NAME) %>%
    slice_max(area_ha, n = 1) %>%
    ungroup() %>% 
    st_centroid()
  
  coords <- st_coordinates(lev3c_centroids)
  text(coords[,1], coords[,2], labels = lev3c_biggest$NA_L3CODE, cex = 1.2)
}



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

st_write(df_region2c, 'data_processed/ecoregions/regions_r1.1.gpkg',
         overwrite = TRUE, append = FALSE)
