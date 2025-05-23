# Purpose: Compare stepwat and RAP biomass, at the 200 simulation sites 
# so that some sort of quantile matching relationship can determined
# so that stepwat biomass is translated to biomass levels used for 
# fitting the Holdrege et al. 2024 fire equation

# Author: Martin Holdrege

# Oct 3, 2024

# params ------------------------------------------------------------------

version_interp <- 'v3' # interpolation version (see 03_interpolate.R)
# the first run in the vector is the 'main' one, ( and compared to the
# no fire run in some figures)
runs <- c("fire1_eind1_c4grass1_co20_2311")

# matching quality cutoff, this decides
# the cells from which both stepwat and rap data will be pulled from
qual_cutoff <- 0.5 # 1 # decided to use 0.5, (results are similar to using 1)
RCP <- 'Current'
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
source("src/general_functions.R")
source("src/mapping_functions.R")
# read in data ------------------------------------------------------------

# data from the Holdrege et al. 2024 publication
rap_df1 <- read_csv("../cheatgrass_fire/data_processed/data_publication/Wildfire_Climate_Biomass_Data.csv")

cell_nums <- rast("../cheatgrass_fire/data_processed/data_publication/cell_nums.tif")

sites1 <- read_csv("data_raw/site_locations.csv")

# * interpolated stepwat --------------------------------------------------
PFTs <- c('Aherb', 'Pherb')

suffix <- paste0(version_interp, ".*", RCP, ",*.tif")
pattern <- paste0("(", paste0(runs, collapse = '|'), ")", suffix)

paths <- list.files("data_processed/interpolated_rasters/biomass/",
                    pattern = pattern,
                    full.names = TRUE)
basename(paths)
r_sw1 <- rast(paths)

# for filter which data to use
r_qual <- rast(paste0("data_processed/interpolation_quality/matching_quality",
                 version_interp, ".tif"))

# for later reproducibility (i.e., know which sites were used, where)
# gives the site number, of the site that was matched to the given pixel
r_site <- rast(paste0("data_processed/interpolation_data/interp_locations_200sites_",
                     version_interp, ".tif"))

# vectors ----------------------------------------------------------------

pft_lookup <- c('pfgAGB' = 'Pherb','afgAGB' = 'Aherb')

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
  mutate(PFT = pft_lookup[PFT],
         dataset = 'site level') 


# pulling together interpolation level data -------------------------------

# identify stepwat grid-cells w/ good matching, and that are also cells
# included in the fire manuscript
cell_nums2 <- extend(cell_nums, r_qual)
mask0 <- r_qual < qual_cutoff & !is.na(cell_nums2)
mask <- classify(mask0, rcl = cbind(FALSE, NA))

r_site2 <- mask(r_site, mask)
names(r_site2) <- 'site'
r_sw2 <- c(mask(r_sw1, mask), r_site2)


cell_nums_m <- mask(cell_nums2, mask)
cell_nums_mv <- as.numeric(cell_nums_m[!is.na(cell_nums_m)]) # vector

rap_df_m <- rap_df2 %>% 
  filter(cell_num %in% cell_nums_mv) %>% 
  pivot_longer(cols = c("pfgAGB", 'afgAGB'),
               names_to = 'PFT',
               values_to = 'biomass') %>% 
  mutate(PFT = pft_lookup[PFT],
         dataset = 'interpolated')

rap_comb <- bind_rows(rap_df_m, rap_sub2) # all the rap data

df_sw1 <- as_tibble(r_sw2) 
df_sw1$cell_num <- terra::cells(r_sw2)

info1 <- create_rast_info(r_sw1)

df_sw2 <- df_sw1 %>% 
  pivot_longer(cols = -c(cell_num, site),
               names_to = 'id',
               values_to = 'biomass') %>% 
  left_join(info1, by = 'id') %>% 
  mutate(fire = str_extract(run, 'fire\\d'))

tmp1 <- df_sw2 %>% 
  select(run, fire, biomass, PFT, graze, cell_num, site) %>% 
  mutate(dataset = 'interpolated') 

sw_comb <- tmp1%>% 
  # remove the (interpolation) version number
  # so run names are the same between site level and interpolated
  mutate(run = str_replace(run, 'v\\d$', ""))

# figure of inerpolation area ---------------------------------------------

study <- cell_nums2
study[study > 0] <- 1
mask_tmp <- mask
mask_tmp[is.na(mask_tmp)] <- 0
study2 <- study + mask_tmp

# percent of fire ecology study area covered with strict interpolation
pcent <- round(sum(values(mask), na.rm = TRUE)/sum(values(study), na.rm = TRUE)*100)

cls <- data.frame(id=1:2, name=c("Study area used in Fire Ecology",
                                  'Subset of area used for quantile mapping'))
levels(study2) <- cls

bb <- st_bbox(terra::trim(study2))
m <- plot_map2(terra::trim(study2), expand_bbox = FALSE) +
  scale_fill_manual(values = c('grey', 'blue'), na.value = 'transparent',
                    name = NULL, na.translate = FALSE) +
  geom_sf(data = sites2, size = 1) +
  ggplot2::coord_sf(
    xlim = bb[c("xmin", "xmax")],
    ylim = bb[c("ymin", "ymax")],
    crs = crs_scd) +
  theme(legend.position = 'bottom',
        legend.direction = 'vertical') + 
  labs(
    subtitle = paste('Area used for interpolation for quantile mapping\n',
                     'with STEPWAT2 sites also shown'),
    caption = paste0(
      pcent, "% of the Holdrege et al (Fire Ecology) study area is covered by\n",
      " stepwat interpolated data when a strict (", qual_cutoff, ") matching quality",
      " cutoff is used."
    )
    )

jpeg(paste0("figures/bio_matching/interpolation_area_qm",  qual_cutoff,".jpg"),
     res = 800, width = 6, height = 6.5, units = 'in') 
print(m)
dev.off()

# saving output -----------------------------------------------------------

out <- list(
  runs = runs,
  sw_comb = sw_comb,
  rap_comb = rap_comb,
  qual_cutoff = qual_cutoff
)

saveRDS(out, 'data_processed/temp_rds/rap_sw_matching.rds')

# saving these csv's for reproducibility
if(qual_cutoff == 0.5) {
  sw <- sw_comb %>% 
    filter(dataset == 'interpolated')
  
  rap <- rap_comb %>% 
    filter(dataset == 'interpolated') %>% 
    select(-site)
  
  write_csv(sw, paste0('data_processed/qm/stepwat_for_qm_', qual_cutoff, 
                       'match_v2.csv'))
  
  write_csv(rap, paste0('data_processed/qm/rap_for_qm_', qual_cutoff, 
                       'match_v2.csv'))
}