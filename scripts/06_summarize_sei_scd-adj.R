# Purpose:
# calculate five number summaries of Qs and SEI (remote sensing adjusted) values
# for creating boxplots and related figures for each ecoregion

# Author: Martin Holdrege
# Started: May 12, 2025

# params ------------------------------------------------------------------

source("src/params.R")
test_run <- FALSE
runv <- paste0(run, v_interp)
fire_breaks <- c(seq(0, 2.5, by = 0.25), 3, 4) # intervals to cut fire probability into

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')
source("src/SEI_functions.R")

# read in data ------------------------------------------------------------

 # file create in 05_interpolated_summarize_sei_scd-adj.R
r_qsei1 <- rast(file.path("data_processed/interpolated_rasters/", v_interp,
                     paste0(runv, "_q-sei_scd-adj_summary.tif")))

eco1 <- load_wafwa_ecoregions(wafwa_only = FALSE, total_region = TRUE)
eco2 <- vect(eco1['ecoregion'])

r_eco1 <- load_wafwa_ecoregions_raster(wafwa_only = FALSE)

size <- load_cell_size()

# * fire ------------------------------------------------------------------
# this is for summarizing SEI at different bins of fire probability
# created in "scripts/05_interpolated_maps_fire.R"
r_fire1 <- rast(file.path("data_processed/interpolated_rasters", 
          paste0(run, "_fire-prob_future_summary_across_GCMs.tif")))


# * setup -----------------------------------------------------------------


if(test_run) {
  r_qsei1 <- downsample(r_qsei1)
  r_fire1 <- downsample(r_fire1)
  r_eco1 <- downsample(r_eco1)
  size <- downsample(size)
}


# names of layers --------------------------------------------------

into <- c("group", "type", "RCP", "years", "graze", 'summary')

info1 <- create_rast_info(r_qsei1, into = into) 

# calculate means by ecoregions -------------------------------------------

# taking  unweighted summary stats is ok, because we're using an equal area projection
sixnum <- function(x, na.rm = TRUE) {
  five <- fivenum(x, na.rm = na.rm)
  c(five, mean(x, na.rm = na.rm))
}
# 5 number summary by ecoregion
eco_smry_l1 <- map(info1$id, function(id) {
  tmp <- terra::extract(r_qsei1[[id]], eco2, sixnum,  na.rm = TRUE)
  tmp <- as_tibble(tmp)
  colnames(tmp) <- c('ID', 'min', 'lower', 'middle', 'upper', 'max', 'mean')
  tmp$'id' <- id
  tmp$region <- values(eco2)$ecoregion
  tmp$ID <- NULL
  tmp
})

eco_smry1 <- eco_smry_l1 %>% 
  bind_rows() %>% 
  left_join(info1, by = 'id')

# calculating % CSA & GOA ------------------------------------------------------

info_sei1 <- info1[info1$type == 'SEI', ]
r_sei <- r_qsei1[[info_sei1$id]]

# % core
region_area0 <- terra::extract(size, eco1, sum,  na.rm = TRUE)
region_area <- region_area0 %>% 
  mutate(region = as.character(eco1$ecoregion)) %>% 
  select(-ID) %>% 
  rename(area_region = area)

sei_core1 <- terra::extract(r_sei, eco1, percent_csa,  na.rm = TRUE)
sei_core2 <- pivot_longer_extracted(sei_core1, as.character(eco1$ecoregion), 
                                    values_to = 'percent_csa') 
# % growth opportunity area
sei_goa1 <- terra::extract(r_sei, eco1, percent_goa,  na.rm = TRUE)
sei_goa2 <- pivot_longer_extracted(sei_goa1, as.character(eco1$ecoregion), 
                                    values_to = 'percent_goa') 

sei_pcent <- left_join(sei_core2, sei_goa2, by = c('region', 'id')) %>% 
  left_join(region_area, by = 'region') %>% 
  left_join(info1, by = 'id') %>% 
  select(-run2)


sei_pcent_long <- sei_pcent %>% 
  mutate(percent_ora = 100 - percent_goa - percent_csa) %>% 
  pivot_longer(cols = matches('percent_'),
               names_to = 'c3',
               values_to = 'c3_percent') %>% 
  mutate(c3_area = c3_percent/100*area_region,
         c3 = str_to_upper(str_replace(c3, 'percent_', ''))) %>% 
  select(-group, -type, -id)

# calculating C9 transition -----------------------------------------------
# within a grazing level, and for a given scenario the area falling into
# each of 9 possible SEI class changes


r_c3 <- sei2c3(r_sei)

info_sei2 <- info_sei1 %>% 
  filter(RCP == 'Current') %>% 
  select(-RCP, -years, -run2, -summary, -type) %>% 
  right_join(filter(info_sei1, RCP != 'Current'),
             by = join_by(run, group, graze),
             suffix = c('_cur', '_fut')) 

r_c9a <- c3toc9(current = r_c3[[info_sei2$id_cur]],
                future = r_c3[[info_sei2$id_fut]])

r_c9b <- c(r_eco1, size, r_c9a)

# * area of c9 transition -------------------------------------------------

c9_area1 <- as.data.frame(r_c9b) %>% 
  rename(region = ecoregion) %>% 
  mutate(cell_num = 1:n()) %>% 
  pivot_longer(matches('^fire.*eind'),
               values_to = 'c9',
               names_to = "id") %>% 
  mutate(c9 = c9_factor(c9)) %>% 
  group_by(id, region, c9) %>% 
  summarise(area = sum(area),
            .groups = 'drop_last') %>% 
  mutate(area_perc = area/sum(area)*100)


# adding in entire study region
c9_area2 <- c9_area1 %>% 
  group_by(id, c9) %>% 
  summarise(area = sum(area), .groups = 'drop_last') %>% 
  mutate(area_perc = area/sum(area)*100,
           region = levels(eco1$ecoregion)[1]) %>% 
  bind_rows(c9_area1) %>% 
  ungroup() %>% 
  left_join(info1, by = 'id') %>% 
  mutate(type = 'c9') %>% 
  select(-run2)

# check
test <- c9_area2 %>% 
  group_by(id, region) %>% 
  summarize(total = sum(area_perc)) %>% 
  pull(total)

stopifnot(abs(test - 100) < 0.0001)

# mean SEI by fire probability bin ----------------------------------------

summary4bins <- 'median' # summary across GCMs (pixelwise) used
info_sei <- info1 %>% 
  filter(summary == summary4bins, group == 'SEI') %>% 
  select(-run, -run2, -group, -type)

info_fire1 <- create_rast_info(r_fire1,
                               into = c("type", "RCP", "years", "graze", 'summary'))
info_seifire <- info_fire1 %>% 
  select(-type) %>% 
  right_join(info_sei, by = c('RCP', 'years', 'graze', 'summary'),
             suffix = c('_fire', '_sei'))

r_fire2 <- r_fire1[[info_seifire$id_fire]]
r_sei1 <- r_qsei1[[info_seifire$id_sei]]

# * binning fire ----------------------------------------------------------
r_comb <- c(r_eco1, r_sei1, r_fire2)

df_seifire1 <- as.data.frame(r_comb)

stopifnot(max(fire_breaks) > max(unlist(minmax(r_fire2))))

df_seifire2 <- df_seifire1 %>% 
  rename(region = ecoregion) %>% 
  mutate(cell_num = 1:n()) %>% 
  pivot_longer(- c(region, cell_num), names_to = 'id') %>% 
  left_join(bind_rows(info_sei, info_fire1[names(info_sei)]), by = 'id') %>% 
  mutate(group = str_extract(id, 'SEI|fire-prob')) %>% 
  select(-id) %>% 
  pivot_wider(values_from = 'value',names_from = "group") %>% 
  mutate(fire_bin = cut(`fire-prob`, fire_breaks, right = FALSE)) %>% 
  # some pixels outside of the regions (from the fire raster)
  # discarding these
  filter(!is.na(region))

summarize_sei <- function(df) {
  summarize(df,
            n_pixels = n(),
            percent_csa = percent_csa(SEI),
            percent_goa = percent_goa(SEI),
            SEI = mean(SEI),
            .groups = 'drop')
}

df_seifire_regions <- df_seifire2 %>% 
  group_by(region, RCP, years, graze, summary, fire_bin) %>% 
  summarize_sei()

df_seifire_entire <- df_seifire2 %>% 
  select(-region) %>% 
  group_by(RCP, years, graze, summary, fire_bin) %>% 
  summarize_sei() %>% 
  mutate(region = eco1$ecoregion[1]) # entire study area

df_seifire3 <- bind_rows(df_seifire_entire, df_seifire_regions)

# save output -------------------------------------------------------------

prefix <- if(test_run) 'test' else runv

write_csv(eco_smry1, paste0('data_processed/raster_means/', prefix, 
                            '_q-sei_scd-adj_summaries_by-ecoregion.csv'))

write_csv(sei_pcent , paste0('data_processed/raster_means/', prefix, 
                            '_sei-class-pcent_scd-adj_summaries_by-ecoregion.csv'))

write_csv(sei_pcent_long, paste0('data_processed/raster_means/', prefix, 
                             '_sei-class-pcent-long_scd-adj_summaries_by-ecoregion.csv'))

write_csv(df_seifire3, paste0('data_processed/raster_means/', prefix, 
                              '_sei-by-fire-bin_scd-adj_summaries_by-ecoregion.csv'))

write_csv(c9_area2, paste0('data_processed/raster_means/', prefix, 
                              '_c9-area_scd-adj_summaries_by-ecoregion.csv'))
