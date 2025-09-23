# Purpose: Calculate trade-off slopes (of normalized sei and burn probability)
# then plot the trade-off slope against predictor variables, for each GCM

# Author: Martin Holdrege

# Started: September 18, 2025


# params ------------------------------------------------------------------

source("src/params.R")
runv <- opt$runv
years <- opt$years
yr_lab <- opt$yr_lab
vr <- opt$vr
vr_name <- opt$vr_name

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')
source('src/stats_functions.R')
source('src/SEI_functions.R')
options(readr.show_col_types = FALSE)
theme_set(theme_custom1())
# read in data ------------------------------------------------------------

# file created in 08_fire_area_figs.R
c3eco1 <- read_csv(
  paste0('data_processed/raster_means/',  runv, '_', vr, "_", years, 
         '_sei-mean_ba_by-GCM-region-c3.csv'))

# means of drivers of fire
# created in "scripts/06_summarize_fire_drivers.R"
driver1 <- read_csv(paste0('data_processed/raster_means/', run, '_', vr,
                           '_fire-driver-means_by-ecoregion-c3.csv'))

# *gridded fire & sei -----------------------------------------------------
# used for normalizing the tabular data (need mean and sd)

# fire probability by grazing level for current conditions
# from "scripts/03_interpolate.R"
fire_files <- list.files(
  file.path("data_processed/interpolated_rasters/fire/", v_interp),
  pattern = paste0(run, '.*Current'), full.names = TRUE) 

stopifnot(length(fire_files) == 4) # 4 grazing levels
r_fire1 <- rast(fire_files)

# sei
# from "scripts/05_interpolated_summarize_sei_scd-adj.R"
r_sei1 <- rast(file.path("data_processed/interpolated_rasters/", v_interp,
          paste0(runv, yr_lab, "_q-sei_scd-adj_summary.tif")))

# summarize gridded product -----------------------------------------------

r_sei2 <- r_sei1[[str_subset(names(r_sei1), 'SEI_Current')]]
stopifnot(nlyr(r_sei2) == 4)

# the directly interpolated rasters have a few pixels outside the ecoregions
r_fire2 <- mask(r_fire1, r_sei2[[1]])

fire_stats <- global_smry(r_fire2)
sei_stats <- global_smry(r_sei2)


# calculate slopes --------------------------------------------------------

normalize <- function(x, stats) {
  stopifnot(c('mean', 'sd') %in% names(stats))
  (x - stats[['mean']])/stats[['sd']]
}

c3eco2 <- c3eco1 %>% 
  df_factor() %>% 
  mutate(SEI_norm = normalize(SEI_mean, stats = sei_stats),
         # convert % to proportion
         ba_norm = normalize(expected_ba_perc, stats = fire_stats))
# hist(c3eco2$ba_norm)
# hist(c3eco2$SEI_norm)
# 
# summary(c3eco2$SEI_norm)
# with(filter(c3eco2, RCP == 'Current'), 
#      weighted.mean(SEI_norm, w = area))

# calculating slopes (s)
slopes1 <- c3eco2 %>% 
  select(-rcp_year_c3_gcm) %>% 
  group_by(run, RCP, years, rcp_year, GCM, region, c3) %>% 
  nest() %>% 
  mutate(slope_deg = map_dbl(data, .f = \(df) {
    calc_slope_deg(x = df$SEI_norm, y = df$ba_norm)
  })) %>% 
  select(-data) %>% 
  ungroup()


# combine data ------------------------------------------------------------
slopes2 <- driver1 %>% 
  # climate variables don't have separate values per grazing level
  # and graze is NA, otherwise just choosing moderate b/ it's the reference
  # grazing level, could also have taken an average across grazing trms
  filter(is.na(graze) | graze == 'Moderate') %>% 
  select(-graze) %>% 
  right_join(slopes1, by = join_by(years, RCP, GCM, region, c3)) %>% 
  df_factor() %>% 
  mutate(GCM = factor(GCM, levels = c('Current', names(cols_GCM1)),
                      labels = c('Historical', names(cols_GCM1))),
         variable = driver2factor(variable, include_sagebrush = TRUE))

# figures -----------------------------------------------------------------

rcps <- c('RCP45', 'RCP85')


for(rcp in rcps) {
  df <- slopes2 %>% 
    filter(RCP %in% c('Current', rcp))
  g <- ggplot(df, aes(mean, slope_deg)) +
    geom_smooth(aes(linetype = c3), se = FALSE,
                alpha = 0.7) +
    geom_point(aes(shape = c3, color = GCM)) +
    facet_grid(region~ variable, scales = 'free_x', switch = 'x') +
    scale_color_manual(values = cols_GCM2,
                       name = 'GCM (or historical)') +
    scale_linetype_manual(name = 'SEI class',
                          values = c('solid', '11', '41')) +
    scale_shape(name = 'SEI class') +
    theme(strip.placement.x = 'outside') +
    labs(x = NULL,
         y = 'Trade-off slope (degrees)',
         subtitle = paste0(rcp_label(rcp, years, include_parenth = FALSE)))
  ggsave(
    paste0('figures/sei/tradeoff/slope_by-region-c3_',
           vr, '_', rcp, '_', years, '_', run, '.png'),
    plot = g, width = 11, height = 9, dpi = 600
    )
}



