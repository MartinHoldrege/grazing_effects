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
  mutate(graze_num = as.numeric(graze)) %>% 
  group_by(run, RCP, years, rcp_year, GCM, region, c3) %>% 
  nest() %>% 
  mutate(
    slope_SEI = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$SEI_norm)),
    slope_ba = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$ba_norm)),
    slope_deg0 = slope_to_deg(-1*slope_ba/abs(slope_SEI)),
    range_SEI = map_dbl(data, \(df) max(df$SEI_norm) - min (df$SEI_norm)),
    range_ba = map_dbl(data, \(df) max(df$ba_norm) - min(df$ba_norm))
  ) %>% 
  # select(-data) %>% 
  # ungroup() %>% 
  group_by(.row = 1:n()) %>% 
  # maximum absolute range in sei or ba in units of SD
  mutate(range_max = max(c(abs(range_SEI), abs(range_ba))),
         slope_max = max(c(abs(slope_SEI), abs(slope_ba)))) %>% 
  ungroup() %>% 
  select(-.row)


delta_cutoff <- 0.05
hist(slopes1$slope_max, breaks = 20) # 'low' point in the histogram at about 0.05
slopes2 <- slopes1 %>% 
  mutate(
    slope_deg = ifelse(slope_max < delta_cutoff, NA, slope_deg0) 
         
  )
hist(slopes1$slope_deg0)
hist(slopes2$slope_deg)
# combine data ------------------------------------------------------------
slopes3 <- driver1 %>% 
  # climate variables don't have separate values per grazing level
  # and graze is NA, otherwise just choosing moderate b/ it's the reference
  # grazing level, could also have taken an average across grazing trms
  filter(is.na(graze) | graze == 'Moderate') %>% 
  select(-graze) %>% 
  right_join(slopes2, by = join_by(years, RCP, GCM, region, c3)) %>% 
  df_factor() %>% 
  mutate(GCM = factor(GCM, levels = c('Current', names(cols_GCM1)),
                      labels = c('Historical', names(cols_GCM1))),
         variable = driver2factor(variable, include_sagebrush = TRUE))


# figures -----------------------------------------------------------------

rcps <- c('RCP45', 'RCP85')

cap1 <- paste0('slope = -1*(delta burned area)/abs(delta SEI),\n',
               'where the delta values are change units (standard deviations) ',
               'per grazing intensity level increase.\n',
               'Slope is NA if neither delta was > |', delta_cutoff, '|')
for(rcp in rcps) {
  df <- slopes3 %>% 
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
         subtitle = paste0(rcp_label(rcp, years, include_parenth = FALSE)),
         caption = cap1)
  g
  ggsave(
    paste0('figures/sei/tradeoff/slope_by-region-c3_',
           vr, '_', rcp, '_', years, '_', run, '.png'),
    plot = g, width = 11, height = 9, dpi = 600
    )
}



# conceptual figure -------------------------------------------------------
# conceptual figure showing how the trade-off slope is calculated,
# using artificial data

z <- function(x) {
  m <- mean(x)
  sd = sd(x)
  (x - m)/sd
}
df <- tibble(
  graze = c(
    'Light', 'Moderate', 'Heavy', 'Very Heavy',
    'Light', 'Moderate', 'Heavy', 'Very Heavy',
    'Light', 'Moderate', 'Heavy', 'Very Heavy'
  ),
  SEI = c(
    c(-0.3, -0.1, -0.1, -0.3) ,
    rev(c(-0.3, 0, 0.3, 0.6)),
    rev(c(0.2, 0.4, 0.6, 0.8))
  ),
  ba = c(
    1, 0.7, 0.3, 0,
    0.7, 0.2, -0.1, -0.4,
    -0.5, -0.5, -0.5, -0.5
  ),
  c3 = c(
    rep('ORA', 4),
    rep('GOA', 4),
    rep('CSA', 4)
  )
)

df2 <- df %>% 
  df_factor() %>% 
  mutate(ba = z(ba),
         SEI = z(SEI),
         graze_num = as.numeric(graze)) 

range <- range(c(df2$ba, df2$SEI))
g <- ggplot(df2, aes(SEI, ba, group = c3)) +
  geom_path() +
  geom_point(aes(color = graze)) +
  scale_color_graze() +
  labs(y = 'Burned area (normalized)',
       x = 'Mean SEI (normalized)') +
  coord_cartesian(xlim = range,
                  ylim = range)

ggsave('figures/sei/tradeoff/conceptual/SEI_vs_ba_conceptual.png', 
       plot = g,
       width = 6, height = 4, dpi = 300)

s1 <- df2 %>% 
  group_by(c3) %>% 
  nest() %>% 
  mutate(slope_SEI = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$SEI)),
         slope_ba = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$ba)),
         slope_deg = slope_to_deg(-1*slope_ba/abs(slope_SEI))) %>% 
  select(-data)
  
df_l <- df2 %>% 
  group_by(c3) %>% 
  group_split()

plot_graze <- function(df, y, y_title) {
  ggplot(df, aes(graze_num, .data[[y]])) +
    geom_smooth(method = 'lm', se = FALSE, color = 'gray') +
    geom_point(aes(color = graze)) +
    coord_cartesian(ylim = range) +
    theme(axis.text = element_blank()) +
    labs(x = 'Grazing',
         y = y_title) +
    scale_color_graze() +
    theme(legend.position = 'none')
}


map(df_l, \(df) {
  y <- 'SEI'
  g <- plot_graze(df, y = y, y_title = y)
  file <- paste0('figures/sei/tradeoff/conceptual/',
                 y, '_vs_graze_', unique(df$c3), '.png')
  ggsave(file, 
         plot = g,
         width = 1.2, height = 1.2, dpi = 300)
})

map(df_l, \(df) {
  y <- 'ba'
  g <- plot_graze(df, y = y, y_title = 'Burned Area')
  file <- paste0('figures/sei/tradeoff/conceptual/',
                 y, '_vs_graze_', unique(df$c3), '.png')
  ggsave(file, 
         plot = g,
         width = 1.2, height = 1.2, dpi = 300)
})
