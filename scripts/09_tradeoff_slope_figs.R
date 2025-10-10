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
    facet_grid(region~ variable, scales = 'free_x', switch = 'x',
               labeller = labeller(variable = driver_labeller())) +
    scale_color_manual(values = cols_GCM2,
                       name = 'GCM (or historical)') +
    scale_linetype_manual(name = 'SEI class',
                          values = c('solid', '11', '41')) +
    scale_shape(name = 'SEI class') +
    theme(strip.placement.x = 'outside',
          strip.text.x = ggtext::element_markdown(),
          axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
    labs(x = NULL,
         y = 'Trade-off slope (degrees)',
         subtitle = paste0(rcp_label(rcp, years, include_parenth = FALSE))#,
         #caption = cap1
         ) +
    scale_y_continuous(breaks = c(-45, 0, 45, 90)) 
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


# * functions ---------------------------------------------------------------

z <- function(x) {
  m <- mean(x)
  sd = sd(x)
  (x - m)/sd
}

calc_intercept = function(x, y, slope) {
  stopifnot(length(slope) == 1)
  xm <- mean(x)
  ym <- mean(y)
  intercept = ym - (slope*xm) # solving for the slope
  intercept
}

# endpoints of a line through a cloud of data given it's slope and intercept
# while keeping the line from extending past a 'box' drawn aroudn the data
calc_segment <- function(xmin, ymin, xmax, ymax, slope, intercept, tol = 1e-12) {
  a <- intercept; b <- slope
  calc_y <- \(x) a + b * x
  calc_x <- \(y) (y - a) / b
  
  # Horizontal line
  if (abs(b) <= tol) {
    if (a >= ymin - tol && a <= ymax + tol) {
      return(tibble(x_start = xmin, y_start = a, x_end = xmax, y_end = a))
    } else {
      return(tibble(x_start = NA_real_, y_start = NA_real_, x_end = NA_real_, y_end = NA_real_))
    }
  }
  
  # Near-vertical guard: treat as vertical segment at x = x0
  if (abs(b) >= 1/tol) {
    # pick x0 from the line using a y inside the box (midpoint is stable)
    y_mid <- 0.5 * (ymin + ymax)
    x0 <- calc_x(y_mid)
    x0 <- min(max(x0, xmin), xmax)  # clip to box in case of numerical drift
    return(tibble(x_start = x0, y_start = ymin, x_end = x0, y_end = ymax))
  }
  
  # Robust interval intersection:
  # parameterize with x = t, y = a + b t; require t in [xmin,xmax] and y in [ymin,ymax]
  t_ymin <- (ymin - a) / b
  t_ymax <- (ymax - a) / b
  t_lo <- max(xmin, min(t_ymin, t_ymax))
  t_hi <- min(xmax, max(t_ymin, t_ymax))
  
  if (t_lo > t_hi + tol) {
    return(tibble(x_start = NA_real_, y_start = NA_real_, x_end = NA_real_, y_end = NA_real_))
  }
  
  x_start <- t_lo
  x_end   <- t_hi
  y_start <- calc_y(x_start)
  y_end   <- calc_y(x_end)
  
  # Final tiny clamp for plotting safety (handles floating-point nicks)
  y_start <- min(max(y_start, ymin), ymax)
  y_end   <- min(max(y_end,   ymin), ymax)
  
  tibble(x_start, y_start, x_end, y_end)
}


# * end functions ---------------------------------------------------------


if (FALSE) {
  # selecting a value to get a 0 slope, but not along a straight line
  fn <- function(a) {
    y <- c(-0.5, -0.44, -0.55, a)
    x <- 1:4
    slope <- calc_slope(x, y)
    abs(slope)
  }
  
  l <- optim(par = list(a = -0.5), fn = fn, method = 'Brent', lower = -1, upper = 1)
  l$par
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
    -0.5, -0.44, -0.55, -0.4633333 # value from optimizer above
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

df3 <- df2 %>% 
  group_by(c3)
s1 <- df3 %>% 
  nest() %>% 
  mutate(slope_SEI = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$SEI)),
         slope_ba = map_dbl(data, .f = \(df) calc_slope(df$graze_num, df$ba)),
         slope_ratio = -1*slope_ba/abs(slope_SEI),
         slope_deg = slope_to_deg(slope_ratio),
         intercept = map_dbl(data, .f = \(df) calc_intercept(x = df$SEI,
                                                             y = df$ba,
                                                             slope = slope_ratio)))


s2 <- df3 %>% 
  summarize(across(c(SEI, ba), .fns = list(min = min, max = max))) %>% 
  right_join(s1, by = 'c3') %>% 
  select(-data)

segments <- s2 %>% 
  group_by(c3) %>% 
  nest() %>% 
  mutate(segments = map(data, \(df) {
    calc_segment(xmin = df$SEI_min, ymin = df$ba_min, xmax = df$SEI_max, 
                          ymax = df$ba_max, slope = df$slope_ratio, 
                          intercept = df$intercept)
  })) %>% 
  select(-data) %>% 
  unnest(cols = segments)

# calculate the 'equation' of the trade-off line, by using
# the calculated slope and forcing it through the mean of x and y


range <- range(c(df2$ba, df2$SEI))
g <- ggplot(df2, aes(SEI, ba, group = c3)) +
  geom_path() +
  geom_point(aes(color = graze)) +
  scale_color_graze() +
  labs(y = 'Burned area (normalized)',
       x = 'Mean SEI (normalized)') +
  coord_cartesian(xlim = range,
                  ylim = range) +
  geom_segment(data = segments,
               aes(x = x_start, xend = x_end,
                   y = y_start, yend = y_end),
               linetype = 2, alpha = 0.5)
g
ggsave('figures/sei/tradeoff/conceptual/SEI_vs_ba_conceptual.png', 
       plot = g,
       width = 6, height = 4, dpi = 300)

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

print(select(s1, -data))
# c3    slope_SEI      slope_ba slope_ratio   slope_deg intercept
# <fct>     <dbl>         <dbl>       <dbl>       <dbl>     <dbl>
# 1 ORA   -4.97e-17 -0.628           1.26e+16 90           1.14e+16
# 2 GOA   -7.70e- 1 -0.665           8.64e- 1 40.8         1.16e- 1
# 3 CSA   -5.13e- 1  0.0000000185   -3.60e- 8 -0.00000206 -9.71e- 1