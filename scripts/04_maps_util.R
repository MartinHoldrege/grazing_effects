# Purpose: load utilization rasters, and combine
# and calculate delta utilization


# parameters --------------------------------------------------------------

source("src/params.R")
summary <- 'median'

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source("src/general_functions.R")
source("src/mapping_functions.R")
source("src/fig_params.R")
source("src/fig_functions.R")

# read in data ------------------------------------------------------------

folder <- file.path("data_processed/interpolated_rasters/util", v_interp)
regex <- paste0(run, '.*', summary, '.tif$')
paths <- list.files(folder, regex, full.names = TRUE)
r_util1 <- rast(paths)


# raster info -------------------------------------------------------------

into <- c("PFT", "type", "RCP", "years", "graze", "summary")
info_util1 <- create_rast_info(r_util1, into = into)


# calculate differences ---------------------------------------------------

type_absolute <- 'util'
type_diff <- 'util-diff-cref'
type_diff_perc <- paste0(type_diff, '-perc')
match_vars <- c('run', 'PFT', 'type', 'graze', 'summary')

r_util_diff1 <- rast_diff(r = r_util1, info = info_util1,
                          type_absolute = type_absolute,
                          type_diff = type_diff,
                          match_vars = match_vars,
                          include_percent = TRUE)

# 20 panel maps -----------------------------------------------------------

info_20panel <- create_rast_info(r_util_diff1, into = into) %>% 
  bind_rows(info_util1) %>% 
  filter(years != '2030-2060' & years != '2031-2060') 

args <- list(
  PFT = c('Aherb', 'Pherb', 'Sagebrush'),
  titles = c(Aherb = 'Annual herbacious utilization',
             Pherb = 'Perennial herbacious utilization',
             Sagebrush = 'Sagebrush utilization')
)

r <- c(r_util1, r_util_diff1)


for (pft in args$PFT) {
  info_tmp <- info_20panel %>% 
    filter(PFT == pft)
  
  title <- args$titles[pft]
  
  # absolute change
  g <- plot_map_20panel(
    r = r,
    info = info_tmp[info_tmp$type != type_diff_perc, ],
    type_absolute = type_absolute,
    type_diff =  type_diff,
    title = title,
    name4absolute = "Utilization",
    legend_title_absolute = lab_util0,
    legend_title_diff = lab_util1
  )
  
  filename <- paste0('figures/util/maps/20panel_', pft, '_util_', run, ".png")
  png20panel(filename)
  print(g)
  dev.off()
  
  # % change 
  lims_diff <- range_raster(r_util_diff1[[info_tmp$id[str_detect(info_tmp$type, '-perc$')]]],
                            absolute = TRUE)
  
  if(lims_diff[1] < -100) {
    lims_diff[1] <- -100
  } 
  
  if(lims_diff[2] > 300) {
    lims_diff[2] <- 300
  }
  g <- plot_map_20panel(
    r = r,
    info = info_tmp[info_tmp$type != type_diff, ],
    type_absolute = type_absolute,
    type_diff = type_diff_perc,
    title = title,
    name4absolute = "Utilization",
    legend_title_absolute = lab_util0,
    legend_title_diff = expression(~Delta * " Utilization (%)"),
    lims_diff = lims_diff,
    midpoint_diff = 0
  )
  
  filename <- paste0('figures/util/maps/20panel_', pft, '_util-perc_', run, ".png")
  png20panel(filename)
  print(g)
  dev.off()
  
}


