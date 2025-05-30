# Purpose: Create maps of Q and SEI, these are values where delta's from 
# stepwat were applied to SCD cover layers (aggregated to 1km)

# Author: Martin Holdrege

# Started: May 8, 2025

# parameters --------------------------------------------------------------

source('src/params.R')
runv <- paste0(run, v)
groups <- c('Sagebrush', 'Pherb', 'Aherb', 'SEI')
test_run <- FALSE

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source("src/general_functions.R")
source("src/fig_params.R")
source("src/mapping_functions.R")
source('src/fig_functions.R')

# read in data ------------------------------------------------------------

# actual (non difference values)
# cover
r_cov1 <- rast(
  file.path("data_processed/interpolated_rasters", v,
            paste0(runv, "_cover_scd-adj_summary.tif"))
)
# q and sei
r_qsei1 <- rast(
  file.path("data_processed/interpolated_rasters", v,
            paste0(runv, "_q-sei_scd-adj_summary.tif"))
)

# difference from current conditions
# cover
r_cov_d1 <- rast(
  file.path("data_processed/interpolated_rasters", v,
            paste0(runv, "_cover-rdiff-cref_scd-adj_summary.tif"))
)
# q and sei
r_qsei_d1 <- rast(
  file.path("data_processed/interpolated_rasters", v,
            paste0(runv, "_q-sei-rdiff-cref_scd-adj_summary.tif"))
)


# prepare data ------------------------------------------------------------

if(test_run) {
  size <- 100
  r_cov1 <- spatSample(r_cov1, size = size, method = 'regular',
                      as.raster = TRUE)
  r_qsei1 <- spatSample(r_qsei1, size = size, method = 'regular',
                         as.raster = TRUE)
  r_cov_d1 <- spatSample(r_cov_d1, size = size, method = 'regular',
                       as.raster = TRUE)
  r_qsei_d1 <- spatSample(r_qsei_d1, size = size, method = 'regular',
                       as.raster = TRUE)
}


r_comb1 <- c(r_cov1, r_cov_d1, r_qsei1, r_qsei_d1)
into = c("group", "type", "RCP", "years", "graze", "summary")

info1 <- create_rast_info(r_comb1, into = into) %>% 
  filter(years != '2030-2060' & years != '2031-2060',
         summary == 'median') 


# 20 panel map (absolute delta) -------------------------------------------


args <- list(
  groups = c('Sagebrush', 'Pherb', 'Aherb', 
             'Sagebrush', 'Pherb', 'Aherb', 
             'SEI'),
  type_absolute = c('cover', 'cover', 'cover',
                    'Q', 'Q', 'Q',
                    'SEI'),
  full_names = c('Sagebrush' = 'Sagebrush',
                 'Pherb' = 'Perennial herbacious',
                 'Aherb' = 'Annual herbacious',
                 'SEI' = 'Sagebrush ecological integrity'),
  type_suffix = c('cover' = 'cover \n(modeled changes combined with remote sensing estimates)',
                  'Q' = 'quality score \n(modeled changes combined with remote sensing estimates)',
                  'SEI' = '\n(modeled changes combined with remote sensing estimates)'),
  palette_absolute = list(
    'cover' = cols_map_bio(10),
    'Q' = cols_map_bio(10),
    'SEI' = cols_seicont
  ),
  values_absolute = list(
    'cover' = NULL,
    'Q' = NULL,
    'SEI' = vals_seicont
  ),
  legend_title_absolute = c('cover' = 'cover (%)',  'Q' = 'Q', 'SEI' = 'SEI'),
  legend_title_diff = c(cover = "\u0394 cover (percentage points)",
                        Q = '\u0394 Q',
                        SEI = '\u0394 SEI')

)

args$type_diff <- unique(args$type_absolute) %>% 
  str_replace('cover', 'cov') %>% 
  paste0(., "-rdiff-cref") %>% 
  setNames(unique(args$type_absolute))

map2(args$groups, args$type_absolute, function(group, type) {

  info_tmp <- info1 %>% 
    filter(.data$group == !!group,
           .data$type == !!type | .data$type == args$type_diff[[!!type]])
  
  title <- paste(args$full_names[[group]], args$type_suffix[[type]])
  
  g <- plot_map_20panel(
    r = r_comb1,
    info = info_tmp,
    type_absolute = type,
    type_diff =  args$type_diff[[type]],
    title = title,
    name4absolute = type,
    palette_absolute =  args$palette_absolute[[type]],
    palette_diff = cols_map_bio_d2,
    values_absolute =  args$values_absolute[[type]],
    legend_title_absolute = args$legend_title_absolute[[type]],
    legend_title_diff = args$legend_title_diff[[type]]
  )
  
  prefix <- if(group == "SEI") {
    group
  } else {
    paste0(group, "-", type)
  }
  
  filename <- paste0('figures/sei/maps/20panel_', prefix, '_scd-adj_', runv, ".png")
  png20panel(filename)
  print(g)
  dev.off()
  
})

# SEI class and change maps -----------------------------------------------


c3_graze <- c('M' = 'Moderate', 'VH' = 'Very Heavy')
c3_rcp <- 'RCP45' # future scenario shown on map


c3_info <- info1 %>% 
  filter(type == 'SEI', graze %in% c3_graze,
         RCP %in% c('Current', c3_rcp)) %>% 
  select(-run2)

c3_info2 <- c3_info %>% 
  filter(RCP == 'Current') %>% 
  select(-RCP, -years) %>% 
  left_join(c3_info[c3_info$RCP != 'Current', ],
            by = join_by(run, group, type, graze, summary),
            suffix = c('_cur', '_fut'))

r_c3 <- sei2c3(r_comb1[[c3_info$id]])

r_c9 <- c3toc9(current = r_c3[[c3_info2$id_cur]],
               future = r_c3[[c3_info2$id_fut]])



g <- plot_4panel_c3c9(r_c3 = r_c3, r_c9 = r_c9,
                 info_c3 = c3_info[c3_info$RCP == 'Current', ],
                 info_c9 = c3_info[c3_info$RCP != 'Current', ]
                 )
filename <- paste0('figures/sei/maps/4panel_c3c9_RCP45_MVH_scd-adj_', runv, ".png")
png(filename, res = 600, height = 3*1.8 + 1, width = 3*1.8 - 0.3, units = 'in')
g
dev.off()
