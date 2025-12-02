# Purpose: create maps showing what's driving the 
# spatial pattern in wildfire 


# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source('src/general_functions.R')
source('src/fig_params.R')
source('src/mapping_functions.R')
source('src/fig_functions.R')
theme_set(theme_custom1())
# params ------------------------------------------------------------------

source('src/params.R')
runv <- opt$runv
ref_graze <- opt$ref_graze
v_interp <- opt$v_interp
test_run <- opt$test_run

rcp_yr <- 'Current_Current'
# read in data ------------------------------------------------------------
# this rasters created in 03_interpolation.R
# but using inputs from 06_fire_attribion_sref.R (i.e.,
# unlike the normal approach,
# the 03_ needs to re-run after the 06_ to create these files)

files <- list.files(
  file.path("data_processed/interpolated_rasters/fire_attrib_sref", v_interp),
  pattern = paste0(runv, ".*", rcp_yr, '.*', ref_graze),
  full.names = TRUE)

# right now the code is setup for running for reference
# conditions (historical climate, ref graze) only
stopifnot(length(files) == 5)

r1 <- rast(files)

name_suffix <- paste0('_', rcp_yr, '_', ref_graze, '_Current.tif')

# read in biomass
pfts <- c('Aherb', 'Pherb')
files_bio <- paste0("data_processed/interpolated_rasters/biomass/", v_interp,
                    '/', runv, '_', pfts,  '_biomass', name_suffix)

r_bio <- rast(files_bio)
names(r_bio) <- pfts

# read in climate
clim_vars <- c('MAT', 'MAP', 'psp')

files_clim <- paste0("data_processed/interpolated_rasters/climate/", v_interp,
                    '/', clim_vars, '_climate_', rcp_yr, '_Current_20250228_', 
                    v_interp, '.tif')

r_clim <- rast(files_clim)
names(r_clim) <- str_to_upper(clim_vars)

r_driv1 <- c(r_clim, r_bio)

if(test_run) {
  r1 <- downsample(r1)
  r_driv1 <- downsample(r_driv1)
}

# get info ----------------------------------------------------------------

# the 'summary' part of the filename doesn't really mean anything here
info1 <- create_rast_info(r1,
  into = c('type', 'variable', 'RCP', 'years', 'graze','summary')
)


# create maps -------------------------------------------------------------
r <- r_driv1[[1]]
plot_delta_median <- function(r, var, tag = '') {
  stopifnot(nlyr(r) == 1)
  med <- global(r, median, na.rm = TRUE)[[1]]
  diff <- r - med
  limits <- range_raster(diff, absolute = TRUE)
  scale_name <-  labs_pred1[[as.character(var)]]
  p <- plot_map_inset(r = diff,
                      colors = rev(cols_map_bio_d2),
                      tag_label = tag,
                      limits = limits,
                      scale_name = scale_name
  ) 
  p&theme(legend.position = 'left')
  
}


range <- range_raster(r1)

lets <- matrix(fig_letters[1:(nrow(info1)*2)], 
               ncol = 2, byrow = TRUE)

info2 <- info1 %>% 
  df_factor() %>% 
  mutate(variable = driver2factor(variable)) %>% 
  arrange(variable) %>% 
  mutate(tag = lets[, 2],
         tag = paste(tag, 'spatial effect of', variable, 'on fire')) 

limits <- c(-2, 2)

if(max(abs(range)) > max(abs(limits))) {
  message('limits of map values trimmed for easier reading')
}
plots1 <- map2(info2$id, info2$tag, function(id, tag) {
  plot_delta_fire(r = r1[[id]],
                  panel_tag = tag, limits = limits,
                  add_inset = TRUE) +
    theme(legend.position = 'bottom')
})


plots2 <- map2(as.character(info2$variable), lets[,1], function(var, let) {
  tag <- paste(let, 'difference from median', var)
  plot_delta_median(
    r = r_driv1[[var]],
    var = var,
    tag = tag
  )
})

left <- wrap_plots(plots2, ncol = 1)
right <- wrap_plots(plots1, ncol = 1) +
  plot_layout(guides = 'collect')

g <- wrap_plots(left, right, ncol = 2, widths = c(1, 1))

filename <- paste0('figures/fire_attribution/maps/delta-fire-attrib_sref_graz',
                   words2abbrev(ref_graze), '_', rcp_yr, '_', runv, '.png'
)

ggsave(filename,
       plot = g,
       dpi = 600,
       width = 8, 
       height = 13.5)
