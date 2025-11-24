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

# read in data ------------------------------------------------------------
# this rasters created in 03_interpolation.R
# but using inputs from 06_fire_attribion_sref.R (i.e.,
# unlike the normal approach,
# the 03_ needs to re-run after the 06_ to create these files)

files <- list.files(
  file.path("data_processed/interpolated_rasters/fire_attrib_sref", v_interp),
  pattern = paste0(runv, ".*", ref_graze),
  full.names = TRUE)

# right now the code is setup for running for reference
# conditions (historical climate, ref graze) only
stopifnot(length(files) == 5)

r1 <- rast(files)

if(test_run) {
  r1 <- downsample(r1)
}

# get info ----------------------------------------------------------------

# the 'summary' part of the filename doesn't really mean anything here
info1 <- create_rast_info(r1,
  into = c('type', 'variable', 'RCP', 'years', 'graze','summary')
)


# create maps -------------------------------------------------------------

range <- range_raster(r1)

info2 <- info1 %>% 
  df_factor() %>% 
  mutate(variable = driver2factor(variable)) %>% 
  arrange(variable) %>% 
  mutate(tag = fig_letters[1:n()],
         tag = paste(tag, 'spatial effect of', variable)) 

limits <- c(-2, 2)

if(max(abs(range)) > max(abs(limits))) {
  message('limits of map values trimmed for easier reading')
}
plots1 <- map2(info2$id, info2$tag, function(id, tag) {
  plot_delta_fire(r = r1[[id]],
                  panel_tag = tag, limits = limits)
})

design <- '
    ABC
    DEF
  '

plots2 <- c(plots1, list(patchwork::guide_area()))
g <- patchwork::wrap_plots(plots2, design = design) + 
  plot_layout(guides = 'collect')

filename <- paste0('figures/fire_attribution/maps/delta-fire-attrib_sref_graz',
                   words2abbrev(ref_graze), '_', unique(info2$RCP), '_',
                   unique(info2$years), '_', runv, '.png'
)

ggsave(filename,
       plot = g,
       dpi = 600,
       width = 8, 
       height = 6)
