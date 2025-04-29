# Purpose: Create maps of Q and SEI

# Author: Martin Holdrege

# Started: April 22, 2025



# parameters --------------------------------------------------------------

v <- 'v4' # interpolation version
# multiple runs may need to listed if different runs done for different
# grazing levels

run <- 'fire1_eind1_c4grass1_co20_2503'

groups <- c('Sagebrush', 'Pherb', 'Aherb', 'SEI')


# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source("src/general_functions.R")
source("src/fig_params.R")
source("src/mapping_functions.R")
source('src/fig_functions.R')

# read in data ------------------------------------------------------------

# actual (absolute) values
r_abs1 <- rast(
  file.path("data_processed/interpolated_rasters", 
            paste0(run, "_q-sei_future_summary_across_GCMs.tif"))
)

# difference values
r_diff1 <- rast(
  file.path("data_processed/interpolated_rasters", 
            paste0(run, "_q-sei-rdiff-cref_summary.tif"))
)

# prepare data ------------------------------------------------------------

r_comb1 <- c(r_abs1, r_diff1)
into = c("group", "type", "RCP", "years", 
         "graze", "summary")

info1 <- create_rast_info(r_comb1, into = into) %>% 
  mutate(group = ifelse(str_detect(type, 'SEI|sei'), 'SEI', group)) %>% 
  filter(years != '2030-2060' & years != '2031-2060') 
  

# 20 panel map (absolute delta) -------------------------------------------


args <- list(
  groups = c('Sagebrush', 'Pherb', 'Aherb', 'SEI'),
  titles = c(Aherb = 'Annual herbacious quality score',
             Pherb = 'Perennial herbacious quality score',
             Sagebrush = 'Sagebrush quality score',
             SEI = 'Climatic Potential Sagebrush Ecological Integrity'),
  type_absolute = c(Aherb = 'Q',
                    Pherb = 'Q',
                    Sagebrush = 'Q',
                    SEI = 'SEI')
  
)


args$type_diff <- paste0(args$type_absolute, "-rdiff-cref") %>% 
  setNames(names(args$type_absolute))

for(group in args$groups){
  
  info_tmp <- info1 %>% 
    filter(.data$group == !!group)
  
  if(group == 'SEI') {
    palette_absolute = cols_seicont
    values_absolute = vals_seicont
  } else {
    palette_absolute = cols_map_bio(10)
    values_absolute = NULL
  }
  
  g <- plot_map_20panel(
    r = r_comb1,
    info = info_tmp,
    type_absolute = args$type_absolute[[group]],
    type_diff =  args$type_diff[[group]],
    title = args$titles[[group]],
    name4absolute = args$type_absolute[[group]],
    palette_absolute =  palette_absolute,
    values_absolute =  values_absolute,
    palette_diff = cols_map_bio_d2
  )

  filename <- paste0('figures/sei/maps/20panel_', group, '_', run, ".png")
  png20panel(filename)
  print(g)
  dev.off()


}


# 20 panel % change --------------------------------------------------------------

# g <- plot_map_20panel_perc(
#   r = r_comb1,
#   info = info_tmp,
#   type_absolute = args$type_absolute[[group]],
#   type_diff =  args$type_diff[[group]],
#   title = args$titles[[group]],
#   name4absolute = args$type_absolute[[group]],
#   palette_absolute = cols_map_bio(10),
#   palette_diff = cols_map_bio_d2,
#   lims_diff = c(-100, 300)
# )
# 
# filename <- paste0('figures/sei/maps/20panel_', group, '-', 
#                    args$type_absolute[[group]], '-perc_', run, ".png")
# png20panel(filename)
# print(g)
# dev.off()

# for(group in args$groups){
#   
#   info_tmp <- info1 %>% 
#     filter(.data$group == !!group)
#   
#   type_absolute <- args$type_absolute[[group]]
#   type_diff =  args$type_diff[[group]]
#   name4absolute <-  args$type_absolute[[group]]
#   name4diff <- paste("\u0394", name4absolute, '(%)')
#   r_perc <- calc_rast_perc_change(r_comb1, info = info_tmp,
#                                   type_absolute = type_absolute,
#                                   type_diff = type_diff)
#   
#   g <- plot_map_20panel(
#     r = r_perc,
#     info = info_tmp,
#     type_absolute = type_absolute,
#     type_diff =  type_diff,
#     title = args$titles[[group]],
#     name4absolute = name4absolute,
#     name4diff = name4diff,
#     palette_absolute = cols_map_bio(10),
#     palette_diff = cols_map_bio_d2,
#     lims_diff = c(-100, 300),
#     midpoint_diff = 0
#   )
#   
#   filename <- paste0('figures/sei/maps/20panel_', group, '-', 
#                      args$type_absolute[[group]], '-perc_', run, ".png")
#   png20panel(filename)
#   print(g)
#   dev.off()
# }
