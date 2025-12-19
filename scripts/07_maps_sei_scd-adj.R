# Purpose: Create maps of Q and SEI, these are values where delta's from 
# stepwat were applied to SCD cover layers (aggregated to 1km)

# Author: Martin Holdrege

# Started: May 8, 2025

# parameters --------------------------------------------------------------

source('src/params.R')
v <- v_interp
runv <- paste0(run, v)
groups <- c('Sagebrush', 'Pherb', 'Aherb', 'SEI')
test_run <- TRUE #  opt$test_run # 
c3_rcps <- c('RCP45', 'RCP85') # future scenario shown on map
ref_graze <-  opt$ref_graze
target_graze_v <- c('Heavy', 'Very Heavy') # comparison grazing levels
years <- opt$years
v_interp <- opt$v_interp
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(patchwork)
source("src/general_functions.R")
source("src/fig_params.R")
source("src/mapping_functions.R")
source('src/fig_functions.R')
source('src/SEI_functions.R')

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

# created in "scripts/06_summarize_sei_scd-adj.R"
r_c12a <-  rast(paste0('data_processed/interpolated_rasters/',v_interp,
                       '/', runv, '_c12_median-pw_',  years, '.tif'))

# prepare data ------------------------------------------------------------

if(test_run) {
  r_cov1 <- downsample(r_cov1)
  r_qsei1 <- downsample(r_qsei1)
  r_cov_d1 <- downsample(r_cov_d1)
  r_qsei_d1 <- downsample(r_qsei_d1)
  r_c12a <- downsample(r_c12a)
}


r_comb1 <-c(r_cov1, r_cov_d1, r_qsei1, r_qsei_d1)
into = c("group", "type", "RCP", "years", "graze", "summary")

info1 <- create_rast_info(r_comb1, into = into) %>% 
  filter(years != !!years,
         summary == 'median') 


info_c12 <- create_rast_info(r_c12a,into =  c("type", "RCP", "years", 
                                              "graze", 'summary')) %>% 
  filter(summary == 'median')



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

# * climate c9 ------------------------------------------------------------
# include map of change in sei

# see pre December 2025 commits for code to make this figure
# four panels, top row is is historical c3, bottomr row is future c9,
# columns are two grazing levels
# filename <- paste0('figures/sei/maps/4panel_c3c9_',
#                    c3_rcp, '_',
#                    paste(names(c3_graze), collapse = ''), '_scd-adj_', 
#                    runv, ".png")



# * grazing c9 --------------------------------------------------------------
# change sei class relative to moderate grazing (within a climate scenario)


target_grazeL <- list(
  c('Light', 'Heavy', 'Very Heavy'),
  'Very Heavy')



for(rcp in c3_rcps) {
  for(target_graze in target_grazeL) {
    
    c3_info <- info1 %>% 
      filter(type == 'SEI',
             graze %in% c(ref_graze, target_graze),
             RCP %in% c('Current', rcp)) %>% 
      select(-run2, -group, -type)
    
    c3_info2 <- c3_info %>% 
      filter(graze == ref_graze) %>% 
      select(-graze) %>% 
      left_join(c3_info[c3_info$graze != ref_graze, ],
                by = join_by(run, RCP, years, summary),
                suffix = c('_ref', '_target'))
    
    r_c3 <- sei2c3(r_comb1[[c3_info$id[c3_info$graze == ref_graze]]])
    c3_target <- sei2c3(r_comb1[[c3_info2$id_target]])
    r_c9 <- c3toc9(r_c3[[c3_info2$id_ref]], c3_target)
    
    info_c3 <- c3_info[c3_info$graze == ref_graze, ]
    info_c9 <- c3_info[c3_info$graze != ref_graze, ]
    
    g <- plot_c3c9_gref(r_c3 = r_c3, r_c9 = r_c9,
                        info_c3 = info_c3, info_c9 = info_c9)
    
    n <- nrow(c3_info)
    filename <- paste0('figures/sei/maps/c3c9gref_', n, 'panel_',
                       rcp, '_',
                       paste0(words2abbrev(target_graze), collapse = ''), '_', 
                       runv, ".png")
    
    width <- if(n == 8) {
      8.5
    } else if (n == 4) {
      4.5
    } else {
      stop('set fig width')
    }
    
    png(filename, res = 800, height = 6.2, width = width, units = 'in')
    print(g)
    dev.off()
  }
}


# * climate & grazing -----------------------------------------------------
# 3 panel maps showing sei class (historical climate, reference grazing)
# as well as one map of change due to grazing and one map of change
# due to climate

target_rcp <- 'RCP45'
target_yr <- years
for (target_graze in target_graze_v) {

  
  g <- plot_c3c9_3panel(r_sei = r_comb1,
                        info = info1,
                        ref_graze = ref_graze,
                        target_rcp = target_rcp,
                        target_yr = target_yr,
                        target_graze = target_graze,
                        include_cgref = FALSE)
  
  filename <- paste0('figures/sei/maps/c3c9cgref_3panel_',
                     target_rcp, '_',target_yr, '_', words2abbrev(target_graze), 
                     '_', runv, ".pdf")
  ggsave(
    filename = filename,
    plot = g,
    height = 6.8, width = 7.4,
    device = cairo_pdf, # so text can be edited
    family = "sans"  
  )
  
  # 4 panel version (that includes cgref--the combined 
  # grazing and climate change effect)
  g <- plot_c3c9_3panel(r_sei = r_comb1,
                        info = info1,
                        ref_graze = ref_graze,
                        target_rcp = target_rcp,
                        target_yr = target_yr,
                        target_graze = target_graze,
                        include_cgref = TRUE)
  filename <- paste0('figures/sei/maps/c3c9cgref_4panel_',
                     target_rcp, '_',target_yr, '_', words2abbrev(target_graze), 
                     '_', runv, ".pdf")
  ggsave(
    filename = filename,
    plot = g,
    height = 6.8, width = 7.4,
    device = cairo_pdf, # so text can be edited
    family = "sans"  
  )
}


# **c12 -------------------------------------------------------------------
# 4 panel showing c12 change (i.e. 'stable' class broken down
# into areas with stable SEI and those with decline SEI)


