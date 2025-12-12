# Purpose:
# calculate five number summaries of Qs and SEI (remote sensing adjusted) values
# for creating boxplots and related figures for each ecoregion

# Author: Martin Holdrege
# Started: May 12, 2025

# params ------------------------------------------------------------------

source("src/params.R")
test_run <- opt$test_run
run <- opt$run
v_interp <- opt$v_interp
runv <- paste0(paste0(opt$run), v_interp)
vr <- opt$vr
vr_name <- opt$vr_name
yr_lab <- opt$yr_lab
fire_breaks <- c(seq(0, 2.5, by = 0.25), 3, 4) # intervals to cut fire probability into
ref_graze <- opt$ref_graze
path_tmp <- file.path('tmp', yr_lab)
rerun <- FALSE
dir <- file.path('tmp', yr_lab)
# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
source('src/general_functions.R')
source('src/mapping_functions.R')
source("src/SEI_functions.R")
source("src/distrib_functions.R")

# read in data ------------------------------------------------------------

 # file create in 05_interpolated_summarize_sei_scd-adj.R
# r_qsei1 <- rast(file.path("data_processed/interpolated_rasters", v_interp,
#                      paste0(runv, yr_lab, "_q-sei_scd-adj_summary.tif")))

r_q1 <- rast(file.path(path_tmp, 'q_scd_adj1.tif')) # q layers for each GCM
r_sei1 <- rast(file.path(path_tmp, 'sei_scd_adj1.tif')) # sei layers for each GCM
r_qsei1 <- c(r_q1, r_sei1)
eco1 <- load_wafwa_ecoregions(wafwa_only = FALSE, total_region = TRUE,
                              v = vr)
eco2 <- vect(eco1['ecoregion'])

r_eco1 <- load_wafwa_ecoregions_raster(wafwa_only = FALSE,
                                       v = vr)

size <- load_cell_size()

# * setup -----------------------------------------------------------------

if(test_run) {
  r_qsei1 <- downsample(r_qsei1)
  r_eco1 <- downsample(r_eco1)
  size <- downsample(size)
}

# names of layers --------------------------------------------------

into <- c("group", "type", "RCP", "years", "graze", 'GCM')

info1 <- create_rast_info(r_qsei1, into = into) %>% 
  # the 'biomass' layers are actually Q scores 
  mutate(type = ifelse(type == 'biomass', 'Q', type))


# calculate change layers -------------------------------------------------
# these are for calculate summary stats needed for histograms

info_sei <- info1 %>% 
  filter(type == 'SEI') %>% 
  select(-group)

# not including r_eco4 in the digest b/ 
# ecoregion summaries are not saved (re-run every time, and this
# way not everything will rerun if use new regions--change
# if code adjusted)
hash <- digest::digest(list(r_qsei1, r_eco1, size))

path_hash <- paste0(dir, "/06_hash", ".txt")

if(!file.exists(path_hash)) {
  dir.create('tmp/', showWarnings = FALSE)
  write_lines(hash, path_hash)
  rerun <- TRUE
} else {
  old_hash <- read_lines(path_hash)
  
  if(old_hash != hash & !test_run) {
    write_lines(hash, path_hash)
    rerun <- TRUE # have to rerun if the hash is different (inputs don't match)
  } else if (test_run) {
    rerun <- TRUE
  }
}

if(tmp_exists('r_sei_cref', rerun, dir = dir)) {
  r_sei_cref <- read_tmp_tif('r_sei_cref', dir = dir)
  r_sei_gref <- read_tmp_tif('r_sei_gref', dir = dir)
} else {
  # change due to climate
  r_sei_cref <- calc_rast_cref(r_qsei1, info = info_sei, 
                               type_from = 'SEI_SEI',
                               type_to = 'SEI-cref',
                               by = c('run', 'type', 'graze'))
  
  # change due to grazing
  r_sei_gref <- calc_rast_gref(r_qsei1, info = info_sei, 
                               ref_graze = ref_graze,
                               type_from = 'SEI_SEI',
                               type_to = 'SEI-gref',
                               by = c('run', 'type', 'RCP', 'years', 'GCM'))
  if(!test_run) {
    r_sei_cref <- writeReadRast(r_sei_cref, 'r_sei_cref', dir = dir)
    r_sei_gref <- writeReadRast(r_sei_gref, 'r_sei_gref', dir = dir)
  }

}


# calculate distribution summaries --------------------------------------------



p1_cref <- paste0('data_processed/raster_means/', runv, '_', vr,
              '_', years, '_sei_hist-data_by-gcm-ecoregion_cref.rds')
p1_gref <- paste0('data_processed/raster_means/', runv, '_', vr,
                  '_', years, 
                  '_sei_hist-data_by-gcm-ecoregion_gref.rds')

if(file.exists(p1_cref) & file.exists(p1_gref) & !rerun) {
  # nothing needs to be done
} else {
  breaks <- seq(-1, 1, length.out = 201)
  hist_sei_cref <- hist_from_rast(r = c(r_eco1, r_sei_cref), breaks = breaks,
                                  by_ecoregion = TRUE)
  hist_sei_gref <- hist_from_rast(c(r_eco1, r_sei_gref), breaks = breaks,
                                  by_ecoregion = TRUE)
  
  if(!test_run) {
    saveRDS(hist_sei_cref, p1_cref)
    
    saveRDS(hist_sei_gref, p1_gref)
  }

}

# calculate means by ecoregions -------------------------------------------

# taking  unweighted summary stats is ok, because we're using an equal area projection
sixnum <- function(x, na.rm = TRUE) {
  five <- fivenum(x, na.rm = na.rm)
  c(five, mean(x, na.rm = na.rm))
}
sixnum_names <- c('min', 'lower', 'middle', 'upper', 'max', 'mean')

# doing several layers at a time 
# because getting memory errors
nms <- names(r_qsei1)
chunk_size <- 50
nms_chunks <- split(nms, ceiling(seq_along(nms) / chunk_size))


if(tmp_exists('tmp_l_sixnum', rerun, dir = dir, suffix = '.rds')) {
  tmp_l_sixnum <- readRDS(file.path(dir, 'tmp_l_sixnum.rds')) 
} else {
  tmp_l_sixnum <- map(nms_chunks, function(nm) {
    as_tibble(terra::extract(r_qsei1[[nm]], eco2, sixnum, na.rm = TRUE))
  })
  
  if(!test_run) {
    saveRDS(tmp_l_sixnum, file.path(dir, 'tmp_l_sixnum.rds')) 
  }
}

tmp <- reduce(tmp_l_sixnum, .f = left_join, by = 'ID')

gcms <- str_replace_all(unique(info1$GCM), '\\-', '\\.')
regex <- '\\.\\d$'
regex <- paste0('(?<=(', paste0(gcms, collapse = ')|('), '))\\.\\d$')

eco_gcm1 <- as_tibble(tmp) %>% 
  pivot_longer(cols = -ID) %>% 
  mutate(id = str_replace(name, regex, ''),
         # output of terra::extract had the '-' replaced by '.', 
         id = str_replace_all(id, '\\.', '\\-'),
         stat_num = str_extract(name, regex),
         stat_num = str_replace(stat_num, '^\\.', ''),
         stat_num = as.numeric(stat_num),
         stat_num = ifelse(is.na(stat_num), 0, stat_num) + 1,
         stat = sixnum_names[stat_num],
         region = eco2$ecoregion[ID]) 

check <- eco_gcm1 %>% 
  group_by(id, region) %>% 
  summarize(n = n(), .groups = 'drop') %>% 
  pull(n) %>% unique()

stopifnot(check == length(sixnum_names),
          all(complete.cases(eco_gcm1)))

eco_gcm2 <- eco_gcm1 %>% 
  select(-stat_num, -ID, -name) %>% 
  pivot_wider(values_from = 'value',
              names_from = 'stat') %>% 
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


# gw summaries ------------------------------------------------------------

# gcm wise summaries 

group_vars <- c('region', 'group', 'type', 'RCP', 'years', 'graze')
gcm_mean <- eco_gcm2 %>% 
  group_by(across(all_of(group_vars))) %>% 
  summarize_across_GCMs(var = 'mean') 

eco_smry_gw0 <- eco_gcm2 %>% 
  right_join(gcm_mean, 
             by = c(group_vars, 'mean')) 

check <- eco_smry_gw0 %>% 
  group_by(across(all_of(group_vars))) %>% 
  summarize(n = n(), .groups = 'drop') %>% 
  pull(n)

stopifnot(all(check == 3)) # check for non-unique GCMs

eco_smry_gw <- eco_smry_gw0%>% 
  rename(gcm_srmy = GCM)

tmp <- eco_smry_gw0 %>% 
  select(region, RCP, years, graze, GCM) %>% 
  distinct()

sei_pcent_gw <- left_join(tmp, sei_pcent) %>% 
  rename(gcm_srmy = GCM)

sei_pcent_long_gw <- left_join(tmp, sei_pcent_long) %>% 
  rename(gcm_srmy = GCM)

# save output -------------------------------------------------------------

prefix <- if(test_run) 'test' else paste0(runv, vr_name, yr_lab)

write_csv(eco_smry_gw, paste0('data_processed/raster_means/', prefix, 
                            '_q-sei_scd-adj_smry-gw_by-ecoregion.csv'))

write_csv(sei_pcent_gw , paste0('data_processed/raster_means/', prefix, 
                            '_sei-class-pcent_scd-adj_smry-gw_by-ecoregion.csv'))

write_csv(sei_pcent_long_gw, paste0('data_processed/raster_means/', prefix, 
                             '_sei-class-pcent-long_scd-adj_smry-gw_by-ecoregion.csv'))



# see pre Dec 2025 commits for code to create this output
# write_csv(df_seifire3, paste0('data_processed/raster_means/', prefix, 
#                               '_sei-by-fire-bin_scd-adj_summaries_by-ecoregion.csv'))

# see pre Dec 2025 commits for code to create
# write_csv(c9_area2, paste0('data_processed/raster_means/', prefix, 
#                               '_c9-area_scd-adj_smry-gw_by-ecoregion.csv'))
