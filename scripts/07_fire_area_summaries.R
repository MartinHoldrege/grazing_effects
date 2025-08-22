# Purpose: Calculate the expected burned area changes,
# these are simple summary stats, some to be used in the manuscript

# Author: Martin Holdrege

# Date started: September 21, 2025

# params ------------------------------------------------------------------

source('src/params.R')

run <- opt$run
runv <- paste0(run, opt$v_interp)
v <- 'v2'
vr_name <- opt$vr_name

options(readr.show_col_types = FALSE)

graze_ref <- 'Moderate' # reference grazing level
# dependencies ------------------------------------------------------------

library(tidyverse)
source('src/general_functions.R')
source('src/mapping_functions.R')

# created in "scripts/06_fire_area.R"
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                            v, vr_name, "_", run, ".csv"))

ba_smry1 <- read_csv(paste0("data_processed/area/expected-burn-area_smry_", 
                           v, vr_name, "_", run, ".csv"))

area_eco <- read_csv(paste0("data_processed/area/ecoregion-area_", v, vr_name, ".csv"))

# change in ba due to grazing ---------------------------------------------
# change in burned area, reference is moderate grazing, and
# this is done within a climate future

ba_gcm2 <- ba_gcm1 %>% 
  select(-run2, -type) %>% 
  rename(ba = area) %>% 
  df_factor()

# b/ comparing across grazing scenarios, 
# comparing within a gcm (not -pixelwise gcm summaries) b/
# doesn't  make sense to compare different grazing scenarios across differnt gcms
# 
ba_gcm_delta_graze1 <- ba_gcm2 %>% 
  filter(graze == graze_ref) %>% 
  rename(ba_ref = ba) %>% 
  select(-graze) %>% 
  right_join(filter(ba_gcm2, graze != graze_ref),
             by = join_by(ecoregion, run, RCP, years, GCM)) %>% 
  mutate(ba_delta = ba - ba_ref,
         ba_delta_perc = ba_delta/ba_ref*100)

ba_smry2 <- ba_smry1 %>%  
  select(-run2, -matches('_perc')) %>% 
  rename_with(.fn = \(x) str_replace(x, 'area_(?!total)', 'ba_')) 

ba_smry_delta_clim1 <- ba_smry2 %>% 
  filter(RCP == 'Current') %>% 
  rename(ba_cur = ba_median) %>% 
  select(ecoregion, run, graze, ba_cur) %>% 
  right_join(filter(ba_smry2, RCP != 'Current'), 
             by = join_by(ecoregion, run, graze)) %>% 
  mutate(across(.cols = c('ba_median', 'ba_high', 'ba_low'),
                .fns = list('delta' = \(x) (x - ba_cur),
                            'delta_perc' = \(x) (x - ba_cur)/ba_cur*100)
                ),
         across(starts_with("ba") & !contains("perc"), ~ round(.x, 0)),
         across(contains("delta_perc"), ~ round(.x, 1)))

  

# summarize across GCMs ---------------------------------------------------

ba_smry_delta_graze1 <- ba_gcm_delta_graze1 %>% 
  group_by(ecoregion, run, RCP, years, graze) %>% 
  summarize(across(.cols = c('ba_delta_perc'),
                   .fns = list(
                     'median' = median,
                     'low' = calc_low,
                     'high' = calc_high
                   ))) 
  
join_vars <- c('ecoregion', 'run', 'RCP', 'years', 'graze')
tmp <- select(ba_gcm_delta_graze1, all_of(join_vars), ba_delta_perc, ba_delta)

# selecting the delta ba (absolute change in area) 
# that corresponds with the low, median, high % changes in area
# (note that the median change need not correspond to the median % change
# and I'm using the % change as the primary variable of interest)
ba_smry_delta_graze2 <- ba_smry_delta_graze1 %>% 
  left_join(rename(tmp, ba_delta_median = ba_delta),
            by = c(join_vars, ba_delta_perc_median = 'ba_delta_perc')) %>% 
  left_join(rename(tmp, ba_delta_high = ba_delta),
            by = c(join_vars, ba_delta_perc_high = 'ba_delta_perc')) %>% 
  left_join(rename(tmp, ba_delta_low = ba_delta),
          by = c(join_vars, ba_delta_perc_low = 'ba_delta_perc')) %>% 
  mutate(across(starts_with("ba") & !contains("perc"), ~ round(.x, 0)),
         across(contains("delta_perc"), ~ round(.x, 1)))


# if check's didn't pass then 1:1 joining above failed
check <- filter(ba_smry_delta_graze2, RCP != 'Current')
stopifnot(nrow(ba_smry_delta_graze1) == nrow(ba_smry_delta_graze2),
          all(!is.na(check)))

# save output -------------------------------------------------------------


write_csv(ba_smry_delta_graze2,
          paste0("data_processed/area/ba_summaries/ba_graze-delta_gcm-wise_", 
                vr, "_", runv, ".csv"))

write_csv(ba_smry_delta_clim1,
          paste0("data_processed/area/ba_summaries/ba_clim-delta_pixel-wise_", 
                 vr, "_", runv, ".csv"))
