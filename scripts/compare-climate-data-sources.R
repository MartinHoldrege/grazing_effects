# Purpose: 
# compare climate data from stepwat with climate data from gridMet/MACA, that
# is part of the new RR data release from DS
# for the 200 stepwat simulations, particularly for determinging why 
# the spatial patterns in the correlation between montly temp and precip
# seem different between the two data sources. If the data don't match well for the 200 
# sites that means it isn't an interpolation problem

# also compare interpolated values to gridded climate product to see errors
# introduced by interpolation

# Author: Martin Holdrege

# script started: May 24, 2024


# dependencies ------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(patchwork)
source('src/paths.R')
source('src/mapping_functions.R')
theme_set(theme_bw())

# read in data ------------------------------------------------------------

loc1 <- read_csv('data_raw/site_locations.csv') # stepwat sites
sw_clim1 <- read_csv("data_processed/site_means/clim_seasonality.csv")

# data from the RR data release, this is the mean correlation between 
# monthly precip and temperature
rr_cur1 <- rast(file.path(path_large, 'newRR3-analysis/examples/tmp_examples',
                          'Preds19__sc1_NA_ambient_1980-2020-clim__value-sim.tif'))

rr_fut1 <- rast(file.path(path_large, 'newRR3-analysis/examples/tmp_examples',
                          'Preds19__RCP45_2064-2099-clim__value-acrmod-med.tif'))

# climate from the weather data base, summarized in the 00_query_weather_db.R script
db_clim1 <- read_csv('data_processed/site_means/dbWeather_200sites.csv')

# mean montly precip, and min and max temp, from daymet.
# calculated in the daymet_normals.js script, in the SEI repository
day1 <- rast("data_raw/daymet_v4_monthly_normals_1991-2020.tif")

# daymet bioclim variables (file created in interpolation_data.R)
period <- '1991-2020'
day_table1 <- read_csv(paste0("data_processed/interpolation_data/clim_for_interp_", 
                       "daymet-v4_", period, ".csv"),
                       col_types = cols(.default = "d"))

template <- rast("data_processed/interpolation_data/cellnumbers.tif")

# * stepwat maps ------------------------------------------------------------------
# interpolated climate data
# [note the interpolated daymet data is daymet v3]
clim_vars <- c('MAP', 'MAT', 'PTcor')
files <- paste0(clim_vars, "_climate_Current_Current_Current_20230919.tif")

sw_r1 <- rast(file.path('data_processed/interpolated_rasters/climate',
                        files))

names(sw_r1) <- names(sw_r1) %>% 
  str_extract("^[[:alpha:]]+")

# interpolated daymet data (interpolated same way as stepwat output)
# matching criteria based on climate envelope of entire study region
files <- paste0(clim_vars, "_daymet-climate_20240604_orig-criteria.tif")
day_interp_v1 <- rast(file.path('data_processed/interpolated_rasters/climate',
                                files))
# interpolation based on matching criteria using climate envelope of 200 sw sites
files <- paste0(clim_vars, "_daymet-climate_20240605_v2.tif")
day_interp_v2 <- rast(file.path('data_processed/interpolated_rasters/climate',
                                files))
day_interp1 <- c(day_interp_v1, day_interp_v2)

clim_version <- function(x) paste(str_extract(x, "^[[:alpha:]]+"), 
                                  str_extract(x, 'v\\d+$'), sep = "_")
names(day_interp1) <- names(day_interp1) %>% 
  str_replace('orig-criteria', 'v1') %>% 
  clim_version()
  
# * grid-met --------------------------------------------------------------

buffer <- 0
# file created in 00_gridmet_normals.R
grid1 <- rast(paste0("data_processed/interpolation_data/gridmet_bioclim-vars_19912020_",
                     buffer, "buffer.tif"))
# functions ---------------------------------------------------------------

range_lab <- function(x, unit = "") {
  as.numeric(minmax(x)) %>% 
    round() %>% 
    paste(collapse = ' to ') %>% 
    paste0('\n(range:', ., unit, ')')
}

# summarize across GCMs ---------------------------------------------------

sw_clim_med <- sw_clim1 %>% 
  filter(RCP != 'RCP85',
         years != '2030-2060') %>% 
  group_by(site, years, RCP) %>% 
  # PTcor2 is the mean of yearly correlations
  summarise(PTcor_stepwat = median(PTcor),
            MAT_stepwat = median(MAT),
            MAP_stepwat = median(MAP)*10) # median across GCMs, convert to mm

# extract clim at points --------------------------------------------------

# convert site locations to sf object
loc2 <- st_as_sf(loc1[, c('site_id', 'X_WGS84', 'Y_WGS84')],
                 crs = 'EPSG:4326',
                 coords = c('X_WGS84', 'Y_WGS84')) %>% 
  st_transform(crs = crs(rr_cur1)) %>% 
  arrange(site_id)


# * PT corr ---------------------------------------------------------------

terra::activeCat(rr_fut1) <- 'CorTempPPT_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'CorTempPPT_mean'

rr_fut2 <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(PTcor2_rr = CorTempPPT_mean_acrmoddistr)
rr_cur2 <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(PTcor2_rr = CorTempPPT_mean)


# * MAP -------------------------------------------------------------------

terra::activeCat(rr_fut1) <- 'PPT_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'PPT_mean'

rr_fut2a <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(MAP_rr = PPT_mean_acrmoddistr)
rr_cur2a <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(MAP_rr = PPT_mean)

# * MAT -------------------------------------------------------------------

terra::activeCat(rr_fut1) <- 'Tmean_mean_acrmoddistr'
terra::activeCat(rr_cur1) <- 'Tmean_mean'

rr_fut2b <- terra::extract(rr_fut1, loc2) %>% 
  mutate(RCP = 'RCP45') %>% 
  rename(MAT_rr = Tmean_mean_acrmoddistr)
rr_cur2b <- terra::extract(rr_cur1, loc2) %>% 
  mutate(RCP = 'Current') %>% 
  rename(MAT_rr = Tmean_mean)


# * database --------------------------------------------------------------

stopifnot(! 'RCP' %in% names(db_clim1)) # db_clim1 is just for current climate
db_clim2 <- db_clim1 %>% 
  select(Site_id, 'MAT_C', 'MAP_mm', matches('CorrTP')) %>% 
  rename( 'MAT_db' = 'MAT_C', 'MAP_db' = 'MAP_mm', PTcor2_db = 'CorrTP2',
          PTcor_db = 'CorrTP') %>% 
  mutate('RCP' = 'Current')

# *gridmet ----------------------------------------------------------------

grid3 <- terra::project(grid1, crs(loc2))

grid4 <- terra::extract(grid3, loc2)
grid5 <- grid4 %>% 
  mutate(RCP = 'Current') %>% 
  rename('MAT_grid' = 'bio1',
         'MAP_grid' = 'bio12') %>% 
  select(-matches('bio\\d+'))

# * daymet ----------------------------------------------------------------

day2 <- terra::extract(day1, loc2)

day3 <- day2 %>% 
  pivot_longer(cols = -ID,
               names_to = c("name", 'month'),
               names_sep = "_") %>% 
  pivot_wider(values_from = 'value',
              names_from = 'name') %>% 
  mutate(tmean = (tmax + tmin)/2) %>% 
  group_by(ID) %>% 
  summarize(PTcor_day = cor(tmean, prcp),
            MAP_day = sum(prcp),
            MAT_day = mean(tmean)) %>% 
  mutate(RCP = 'Current')


# * combine ---------------------------------------------------------------

comb <- sw_clim_med %>% 
  left_join(bind_rows(rr_fut2, rr_cur2), by = c('site' = 'ID', 'RCP'))%>% 
  left_join(bind_rows(rr_fut2a, rr_cur2a), by = c('site' = 'ID', 'RCP')) %>% 
  left_join(bind_rows(rr_fut2b, rr_cur2b), by = c('site' = 'ID', 'RCP')) %>% 
  left_join(db_clim2, by =  c('site' = 'Site_id', 'RCP')) %>% 
  left_join(day3, by = c('site' = 'ID', 'RCP')) %>% 
  left_join(grid5, by = c('site' = 'ID', 'RCP'))

comb_long <- comb %>% 
  #filter(RCP == 'Current') %>% 
  pivot_longer(names_to = c('variable', 'data_source'),
               cols = -c(site, years, RCP),
               names_sep = "_") %>% 
  pivot_wider(names_from = 'data_source') %>% 
  ungroup()


# compare rasters ---------------------------------------------------------
# stepwat interpolated climate data vs daymet v4 climate data


# *prep daymet data --------------------------------------------------------
day1 <- extend(day1, sw_r1) %>% 
  crop(sw_r1) %>% 
  project(crs(sw_r1))
nms <- names(day1) %>% 
  sort()

tmean <- (day1[[str_subset(nms, '^tmax')]] + day1[[str_subset(nms, '^tmin')]])/2
prcp <- day1[[str_subset(nms, '^prcp')]]

# P-T correlation (type 1)
day_ptcor <- app(c(tmean, prcp), fun = function(x) {
  # first 12 layers are temp, remains ones are prcp
  cor(x[1:12], x[13:24])
})

day_map <- app(prcp, fun = sum)
day_mat <- app(tmean, fun = mean)

# * differences -----------------------------------------------------------

# stewpat minus daymet differences
map_diff_sd <- sw_r1[['MAP']] - day_map
mat_diff_sd <- sw_r1[['MAT']] - day_mat
ptcor_diff_sd <- sw_r1[['PTcor']] - day_ptcor

# daymet (interpolated) minus daymet
day_interp2 <- crop(day_interp1, sw_r1)
map_diff_dd <- day_interp2[[c('MAP_v1', "MAP_v2")]]  - day_map
mat_diff_dd <- day_interp2[[c('MAT_v1', "MAT_v2")]]  - day_mat
ptcor_diff_dd <- day_interp2[[c('PTcor_v1', "PTcor_v2")]]  - day_ptcor

diff_dd <- c(map_diff_dd, mat_diff_dd, ptcor_diff_dd)
# figures -----------------------------------------------------------------

# *comparing stepwat weather to gridmet -----------------------------------

cap2 <- "'type 2' is the correlation between monthly values in a given yr, then avaraged\n type 1 is correlation of monthly means"


g0 <- ggplot(comb, aes(PTcor_stepwat)) +
  geom_histogram() +
  facet_wrap(~RCP)+
  labs(x = 'PT correlation (STEPWAT2, type 1)')

g1 <- ggplot(comb, aes(PTcor2_rr)) +
  geom_histogram() +
  facet_wrap(~RCP) +
  labs(x = 'PT correlation (RR dataset, type 2)')

g3 <- g0/g1

g4 <- ggplot(comb, aes(PTcor2_rr, PTcor2_db)) +
  geom_point() +
  labs(x = 'GridMet/MACA (type 2)',
       y = 'Weather database (from 200 sites, type 2)',
       caption = cap2,
       subtitle = 'Correlation (type 2) between monthly temperature and precip,\nfor 200 sites at which STEPWAT simulations are conducted') +
  geom_abline(slope = 1)

g5 <- ggplot(comb, aes(PTcor_db, PTcor2_db)) +
  geom_point() +
  labs(x = 'weather database (type 1)',
       y = 'weather database (type 2)',
       caption = cap2,
       subtitle = 'Correlation (type 2) between monthly temperature and precip,\nfor 200 sites at which STEPWAT simulations are conducted') +
  geom_abline(slope = 1)

# *comparing weather db to stepwat output ----------------------------------
# the climate from the weather db should be closely comparable to what comes
# out of stepwat (which has gone through the wx generator)

base <- function() {
  list(
    geom_abline(slope = 1),
    geom_point(),
    labs(subtitle = 'comparing climate metrics from STEPWAT output\n vs directly from the weather database')
  )
}

pdf('figures/climate/clim_comparison-between-sources.pdf',
    width = 8, height = 8)

g3
g4
g5

ggplot(comb, aes(PTcor_db, PTcor_stepwat))+
  base() +
  labs(caption = cap2,
       x = 'PT cor from database (type 1)',
       y = 'PT cor from stepwat (type 1)')

ggplot(comb, aes(PTcor2_db, PTcor_stepwat))+
  base() +
  labs(caption = cap2,
       x = 'PT cor from database (type 2)',
       y = 'PT cor from stepwat (type 1)')


ggplot(comb, aes(PTcor_db, PTcor2_db))+
  base() + 
  labs(caption = cap2,
     x = 'PT cor from database (type 1)',
     y = 'PT cor from database (type 2)',
     subtitle = 'weather db data (comparing summary methods)')

# * compare all datasets ------------------------------------------

comb_long %>% 
  select(variable, stepwat, day, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(day, stepwat)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free',
             ncol = 2) +
  geom_abline(slope = 1) +
  labs(x = 'Daymet V4',
       y = 'STEPWAT')

comb_long %>% 
  select(variable, stepwat, rr, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(rr, stepwat)) +
  geom_point() +
  facet_wrap(RCP~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'RR dataset',
       y = 'STEPWAT')

comb_long %>% 
  select(variable, stepwat, grid, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(grid, stepwat)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'GridMet',
       y = 'STEPWAT')

comb_long %>% 
  select(variable, stepwat, db, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(db, stepwat)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'Database for STEPWAT',
       y = 'STEPWAT')

comb_long %>% 
  select(variable, day, db, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(day, db)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'DayMet V4',
       y = 'Database for STEPWAT')

comb_long %>% 
  select(variable, grid, day, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(day, grid)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'DayMet V4',
       y = 'GridMet')

comb_long %>% 
  select(variable, rr, day, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(rr, day)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'RR dataset',
       y = 'Daymet V4')

comb_long %>% 
  select(variable, rr, grid, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(rr, grid)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'RR dataset',
       y = 'GridMet')


# * maps ------------------------------------------------------------------
cols_diff <- rev(RColorBrewer::brewer.pal(9, 'RdYlBu'))
cols_diff[5] <- 'grey'

diff_maps <- function(r_map, r_mat, r_ptcor) {
  m1 <- plot_map_inset(r = r_map,
                       colors = cols_diff,
                       limits = c(-300, 300),
                       scale_name = 'mm',
                       tag_label = paste('Interpolated minus DayMet MAP',
                                         range_lab(r_map, 'mm')))
  
  m2 <- plot_map_inset(r = r_mat,
                       colors = cols_diff,
                       limits = c(-5, 5),
                       scale_name = 'C',
                       tag_label = paste('Interpolated minus DayMet MAT',
                                         range_lab(r_mat, 'C')))
  
  m3 <- plot_map_inset(r = r_ptcor,
                       colors = cols_diff,
                       limits = c(-1, 1),
                       scale_name = '',
                       tag_label = paste('Interpolated minus DayMet P-T correlation',
                                         range_lab(r_ptcor)))
  design <- "
  12
  3#
"
  out <- m1 + m2 + m3 + plot_layout(design = design)
  out
}

diff_maps(map_diff_sd, mat_diff_sd, ptcor_diff_sd) + 
  plot_annotation(subtitle = 'Interpolated STEPWAT2 output',
                  caption = " ")

diff_maps(diff_dd[['MAP_v1']], diff_dd[['MAT_v1']], diff_dd[['PTcor_v1']]) + 
  plot_annotation(subtitle = 'Interpolated Daymet V4 values',
                  caption = 'Matching criteria based on entire study area')

diff_maps(diff_dd[['MAP_v2']], diff_dd[['MAT_v2']], diff_dd[['PTcor_v2']]) + 
  plot_annotation(subtitle = 'Interpolated Daymet V4 values',
                  caption = 'Matching criteria based on 200 sites')

dev.off()

# climate envelope --------------------------------------------------------
# figures showing study are climate envelope along with the locations of the 200
# points (to see how representative they are)

day_table2 <- day_table1 %>% 
  # prop_scd is the proportion of the 1km grid cell fill w/ 30 cm scd cells
  dplyr::select(bio1, bio12, ptcor, site_id, matches('^prop')) %>% 
  rename(MAT = bio1, MAP = bio12, PTcor = ptcor)

day_sites1 <- day_table2 %>% 
  filter(!is.na(site_id))


# * figures ---------------------------------------------------------------

b <- c(0, 1e-7, 0.001, 0.01, 0.1, 0.2, 0.3,0.5,1)
base_density <- function(breaks = b) {

  colors <- c('transparent',
              RColorBrewer::brewer.pal(length(breaks), 'Greens')[-c(1:2)])
  
  list(geom_density_2d_filled(breaks = breaks,
                              # density is fraction of max
                              contour_var = 'ndensity'),
       scale_fill_manual(values = colors),
       theme(legend.position = 'none')
       )
}

mask_descripts <- 
  expand_grid(sage_gt = c(NA, 1, 2, 5),
              prop_scd = c(0, 0.1, 0.5)) %>% 
  mutate(description = ifelse(is.na(sage_gt),
                              paste0('1 km pixels containing \n>', 
                                     prop_scd*100, '% 30m SCD pixels.'),
                              paste0('1 km pixels containing \n>', prop_scd*100, 
                                     '% of 30m SCD pixels that have >', 
                                     sage_gt, '% sage cover.')))

# '1 km cells containing > 10% SCD 30 m cells'
figs <- pmap(mask_descripts, function(prop_scd, sage_gt, description) {
  prop <- prop_scd
  study_area <- description
 if(is.na(sage_gt)) {
    study_area <- paste0('1 km pixels contaning >', prop*100, '%\n 30m SCD pixels.')
    filter_var <- "prop_scd_all"
  } else {
    study_area <- paste0('1 km pixels contaning >', prop*100, 
                         '%\n of 30m SCD pixels that have >', sage_gt, '% sage cover.')
    filter_var <- paste0('prop_sage_gt_', sage_gt)
  }

  df <- day_table2 %>% 
      filter(.data[[filter_var]] > prop)

  out <- list()
  a <- ggplot(df, aes(MAT, MAP)) +
    base_density() +
    geom_point(data = day_sites1, alpha = 0.3) +
    coord_cartesian(xlim = range(day_table2$MAT, na.rm = TRUE),
                    ylim = range(day_table2$MAP, na.rm = TRUE))

  b <- ggplot(df, aes(MAT, PTcor)) +
    base_density(breaks = b*3) +
    geom_point(data = day_sites1, alpha = 0.3) +
    labs(y = 'Correlation between monthly precip and temp')+
    coord_cartesian(xlim = range(day_table2$MAT, na.rm = TRUE),
                    ylim = range(day_table2$PTcor, na.rm = TRUE))

  c <- ggplot(df, aes(MAP, PTcor)) +
    base_density() +
    geom_point(data = day_sites1, alpha = 0.3)  +
    labs(y = 'Correlation between monthly precip and temp')+
    coord_cartesian(xlim = range(day_table2$MAP, na.rm = TRUE),
                    ylim = range(day_table2$PTcor, na.rm = TRUE))
  design <- "
    12
    3#
    "
  
  caption = paste('Points show 200 STEPWAT2 sites. Data from DayMet V4.',
                       '\nshading shows climate envelope of the',
                       study_area)
  a + b + c + plot_layout(design = design) +
    plot_annotation(caption = caption) +
    theme(plot.caption = element_text(size = rel(0.6)))

})

g5 <- as.data.frame(grid3[[c('bio1', 'bio12')]]) %>% 
  ggplot(aes(bio1, bio12)) +
  base_density() +
  geom_point(data = grid5, aes(MAT_grid, MAP_grid), alpha = 0.3) +
  coord_cartesian(xlim = range(day_table2$MAT, na.rm = TRUE),
                  ylim = range(day_table2$MAP, na.rm = TRUE)) +
  labs(subtitle = paste0('Gridmet data, from entire sagebrush polygon (', 
                         buffer, ' km buffer)'),
       x = 'MAT', y = 'MAP')

pdf('figures/climate/climate-envelope_v3.pdf', height = 8, width = 8)
figs
g5
dev.off()

# maps of masks -----------------------------------------------------------

figs <- pmap(mask_descripts, function(prop_scd, sage_gt, description) {
  prop <- prop_scd
  study_area <- description
  if(is.na(sage_gt)) {
    study_area <- paste0('1 km pixels contaning >', prop*100, '%\n 30m SCD pixels.')
    filter_var <- "prop_scd_all"
  } else {
    study_area <- paste0('1 km pixels contaning >', prop*100, 
                         '%\n of 30m SCD pixels that have >', sage_gt, '% sage cover.')
    filter_var <- paste0('prop_sage_gt_', sage_gt)
  }
  
  df <- day_table1 %>% 
    filter(.data[[filter_var]] > prop)
  
  r <- fill_raster(df = df[, c('cellnumber', filter_var)],
              template = template)

  g <- plot_map_inset(r, colors = cols_diff,
                 limits = c(0, 1),
                 scale_name = 'proportion',
                 tag_label = 'Proportion of 1km cell filled with 30 m cells')
  g + plot_annotation(caption = study_area)
  
})

pdf('figures/study_area/potential_masks.pdf', height = 8, width = 8)
figs
dev.off()