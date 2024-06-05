# Purpose: 
# compare climate data from stepwat with climate data from gridMet/MACA, that
# is part of the new RR data release from DS
# for the 200 stepwat simulations, particularly for determinging why 
# the spatial patterns in the correlation between montly temp and precip
# seem different between the two data sources. If the data don't match well for the 200 
# sites that means it isn't an interpolation problem

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
day1 <- rast("data_raw/daymet_monthly_normals_1981-2010.tif")

# daymet bioclim variables (file created in interpolation_data.R)
day_table1 <- read_csv("data_processed/interpolation_data/clim_for_interpolation.csv",
                       col_types = cols(.default = "d"))

# * stepwat maps ------------------------------------------------------------------
# interpolated climate data
files <- paste0(c('MAP', 'MAT', 'PTcor'), "_climate_Current_Current_Current_20230919.tif")

sw_r1 <- rast(file.path('data_processed/interpolated_rasters/climate',
                        files))

names(sw_r1) <- names(sw_r1) %>% 
  str_extract("^[[:alpha:]]+")


# * grid-met --------------------------------------------------------------
# ~ 4 km gridded data (from Daniel)
grid1 <- rast(file.path(
  path_large,
  "climate_data/gridmet-clim",
  c('gridmet_ppt_annualClimatologies.nc', 'gridmet_tmmn_annualClimatologies.nc',
    'gridmet_tmmx_annualClimatologies.nc')
))


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

# some how x and y coords got reversed
grid2 <- terra::project(grid1, crs(grid1)) %>% 
  t() %>% 
  flip(direction = 'vertical') %>% 
  flip(direction = 'horizontal')
grid3 <- terra::project(grid2, crs(loc2))

grid4 <- terra::extract(grid3, loc2)
grid5 <- grid4 %>% 
  mutate(MAT_grid = (tmmn + tmmx)/2 -273.15,
         RCP = 'Current') %>% 
  rename('MAP_grid' = 'ppt') %>% 
  select(-tmmn, -tmmx) 

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
# stepwat interpolated climate data vs daymet v3 climate data


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

# Next--check units on sw_r1 MAP and MAT
# create 3 maps

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
  labs(x = 'Daymet V3',
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
  labs(x = 'DayMet V3',
       y = 'Database for STEPWAT')

comb_long %>% 
  select(variable, grid, day, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(day, grid)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'DayMet V3',
       y = 'GridMet')

comb_long %>% 
  select(variable, rr, day, RCP) %>% 
  drop_na() %>% 
  ggplot(aes(rr, day)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_abline(slope = 1) +
  labs(x = 'RR dataset',
       y = 'Daymet V3')

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

m1 <- plot_map_inset(r = map_diff_sd,
               colors = cols_diff,
               limits = c(-300, 300),
               scale_name = 'mm',
               tag_label = paste('Interpolated minus DayMet MAP',
                                 range_lab(map_diff_sd, 'mm')))

m2 <- plot_map_inset(r = mat_diff_sd,
               colors = cols_diff,
               limits = c(-5, 5),
               scale_name = 'C',
               tag_label = paste('Interpolated minus DayMet MAT',
                                 range_lab(mat_diff_sd, 'C')))

m3 <- plot_map_inset(r = ptcor_diff_sd,
               colors = cols_diff,
               limits = c(-1, 1),
               scale_name = '',
               tag_label = paste('Interpolated minus DayMet P-T correlation',
                                 range_lab(ptcor_diff_sd)))
design <- "
  12
  3#
"
m1 + m2 + m3 + plot_layout(design = design)
dev.off()


# climate envelope --------------------------------------------------------
# figures showing study are climate envelope along with the locations of the 200
# points (to see how representative they are)

day_table2 <- day_table1 %>% 
  dplyr::select(bio1, bio12, ptcor, site_id) %>% 
  rename(MAT = bio1, MAP = bio12, PTcor = ptcor)

day_sites1 <- day_table2 %>% 
  filter(!is.na(site_id))
#day_table2 <- sample_n(day_table2, 1e5) # for testing

# * figures ---------------------------------------------------------------

b <- c(0, 1e-7, 0.001, 0.01, 0.1, 0.2, 0.3,0.5,1)
base_density <- function(breaks = b) {

  colors <- c('transparent',
              RColorBrewer::brewer.pal(length(breaks), 'Greens')[-c(1:2)])
  
  list(geom_density_2d_filled(breaks = breaks,
                              # density is fraction of max
                              contour_var = 'ndensity'),
       scale_fill_manual(values = colors),
       labs(caption = paste('Points show 200 STEPWAT2 sites,',
                            '\nshading shows climate envelope of the entire',
                            'study are used for the SCD manuscript interpolations.',
                            '\nData from DayMet V3')),
       theme(legend.position = 'none')
       )
}

pdf('figures/climate/climate-envelope_v1.pdf')

ggplot(day_table2, aes(MAT, MAP)) +
  base_density() +
  geom_point(data = day_sites1, alpha = 0.3) 

ggplot(day_table2, aes(MAT, PTcor)) +
  base_density(breaks = breaks*3) +
  geom_point(data = day_sites1, alpha = 0.3) +
  labs(y = 'Correlation between monthly precip and temp')

ggplot(day_table2, aes(MAP, PTcor)) +
  base_density() +
  geom_point(data = day_sites1, alpha = 0.3)  +
  labs(y = 'Correlation between monthly precip and temp')

dev.off()
