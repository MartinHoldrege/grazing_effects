# Purpose: Make figures that will be part of an explainer box
# that shows how SEI is calculated


# parameters --------------------------------------------------------------

source('src/params.R')
v <- v_interp
runv <- paste0(run, v)

# dependencies ------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
source('src/SEI_functions.R')
source('src/fig_functions.R')
source('src/mapping_functions.R')
theme_set(theme_custom1())


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
            

r_eco1 <- load_wafwa_ecoregions_raster(wafwa_only = TRUE)
# selecting sites ---------------------------------------------------------
# want to select three locations: one CSA, GOA and ORA. these
# 3 will be plotted on maps showing cover, then show how this is converted
# the Q score and then to SEI


# * select raster layers ----------------------------------------------------


# the bandnames don't include the year month of simulation in the run nmame
runv2 <- run %>% 
  str_replace('_\\d{4}', '') %>% 
  paste0(., v)

scenario <- 'Current_Current_Moderate_median'
lyrs_cov <- names(r_cov1) %>% 
  str_subset(paste0('_cover_', scenario)) %>% 
  rev()

names(lyrs_cov) <- str_extract(lyrs_cov, '[[:alpha:]]+(?=_cover)')
stopifnot(length(lyrs_cov) == 3)

r_cov2 <- r_cov1[[lyrs_cov]]

r_q2 <- r_qsei1[[str_replace(lyrs_cov, 'cover', 'Q')]]

r_sei1 <- r_qsei1[[str_subset(names(r_qsei1), 
                              paste0('SEI_', scenario))]]


# * create dataframe ------------------------------------------------------


df1 <- as.data.frame(c(r_cov2, r_q2, r_sei1, r_eco1))

df2 <- df1 %>% 
  rename_with(
    .fn = function(x) {
      str_replace(x, runv2, '') %>% 
        str_replace(scenario, '') %>% 
        str_replace_all('(^_)|(_$)', '')
    }
  ) %>% 
  rename(SEI = SEI_SEI)

df3 <- df2 %>% 
  mutate(cellnum = as.numeric(row.names(.))) %>% 
  drop_na() %>% 
  as_tibble() %>% 
  mutate(c3 = sei2c3(SEI))


# * select points ---------------------------------------------------------


df4 <- df3 %>% 
  # calculate SEI for the assumption where trees and human mod of a Q of 1
  mutate(SEI3 = Sagebrush_Q*Pherb_Q*Aherb_Q) %>% 
  # only keeping cells where that is approximately true
  filter(abs(SEI -SEI3) < 0.001) %>% 
  select(-SEI3)

# checking if the value is between the 'middle' percentiles
is_middle <- function(x) {
  r <- rank(x)
  r2 <- (r - 1)/(max(r) - 1) # scale 0-1
  r2 > 0.3 & r2 < 0.7
}

# selecting only Q values that follw in the middle,
# doing this so that the conceptual figure will be easier to follow
df5 <- df4 %>% 
  group_by(c3) %>% 
  mutate(across(matches('_Q'), .fns = list('mid' = is_middle))) %>% 
  filter(Sagebrush_Q_mid & Pherb_Q_mid & Aherb_Q_mid)

# now selecting ones so the core as the highest Qs for all 3,
# and ORA has the lowest Qs for everything, again to make
# the conceptual figures easier to understand
d6 <- df5 %>% 
  group_by(cellnum) %>% 
  mutate(max = max(Sagebrush_Q, Pherb_Q, Aherb_Q),
         min = min(Sagebrush_Q, Pherb_Q, Aherb_Q)) %>% 
  ungroup() %>% 
  group_by(c3) %>% 
  filter(
    # for CSA, don't want any 'low' Q scores
    (c3 == 'CSA' & min == max(min)) |
      (c3 == 'ORA' & max == min(max)) | (c3 == 'GOA' & Pherb_Q == min(Pherb_Q)) 
    ) %>% 
  arrange(c3) %>% 
  select(c3, everything(), -max, -min, -matches('_mid'))


# get coordinates
r_template <- r_sei1[[1]]
crs <- crs(r_template)
coords <- terra::xyFromCell(r_template , d6$cellnum) %>% as_tibble()  # cols x, y
d7 <- st_as_sf(bind_cols(coords, d6),
               coords = c('x', 'y'),
               crs = crs)

d7 %>% 
  select(c3, SEI, matches('_Q'))

# # A tibble: 3 × 6
# c3        SEI Sagebrush_Q Pherb_Q Aherb_Q           geometry
# <fct>   <dbl>       <dbl>   <dbl>   <dbl>        <POINT [m]>
# 1 CSA   0.570        0.634    0.909  0.990  (-1127500 2188500)
# 2 GOA   0.325        0.574    0.885  0.641  (-1760500 2445500)
# 3 ORA   0.00438      0.0792   0.605  0.0929 (-1966500 2165500)

# cover maps --------------------------------------------------------------

bb  <- st_bbox(r_template)
# r <- downsample(r_cov2[[1]])
r <- r_cov2[[1]]
width_map <- 4; height_map <- 4
add_c3_points <- function(df, crs, bb) {
  list(
    ggnewscale::new_scale_fill(),# allows a second fill scale
    geom_sf(data = df, mapping = aes(fill = c3),
              size = 5,
              stroke = 2,
              shape = 21), 
    coord_sf(crs = crs,
               xlim = bb[c("xmin","xmax")],
               ylim = bb[c("ymin","ymax")],
               expand = FALSE),
    scale_fill_manual(values = cols_c3, guide = 'none'),
    theme(
        legend.position.inside = c(0.99, 0.01),   # new syntax
        legend.justification   = c(0.99, 0.01),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7)
      )
  )
}

map2(lyrs_cov, names(lyrs_cov), function(lyr, pft) {
  r <- r_cov2[[lyr]]
  g <- plot_map2(r) +
    scale_fill_gradientn(na.value = 'transparent',
                         colors = cols_map_bio(10),
                         name = 'Cover (%)') +
    add_c3_points(df = d7, crs = crs, bb = bb)
  
  ggsave(paste0('figures/sei/conceptual/cover-map_', pft, '_scd-adj_', runv, '_',
                scenario, '.png'),
         plot = g, width = width_map, height = height_map, dpi = 600)
})


# sei map -----------------------------------------------------------------

r_c3 <- sei2c3(r_sei1)
r_c3 <- downsample(r_c3)

plot_map2(r_c3) +
  scale_fill_c3()

scale_fill_c3 <- function(...) {
  vals <- unname(c3Palette)
  names(vals) <- c("1","2","3")   # <-- keys equal to your raster values
  scale_fill_manual(
    values = vals,
    breaks = c("1","2","3"),
    labels = c("CSA","GOA","ORA"),
    na.value = "transparent",
    na.translate = FALSE,
    name = "SEI class",
    ...
  )
}
s <- stars::st_as_stars(r_c3)
nm <- names(s)
s[[nm]] <- factor(s[[nm]], levels = c(1,2,3))  # just treat 1/2/3 as categories

g <- plot_map2(s) +
  scale_fill_c3() +
  add_c3_points(df = d7, crs = crs, bb = bb)

ggsave(paste0('figures/sei/conceptual/c3_', runv, '_',
              scenario, '.png'),
       plot = g, width = width_map, height = height_map, dpi = 600)

# Q-curve functions -------------------------------------------------------

q_curves_cover2

df <- q_curves_cover2[[1]]

q_long <- map(q_curves_cover2, function(df) {
  pivot_longer(df,
               cols = -cover,
               names_to = 'region',
               values_to = 'Q')
}) 
df <- q_long[[1]]
names(q_long)
xlim_l <- list(
  'Sagebrush' = 0.3,
  'Pherb' = 0.6,
  'Aherb' =  0.8
)

map2(q_long, names(q_long), function(df, pft) {
  points <- d7
  xlim <- xlim_l[[pft]]
  df$cover <- df$cover*100 #convert to percent
  g <- ggplot(df, aes(x = cover, y = Q)) +
    geom_line(mapping = aes(color = region), alpha = 0.5,) +
    scale_color_manual(values = cols_wafwa_region,
                       name = NULL) +
    geom_point(data = d7, 
               mapping = aes(x = .data[[paste0(pft , '_cover')]],
                             y = .data[[paste0(pft , '_Q')]],
                             fill = c3),
               size = 2,
               stroke = 1,
               shape = 21
               ) + 
    scale_fill_manual(values = cols_c3, guide = 'none') +
    labs(x = 'Cover (%)',
         y = 'Q Value') +
    coord_cartesian(xlim = c(0, xlim*100))
  ggsave(
    paste0('figures/sei/conceptual/', 'q-curve_', pft, '.png'),
    g,
    width = 4, height = 2, dpi = 600
    )
})


