# Purpose: Select sites for doing stepwat simulations, with the goal of 
# representing the whole study area well (1 km grid cells containing >50% SCD 30 m cells)

# Author: Martin Holdrege

# Started: June 17, 2024

# parameters --------------------------------------------------------------

scd_prop <- 0.5 # proportion of grid cell that contains scd cells, to be part of study area
sage_gt <- 1 #minimum % cover of sage, for SCD pixel to be counted. 
period <- '1991-2020' # period from which climate normals were calculated

# criteria used in Renne et al 2024 
criteria <- c(bio1 = 1.55, 
                 bio4 = 61.1, 
                 bio9 = 3.53, 
                 bio12 = 84.9, 
                 bio15 = 9.21, 
                 bio18 = 33.8
)
crit_descript <- 'Renne et al matching criteria'
v <- 'v2'
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(patchwork)
theme_set(theme_bw())
library(rMultivariateMatching)
source('src/mapping_functions.R')
source("src/matching_functions.R")

# load data ---------------------------------------------------------------

template <- rast("data_processed/interpolation_data/cellnumbers.tif")

# daymet data
r_day1 <- rast(paste0("data_raw/daymet_v4_monthly_normals_", period, ".tif"))
day1 <- read_csv(paste0("data_processed/interpolation_data/clim_for_interp_daymet-v4_",
                        period, ".csv"),
                 col_types = paste0(rep('n', 26), collapse = ""))

# * grid-met --------------------------------------------------------------

buffer <- 0
# file created in 00_gridmet_normals.R
r_grid1 <- rast(paste0("data_processed/interpolation_data/gridmet_bioclim-vars_19912020_",
                     buffer, "buffer.tif"))

# preparing dataframes ----------------------------------------------------

grid1 <- as.data.frame(r_grid1)
grid1$cellnumber <- row.names(grid1)
sage_var <- paste0("prop_sage_gt_", sage_gt)
day2 <- day1 %>% 
  filter(.data[[sage_var]] > scd_prop) # filter to mask

# normalize variables
grid_norm <- apply_criteria(grid1, criteria = criteria)
day_norm <- apply_criteria(day2, criteria = criteria)

# kpoints -----------------------------------------------------------------
# finding optimum number of sites 

# columns used for evenly spacing initial 'centers' for kmeans
# MAT, MAP, ppt warmest quarter
cols_for_centers <- c("bio1", "bio12", "bio18")

# test <- kpoints2(target = day_norm,
#          potential_sites = grid_norm,
#          criteria_vars = names(criteria),
#          digits = 4,
#          k = 50,
#          verbose = TRUE)

kseq <- seq(25, 600, by = 25)

set.seed(123)
kresults <- map(kseq, \(k) {
  print(k)

  kpoints2(target = day_norm, 
    potential_sites = grid_norm,
    criteria_vars = names(criteria),
    cols_for_centers = cols_for_centers,
    digits = 4,
    k = k,
    verbose = TRUE)
  })

# saving b/ previous step is slow
saveRDS(kresults, 
        paste0('data_processed/site_selection/kpoints_', v, '.rds'))

kresults <- readRDS(paste0('data_processed/site_selection/kpoints_', v, '.rds'))

# * existing coverage -----------------------------------------------------

cuttoff <- 1

cov200 <- prop_coverage(day_norm,
              sites = day_norm[!is.na(day_norm$site_id), ],
              criteria_vars = names(criteria),
              cuttoff = cuttoff)

cov200 <- tibble(coverage = cov200,
                 k = 200)

# * kpoints figs ----------------------------------------------------------

cap1 <- paste0('sites selected from gridmet (', buffer, ' buffer).\n',
               'Study area = 1km pixels with >', scd_prop*100, 
               '% SCD pixels that have >', sage_gt, '% sage.\n',
               crit_descript)
cov1 <- map(kresults, \(l) {
  tibble(coverage = l$coverage,
         converged = l$ifault == 0,
         k = nrow(l$sites))
}) %>% 
  bind_rows()


g <- ggplot(cov1, aes(k, coverage)) +
  geom_point(aes(color = converged)) +
  geom_point(data = cov200, aes(color = "200 published sites")) +
  labs(y = 'Proportion of daymet cells covered by k sites',
       caption = cap1) +
  scale_color_manual(values = c('blue', 'gray', 'black'),
                     name = 'kmeans converged')

ggsave(paste0("figures/site_selection/kpoints_vs_coverage_", v, '.jpg'),
       g)

# * climate envelope and points ----------------------------------------------
k_to_plot <- c(200, 300, 400) # select 3 different site numbers to plot below
# figures showing study are climate envelope along with the locations of the 200
# points (to see how representative they are)

day3 <- day2 %>% 
  rename(MAT = bio1, MAP = bio12)

# cellnumbers of selected sites

sites_grid_l <- map(k_to_plot, function(k) {
  sites_norm_grid <- keep(kresults, \(x) x$k == k)[[1]]$sites
  
  out <- list()
  out$df <- grid1 %>% 
    filter(cellnumber %in% sites_norm_grid$cellnumber) %>% 
    rename(MAT = bio1, MAP = bio12)
  out$k <- k
  out
}) 

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


g <- ggplot(day3, aes(MAT, MAP)) +
  base_density() +
  coord_cartesian(xlim = range(day3$MAT, na.rm = TRUE),
                  ylim = range(day3$MAP, na.rm = TRUE))

g_list <- map(sites_grid_l, \(l) {
  g + geom_point(data = l$df, alpha = 1, color = 'blue', size = 0.5) +
    labs(subtitle = paste(l$k, "gridmet 'sites' selected via kpoints2"))
})

g1 <- g + geom_point(data = day3[!is.na(day3$site_id), ], 
                     alpha = 1, color = 'blue', size = 0.5)+
  labs(subtitle = '200 sites from Palmquist et al (daymet)')

comb <- wrap_plots(c(list(g1), g_list)) +
  plot_layout(ncol = 2) +
  patchwork::plot_annotation(caption = cap1)

pdf(paste0('figures/site_selection/selected-sites_climate-envelope_', v, '.pdf')
    , height = 10, width = 9)
comb
dev.off()

