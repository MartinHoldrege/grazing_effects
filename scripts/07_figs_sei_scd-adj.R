# Purpose: Make figures of Q's and SEI and (SCD adjusted)
# these are weighted by area, showing the entire study area and individual
# sub-regions

# Author: Martin Holdrege
# Script Started: May 13, 2025

# params ------------------------------------------------------------------

source("src/params.R")
runv <- paste0(run, v_interp)
suffix <- paste0(runv, '_v1') # for figures
v_input <- 'v2' # version for some input files

scen_l <- list(RCP = c('RCP45', 'RCP85'),
               years = c('2070-2100', '2070-2100'))
# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
source("src/fig_params.R")
source("src/fig_functions.R")
source("src/general_functions.R")
source("src/SEI_functions.R")

theme_set(theme_custom1())

# Read in data ------------------------------------------------------------

options(readr.show_col_types = FALSE)

# fivenumber summaries (for boxplots) across regions of Q and SEI
# (SCD adjusted) rasters (i.e. these are area weighted summaries)
# created in 06_summarize_sei_scd-adj
qsei1 <- read_csv(paste0('data_processed/raster_means/', runv,  
                         '_q-sei_scd-adj_summaries_by-ecoregion.csv'))

# file created in 05_interpolated_summarize_sei_scd-adj.R
# mean SEI and % core by ecoregion and GCM
sei_byGCM1 <- read_csv(paste0('data_processed/raster_means/', runv, 
                             '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv'))

# average of climatic drivers, per gcm and region
# file output by 04_summarize_fire_drivers.R
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, 
                  '_fire-driver-means_by-ecoregion.csv'))

# expected burned area by gcm
# created in 05_fire_area.R
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                           v_input, "_", run, ".csv"))

# create in "scripts/06_summarize_sei_scd-adj.R"
# percent core and grow areas by region for low, median, and high SEI (across GCMs)
sei_pcent1 <- read_csv(paste0('data_processed/raster_means/', runv, 
                             '_sei-class-pcent_scd-adj_summaries_by-ecoregion.csv'))

# create in "scripts/06_summarize_sei_scd-adj.R"
# area in each of 9 climate caused SEI class change categories
c9_area2 <- read_csv(paste0('data_processed/raster_means/', runv, 
                            '_c9-area_scd-adj_summaries_by-ecoregion.csv'))

# map of ecoregions create in
# "scripts/04_ecoregion_map.R"
g_region1 <- readRDS("figures/ecoregions_v4_all.rds")

# prepare dataframes ------------------------------------------------------

qsei2 <- qsei1 %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE, 
                              add_newline = TRUE)) %>% 
  df_factor() %>% 
  filter_clim_extremes()

# boxplots ----------------------------------------------------------------

pfts <- c("Sagebrush", 'Pherb', 'Aherb')


# * Q and SEI -------------------------------------------------------------

plots <- map(pfts, function(pft) {
  box_fivenum1(df = filter(qsei2, group == pft),
                ylab = lab_q0,
                subtitle = pft) +
    expand_limits(y = c(0, 1))
})
n <- length(plots)
plots[[n]] <- plots[[n]] + theme(strip.text.y = element_blank())

g_sei <- box_fivenum1(df = filter(qsei2, group == 'SEI'),
                       ylab = lab_sei0,
                       subtitle = 'SEI')


g2 <- combine_grid_panels1(plots, remove_y_axis_labels = TRUE)
g3 <- g2 + g_sei + plot_layout(guides = 'collect')

png(paste0("figures/sei/q-sei_scd-adj_weighted_by-region_3pft_boxplot_", suffix, ".png"),
    width = 11, height = 10, units = 'in', res = 600)
g3&theme(legend.position = 'bottom')
dev.off()


# c9 area barcharts -------------------------------------------------------

c9_area3 <- c9_area2 %>% 
  select(-id) %>% 
  pivot_wider(values_from = c('area', 'area_perc'),
              names_from = summary) %>% 
  mutate(across(matches('area'), .fns = \(x) ifelse(is.na(x), 0, x)),
         c9 = c9_factor(c9)) %>% 
  df_factor()

pmap(scen_l, function(RCP, years) {
  
  df <- c9_area3 %>% 
    filter(RCP == !!RCP, 
           years == !! years)

  g <- ggplot(df, aes(c9, y = area_perc_median,fill = c9)) +
    base_c9_area() +
    geom_errorbar(aes(ymin = area_perc_low, ymax = area_perc_high, group = graze),
                  stat = 'identity',
                  width=.3,
                  position=position_dodge(0.9)) +
    guides(pattern = guide_legend(ncol = 1,
                                  override.aes = list(fill = "white", color = 'black',
                                                      size = 0.1))) + 
    labs(x = 'Projected change in SEI class',
         y = '% of region') +
    facet_manual_region(color_strips = TRUE,
                        legend.position.inside = c(0.45, 0.31),
                        # skipping letter used for map
                        region_letters = fig_letters[c(-4)]) 
    # make_legend_small()

  
  g_region2 <- g_region1 +
    theme(legend.position = 'none') +
    labs(subtitle = fig_letters[4])
  
  g2 <- g +
    inset_element(g_region2, 
                  left = -0.05, bottom = 0.05, right = 0.4, top = 0.4,
                  align_to = 'full')


  path <- paste0(
    "figures/sei/c9-bar_scd-adj_by-region_", RCP, "_", years, "_", suffix, ".png"
  )
  ggsave(path, g2, dpi = 600, width = 9, height = 7)
  
})





# GCM level results vs drivers --------------------------------------------

# burned area--attribution ------------------------------------------------

drivers2 <- drivers1 %>% 
  select(-id) %>% 
  rename(region = ecoregion)

sei_byGCM2 <- sei_byGCM1 %>% 
  select(-type, -group)

sei_byGCM3 <- drivers2 %>% 
  filter(type == "climate") %>% 
  select(-graze, -type) %>%  # grazing doesn't apply to climate variables
  pivot_wider(values_from = 'mean',
              names_from = variable) %>% 
  right_join(sei_byGCM2, by = join_by(region, RCP, years, GCM))

sei_byGCM4 <- drivers2 %>% 
  filter(type == "biomass") %>% 
  select(-type) %>% 
  pivot_wider(values_from = 'mean',
              names_from = variable) %>% 
  right_join(sei_byGCM3, by = join_by(region, RCP, years, GCM, graze)) %>% 
  rename(PSP = psp) %>% 
  df_factor() %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE))

driver_vars <- c("MAP", "MAT", "PSP", "Aherb", "Pherb", "Sagebrush")

sei_byGCM5 <- sei_byGCM4 %>% 
  pivot_longer(cols = all_of(driver_vars),
               values_to = 'mean_driver',
               names_to = 'driver') %>% 
  mutate(GCM = factor(GCM, levels = c('Current', names(cols_GCM1))),
         driver = factor(driver, levels = driver_vars))

sei_current <- sei_byGCM5 %>% 
  filter(GCM == 'Current')

rcp_year <- unique(sei_byGCM5$rcp_year) %>% 
  str_subset('RCP')

base_byGCM <- function(df_current) {
  list(
    geom_smooth(aes(linetype = graze), se = FALSE, color = 'gray'),
    geom_point(aes(shape = graze, color = GCM)),
    geom_point(data = df_current,
                 aes(shape = graze, color = 'Historical')),
    facet_grid(region~driver, scales = 'free', switch = 'x'),
      scale_color_manual(values = cols_GCM2,
                         name = 'GCM (or historical)'),
      scale_linetype(name = 'Grazing'),
      scale_shape(name = 'Grazing'),
    theme(strip.placement.x = 'outside')
  )
}

plots1 <- map(rcp_year, function(x) {
  sei_byGCM5 %>% 
    filter(rcp_year == x) %>% 
    ggplot(aes(mean_driver, SEI_mean)) +
    base_byGCM(df_current = sei_current) +
    labs(x = "",
         y = "mean SEI",
         subtitle = x,
         caption = 'Mean of fire probability predictor variable was calculated across
         pixels in ecoregion, for a given GCM')
   
  
})

plots1 <- map(rcp_year, function(x) {
  sei_byGCM5 %>% 
    filter(rcp_year == x) %>% 
    ggplot(aes(mean_driver, SEI_mean)) +
    base_byGCM(df_current = sei_current) +
    labs(x = "",
         y = "mean SEI",
         subtitle = x,
         caption = 'Mean of variables on x-axis calculated across
         pixels in ecoregion, for a given GCM')
  
})

plots2 <- map(rcp_year, function(x) {
  sei_byGCM5 %>% 
    filter(rcp_year == x) %>% 
    ggplot(aes(mean_driver, percent_csa)) +
    base_byGCM(df_current = sei_current) +
    labs(x = "",
         y = "% Core Sagebrush Area",
         subtitle = x,
         caption = 'Mean of variables on x-axis calculated across
         pixels in ecoregion, for a given GCM')
  
})

pdf(paste0("figures/sei/sei_scd-adj_vs_driver_by-GCM_", suffix, '.pdf'),
    width = 14, height = 10)
plots1
plots2
dev.off()
