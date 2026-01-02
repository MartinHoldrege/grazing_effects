# Purpose: Make figures of Q's and SEI and (SCD adjusted)
# these are weighted by area, showing the entire study area and individual
# sub-regions

# Author: Martin Holdrege
# Script Started: May 13, 2025

# params ------------------------------------------------------------------

source("src/params.R")
run <- opt$run
v_interp <- opt$v_interp
vr <- opt$vr
vr_name <- opt$vr_name
runv <- paste0(run, v_interp)
yr_lab <- opt$yr_lab
suffix <- paste0(runv, '_v1') # for figures
entire <- opt$entire
v_input <- paste0('v2', vr_name) # version for some input files
scen_l <- list(RCP = c('RCP45', 'RCP85'),
               years = rep(opt$years, 2))

ref_graze <- opt$ref_graze
rcps <- unique(scen_l$RCP)
years <- opt$years

# default should be true, false if just want to create output tables
create_figs <- TRUE

# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
source("src/fig_params.R")
source("src/fig_functions.R")
source("src/general_functions.R")
source("src/SEI_functions.R")
source('src/mapping_functions.R')
source('src/gw_functions.R')
source('src/distrib_functions.R')
theme_set(theme_custom1())

# Read in data ------------------------------------------------------------

options(readr.show_col_types = FALSE)

# fivenumber summaries (for boxplots) across regions of Q and SEI
# (SCD adjusted) rasters (i.e. these are area weighted summaries)
# created in 06_summarize_sei_scd-adj
# Note: in names 'summaries' was replaced w/ 'smry-gw'
# when switched to gw summaries
qsei1 <- read_csv(paste0('data_processed/raster_means/', runv, vr_name,
                         yr_lab,
                         '_q-sei_scd-adj_smry-gw_by-ecoregion.csv'))


# file created in 05_interpolated_summarize_sei_scd-adj.R
# mean SEI and % core by ecoregion and GCM
sei_byGCM1 <- read_csv(paste0('data_processed/raster_means/', 
                              runv, yr_lab, vr_name,
                             '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv'))

# average of climatic drivers, per gcm and region
# file output by 05_summarize_fire_drivers.R
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                  '_fire-driver-means_by-ecoregion.csv'))

# expected burned area by gcm
# created in 06_fire_area.R
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                           v_input, "_", run, ".csv"))

# map of ecoregions create in
# "scripts/04_ecoregion_map.R"
g_region1 <- readRDS(paste0("figures/ecoregions_", v_interp, "_", vr, ".rds"))

# files for histograms of change
# create in  "scripts/06_summarize_sei_scd-adj.R"
hist_sei_cref1 <- readRDS(paste0('data_processed/raster_means/', runv, '_', vr,
       '_', years,  '_sei_hist-data_by-gcm-ecoregion_cref.rds'))

hist_sei_gref1 <- readRDS(paste0('data_processed/raster_means/', runv, '_', vr,
                                 '_', years,  
                                 '_sei_hist-data_by-gcm-ecoregion_gref.rds'))

# created in 'scripts/06_summarize_sei_scd-adj.R'
c12_area1 <- read_csv(paste0('data_processed/raster_means/', 
                             runv, '_', vr, '_', years,
                             '_c12-area_by-region-smry-gw.csv'))

# prepare dataframes ------------------------------------------------------
nr <- lu(qsei1$region) # number of regions
mr <- nr/5 # multiplier for number of regions (i.e. for figure heights etc. )

qsei2 <- qsei1 %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE, 
                              add_newline = TRUE)) %>% 
  df_factor() %>% 
  filter_clim_extremes(years = years)

c12_area2 <- c12_area1 %>% 
  filter(!is.na(area)) %>%  # no data for entire study area seperately
  group_by(RCP, years, graze, summary, region) %>% 
  mutate(area_perc = area/sum(area)*100,
         c12 = c12_factor(c12),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
  df_factor()

c12_area3 <- c12_area2 %>% 
  # just groups that are CSA or GOA under the given scenario
  # i.e not stable ora or where it became ora
  filter(!str_detect(c12, 'ORA(?!.*becomes)'))
# * gcm-wise calculations -------------------------------------------------


c3_gcm2 <- c3_to_long(sei_byGCM1)

sei_pcent1 <- summarize_c3_gw(sei_byGCM1) %>% 
  filter_clim_extremes(years = years) 

# * SEI class area --------------------------------------------------------

# change in area and % asei_pcent# change in area and % area in each sei class,
# relative to historical climate (within grazing level comparisons)
c3_clim_delta <- sei_pcent1 %>% 
  filter(RCP == 'Current') %>% 
  select(-RCP, -years, -summary,-c3_percent, -rcp_year) %>% 
  right_join(
    filter(sei_pcent1, RCP != 'Current'),
    by = c('region', 'graze', 'c3', 'run'),
    suffix = c('_cur', '')) %>% 
  select(-c3_percent) %>% 
  mutate(delta_area = c3_area - c3_area_cur,
         # % change in area each c3 class
         delta_area_perc = delta_area/c3_area_cur*100)

# * change in area of each sei class relative to reference grazing level
c3_gcm2 <- sei_byGCM1 %>% 
  mutate(percent_ora = 100 - percent_csa - percent_goa) %>% 
  pivot_longer(cols = matches('percent_'),
               names_to = 'c3',
               values_to = 'c3_percent') %>% 
  mutate(c3 = str_to_upper(str_extract(c3, '[a-z]{3}$')),
         c3_area = c3_percent * area) %>% 
  select(-matches('SEI_p'), -area, -group, -type, -id, -c3_percent)

c3_gcm3 <- c3_gcm2 %>% 
  filter(graze == ref_graze) %>% 
  select(-graze) %>% 
  right_join(filter(c3_gcm2, graze != ref_graze),
         by = c('region', 'c3', 'run', 'RCP', 'years', 'GCM'),
         suffix = c('_ref', '')) %>% 
  # across grazing delta
  mutate(delta_area = c3_area - c3_area_ref,
         delta_area_perc = (delta_area)/c3_area_ref*100,
         # using delta continuous SEI for selecting the 'median' 
         # GCM to keep the methods consistent
         delta_SEI = SEI_mean - SEI_mean_ref) 

# summarize across gcms
# gcm-wise summaries b/ comparing within a climate scenario
c3_graze_delta0 <- c3_gcm3 %>% 
  # here grouping by c3 isn't really doing anything (mean SEI )
  # is just repeated
  group_by(region, run, RCP, years, c3, graze) %>% 
  summarise(across(.cols = c('delta_SEI'),
                   .fns = list(
                     'median' = median,
                     'low' = calc_low,
                     'high' = calc_high
                   )),
            .groups = 'drop') %>% 
  pivot_longer(cols = matches('delta_SEI'),
               values_to = 'delta_SEI',
               names_to = 'summary',
               names_prefix = 'delta_SEI_') %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE, 
                              add_newline = TRUE))

# joining in the absolute values that correspond to the % change summaries
# across GCMs [this approach means the summary applies to % change
# the other numbers (absolute area change) are secondary]
c3_graze_delta <- c3_graze_delta0 %>% 
  left_join(c3_gcm3, 
            by = join_by(region, run, RCP, years, c3, graze, delta_SEI)) %>% 
  df_factor()

# GCMs associated with each across gcm summary level for grazing effects
# (for use for filter). This is based on grazing effect on mean SEI
smry_gcms_sei_gref <- c3_graze_delta %>% 
  filter(c3 == 'CSA') %>% 
  select(region, run, RCP, years, rcp_year, graze, summary, GCM)

stopifnot(
  # check for for any multiple to one matches
  nrow(c3_graze_delta) == nrow(c3_graze_delta0),
  # check that all % change values are their 
          all(!is.na(filter(c3_graze_delta, RCP != 'Current'))))

# * distribution data -------------------------------------------------------

hist_sei_cref2 <- hist_dat2df(hist_sei_cref1)
hist_sei_gref2 <- hist_dat2df(hist_sei_gref1) %>% 
  filter(graze != ref_graze)

# gcm-wise summaries


hist_sei_gref3 <- hist_sei_gref2 %>% 
  right_join(smry_gcms_sei_gref, 
             by = join_by(run, RCP, years, graze, GCM, region))

hist_sei_cref3 <- sei_pcent1 %>% 
  select(RCP, years, rcp_year, graze, GCM, region, summary) %>% 
  distinct() %>% 
  # inner join b/ hist_sei_gref2 doesn't have values 
  # for current climate
  inner_join(hist_sei_cref2,
             by = join_by(RCP, years, graze, GCM, region))

# *format for output ------------------------------------------------------
# versions for saving to output as summary stats

c3_clim_delta_wide <- c3_clim_delta %>% 
  mutate(delta_area_perc = round(delta_area_perc, digits = 1),
         across(.cols = c(delta_area, c3_area, c3_area_cur),
                             .fns = round)) %>% 
  pivot_wider(id_cols = c(region, run, graze, c3, RCP, years, c3_area_cur),
              values_from = c(c3_area, delta_area, delta_area_perc),
              names_from = summary) %>% 
  arrange(RCP, years, region, graze, c3)

c3_graze_delta_wide <- c3_graze_delta %>% 
  mutate(delta_area_perc = round(delta_area_perc, digits = 1),
         across(.cols = c(delta_area, c3_area, c3_area_ref),
                .fns = round)) %>% 
  pivot_wider(id_cols = c(region, run, graze, c3, RCP, years, rcp_year),
              values_from = c(c3_area, delta_area, delta_area_perc),
              names_from = summary) %>% 
  arrange(RCP, years, region, graze, c3)

# *save tables ------------------------------------------------------------

# gw summaries across GCMs
p1 <- paste0('data_processed/area/c3/c3_clim_delta-area_gw_', vr, '_', years, '_',
            runv, '.csv')
write_csv(c3_clim_delta_wide, p1)

# gcm-wise (gw)
p2 <- paste0('data_processed/area/c3/c3_graze_delta-area_gw_', vr, '_', years, '_',
             runv, '.csv')
write_csv(c3_graze_delta_wide, p2)

# boxplots ----------------------------------------------------------------
if(create_figs) {
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
                       subtitle = 'SEI') +
  theme(legend.position = 'none')


g2 <- combine_grid_panels1(plots, remove_y_axis_labels = TRUE)
g3 <- g2 + g_sei + plot_layout(guides = 'collect')

png(paste0("figures/sei/q-sei_scd-adj_weighted_by-region_3pft_boxplot_gw", 
           vr_name, yr_lab, "_", suffix, ".png"),
    width = 11, height = 2.3 + 7.7*mr, units = 'in', res = 600)
print(g3&theme(legend.position = 'bottom'))
dev.off()


# c3 area -----------------------------------------------------------------

base_bar <- function() {
  list(geom_bar(aes(x = graze, 
                    fill = c3), width = 0.75,
                stat = 'identity', position = position_dodge2(padding = 0)),
       adjust_graze_x_axis(),
       scale_fill_manual(values = cols_c3,
                         name = lab_c3)
  )
}

hline <- function() {
  geom_hline(yintercept = 0, linetype = 1,
             alpha = 0.5)
}

# dimension of 10 panel figures (2 cols)
width10 <- 5
height10 <- 1.5 + 7.5*mr

# dimensions of 5 panel figures (1 per region)
widthr <- if (nr <= 6) 6.3 else 7
heightr <- if (nr <= 6) 6 else 7

map(rcps, function(rcp) {
  suffix <- paste0(vr_name, '_', rcp, yr_lab, "_", runv)
  df <- sei_pcent1 %>% 
    filter(RCP %in% c('Current', rcp),
           summary == 'median')
  # for labeling
  rcp_year <- unique(df$rcp_year[df$RCP == rcp])
  
  # absolute area
  g <- ggplot(df, aes(y = c3_area)) +
    base_bar() +
    facet_grid(region ~rcp_year, scales = 'free_y') +
    labs(y = lab_area0)
  
  ggsave(paste0('figures/sei/c3_area/c3-bar_abs-area_by-region_gw_',suffix, '.png'),
         plot = g,
         dpi = 600, height = height10, width = width10)
  
  # change in area relative to historical climate
  # (within grazing level)
  df <- c3_clim_delta %>% 
    filter(RCP %in% c('Current', rcp),
           summary == 'median') 
  
  g <- ggplot(df, aes(y = delta_area)) +
    hline()+
    base_bar() +
    facet_manual_region(v = vr,
                        legend.position.inside = c(0.1, 0.2),
                        scales = 'free_y') +
    labs(y = 'Change in area relative to historical climate (ha)',
         subtitle = str_replace(rcp_year, '\n', ' '))
  
  ggsave(paste0('figures/sei/c3_area/c3-bar_clim-delta-area_by-region_gw_', 
                suffix, '.png'),
         plot = g,
         dpi = 600, height = heightr, width = widthr)
  
  g <- ggplot(df, aes(y = delta_area_perc)) +
    hline()+
    base_bar() +
    facet_manual_region(legend.position.inside = c(0.1, 0.2),
                        v = vr) +
    labs(y = 'Change in SEI class area relative to historical climate (%)',
         subtitle = str_replace(rcp_year, '\n', ' '))
  
  ggsave(paste0('figures/sei/c3_area/c3-bar_clim-delta-area-perc_by-region_gw_', 
                suffix, '.png'),
         plot = g,
         dpi = 600, height = heightr, width = widthr)
  
  # change in c3 area relative to moderate grazing 
  # within a climate scenario
  df <- c3_graze_delta_wide %>% 
    filter(RCP %in% c('Current', rcp)) 
  
  # ggsave(paste0('figures/sei/c3_area/c3-bar_graze-delta-area_by-region',
  #               suffix, '.png'),
  #        plot = g,
  #        dpi = 600, height = height10, width = width10)

  g <- ggplot(df, aes(y = delta_area_perc_median)) +
    hline()+
    base_bar() +
    facet_grid(region ~rcp_year) +
    labs(y = paste('Change in SEI class area relative to',
                   str_to_lower(ref_graze), 'grazing (%)'))

  ggsave(paste0('figures/sei/c3_area/c3-bar_graze-delta-area-perc-gw_by-region',
                suffix, '.png'),
         plot = g,
         dpi = 600, height = height10, width = width10)
  
  g <- df %>% 
    filter(region == levels(region)[1]) %>% 
    ggplot(aes(y = delta_area_perc_median)) +
    hline()+
    base_bar() +
    facet_grid(~rcp_year) +
    labs(y = paste('Change in SEI class area relative to',
                   str_to_lower(ref_graze), 'grazing (%)'))
  
  ggsave(paste0('figures/sei/c3_area/c3-bar_graze-delta-area-perc-gw_entire_', 
                rcp, yr_lab, "_", runv, '.png'),
         plot = g,
         dpi = 600, height = 5, width = 7)

})

# *stacked bar ------------------------------------------------------------

# dimensions 
widthr  <- if (nr <= 6) 6.3 else 7
heightr <- if (nr <= 6) 6 else 7

# build & save for each target RCP
df <- sei_pcent1 %>% 
  filter(region != entire,
         c3 != 'ORA') %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE))
purrr::walk(rcps, function(rcp) {
  suffix2 <- paste0(vr, "_", rcp, '_', years,   "_", runv)

  g_stack <- make_stacked_c3_panels(df = df, rcp = rcp, vr = vr)

  ggsave(
    filename = paste0("figures/sei/c3_area/c3-bar-stack_abs-area_by-region_gw_",
                      suffix2, ".pdf"),
    plot = g_stack,
    height = heightr, width = widthr,
    device = cairo_pdf, # so text can be edited
    family = "sans"  # or a specific font like "Arial"
  )

})

# c12 stacked bar ---------------------------------------------------------
# 

purrr::walk(rcps, function(rcp) {

  g_stack <- make_stacked_c12_panels(df = c12_area3, rcp = rcp, vr = vr)
  
  ggsave(
    filename = paste0("figures/sei/c3_area/c12-bar-stack_abs-area_by-region_gw_",
                      vr, "_", rcp, '_', years,   "_", runv, ".pdf"),
    plot = g_stack,
    height = 4.5, width = 4.7,
    device = cairo_pdf, # so text can be edited
    family = "sans"  # or a specific font like "Arial"
  )
  
})
# delta SEI distribution --------------------------------------------------
# histograms of change in SEI

hist_base <- function() {
  list(
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5),
    geom_bar(aes(x = mids, y = density), stat = 'identity'),
    coord_cartesian(xlim = c(-0.3, 0.3)),
    labs(x = lab_sei1)
  )
}

for(rcp in rcps) {
  # grazing effect
  g <- hist_sei_gref3 %>% 
    filter(RCP %in% c('Current', rcp),
           summary == 'median') %>% 
    ggplot() +
    hist_base() +
    facet_grid(region~rcp_year+graze)
  
  filename <- paste0('figures/sei/distributions/hist_sei_gref_by-region_',
                     vr, '_', rcp, '_', years, '_', runv, '.png')
  
  ggsave(filename, 
      plot = g,
      width = 10,
      height = 8,
      dpi = 600)
  # climate effect
  df_cref <- hist_sei_cref3 %>% 
    filter(RCP == rcp, summary == 'median') 
  g <- ggplot(df_cref) +
    hist_base() +
    facet_grid(region~graze) +
    labs(subtitle = str_replace(unique(df_cref$rcp_year), '\n', ' '))
    
  filename <- paste0('figures/sei/distributions/hist_sei_cref_by-region_',
                     vr, '_', rcp, '_', years, '_', runv, '.png')
  ggsave(filename, 
         plot = g,
         width = 7,
         height = 8,
         dpi = 600)
}




# c9 area barcharts -------------------------------------------------------

# see pre-Nov 2025 commits
#   path <- paste0(
#     "figures/sei/c9-bar_scd-adj_by-region", vr_name, 
#     "_", RCP, "_", years, "_", suffix, ".png"
#   )
#   ggsave(path, g2, dpi = 600, width = 9*mr, height = 7*mr)
#   
# })




# GCM level results vs drivers --------------------------------------------

# --attribution ------------------------------------------------

drivers2 <- drivers1 

sei_byGCM2 <- sei_byGCM1 %>% 
  select(-type, -group)

sei_byGCM3 <- drivers2 %>% 
  filter(variable %in% c('MAP', 'MAT', 'psp', 'PSP')) %>% 
  select(-graze, -matches('p\\d+'), -matches('median')) %>%  # grazing doesn't apply to climate variables
  pivot_wider(values_from = 'mean',
              names_from = variable) %>% 
  right_join(sei_byGCM2, by = join_by(region, RCP, years, GCM))

sei_byGCM4 <- drivers2 %>% 
  filter(variable %in% c('Aherb', 'Pherb', 'Sagebrush')) %>% # biomass variables 
  select(-matches('p\\d+'), -matches('median')) %>%  # grazing doesn't apply to climate variables
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

pdf(paste0("figures/sei/sei_scd-adj_vs_driver_by-GCM_",vr, "_", years, "_", 
           suffix, '.pdf'),
    width = 14, height = 10*mr)
print(plots1)
print(plots2)
dev.off()

} # end create figs


