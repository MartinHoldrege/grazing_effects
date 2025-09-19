# Purpose: Create figures of expected burned area and area in different
# fire return inverval categories

# Author: Martin Holdrege

# Started: March 25, 2025


# params ------------------------------------------------------------------

source("src/params.R")

run <- opt$run
v_interp <- opt$v_interp
vr_name <- opt$vr_name
yr_lab <- opt$yr_lab
v <- 'v2' # version of input files (and for now also used in output file names)
suffix <- paste0(v, vr_name, "_", run)
runv <- paste0(run, v_interp)
# create additional figures to explore how age groups/fire probability relate
# (not using actual data)
explanatory_figures <- TRUE 
rcps <- c('RCP45', 'RCP85') # figures seperately made for both RCPs
years <- opt$years
entire <- 'Entire study area' # name of factor level for entire study
create_figs <- TRUE # if false just create intermediate csv output
# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(patchwork)
source("src/general_functions.R")
source("src/SEI_functions.R")
source("src/fig_functions.R")
source("src/fig_params.R")
source("src/probability_functions.R")
theme_set(theme_custom1())
options(readr.show_col_types = FALSE)

# read in data ------------------------------------------------------------

# created in the 06_fire_area.R script
# expected burned area
# includes mid and end century
ba3a <- read_csv(
  paste0("data_processed/area/expected-burn-area_smry_", suffix,".csv"))

# expected burned area by gcm
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                            suffix, ".csv"))

# burned area by current SEI class and ecoregion (pixelwise summaries across gcms) 
ba_c3_smry1 <- read_csv(
  paste0("data_processed/area/expected-burn-area_by-c3-smry_", suffix, ".csv"))

# burned area by current SEI class  and ecoregion (values for each GCM)
ba_c3_gcm1 <- read_csv(
  paste0("data_processed/area/expected-burn-area_by-c3-GCM_", suffix, ".csv"))

# expected area in different years since fire age groups
# (calculated gcm-wise)
area_age_group3 <- read_csv(
  paste0("data_processed/area/area-by-age-group_", suffix, ".csv"))

# expected area in different years since fire age groups
# (calculated pixel-wise)
area_age_group3_pw <- read_csv(
  paste0("data_processed/area/area-by-age-group_pw_", suffix, ".csv"))

# area of our study region in each of 3 ecoregions
area_eco <- read_csv(paste0("data_processed/area/ecoregion-area_", v, vr_name,
                            ".csv"))

# created in the 04_summarize_fire_drivers.R script
# means, by ecoregion, of the drivers of fire probability
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                            '_fire-driver-means_by-ecoregion.csv'))

# create in "scripts/02b_fire_attribution.R"
# site level data on change in fire prob due to individual pred vars
# and dominant drivers of change
driver_l <- readRDS('data_processed/site_means/fire_dominant_drivers.RDS')

# created in 07_bio_figs_weighted.R
sei2 <- readRDS(paste0('data_processed/temp_rds/sei_df', vr_name, yr_lab, '.rds'))

# created in "scripts/06_summarize_sei_scd-adj.R"
sei_pcent1 <- read_csv(paste0('data_processed/raster_means/', 
                              runv, vr_name, yr_lab, 
                              '_sei-class-pcent_scd-adj_summaries_by-ecoregion.csv'))

# created in "scripts/05_interpolated_summarize_sei_scd-adj.R"
sei_pcent_gcm1 <- read_csv(paste0('data_processed/raster_means/', 
                                  runv, yr_lab, vr_name,   
                                  '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv'))

# pixelwise and gcm wise mean sei by ecoregion and current sei class
c3eco_sei2_smry <- read_csv(paste0('data_processed/raster_means/',
                                   runv, yr_lab,vr_name,
                                  '_sei-mean_scd-adj_smry-by-region-c3.csv'))

c3eco_sei2_gcm <- read_csv(paste0('data_processed/raster_means/', 
                                  runv, yr_lab, vr_name,
                                 '_sei-mean_scd-adj_by-GCM-region-c3.csv'))

# interpolation weights from:
# "scripts/04_interpolation_weights.R"
weights1 <- read_csv(paste0('data_processed/interpolation_data/interpolation_weights_', 
                            v_interp, vr_name, '.csv')) %>% 
  df_factor()

# vectors -----------------------------------------------------------------

ecoregions <- region_factor(area_eco$ecoregion) %>% 
  levels() %>% 
  setNames(nm = .)
nr <- length(ecoregions)
mr <- nr/5 # multiplier for figure heights
age_groups <- create_age_groups()

# fig params --------------------------------------------------------------

line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot
line_loc2 <- 1:2 + 0.5

if(vr == 'r1.0') {
  height_3p <- 6 # width for 3 panel figures
  width_3p <- 5 # height for 3 panel figures
} else if (vr == 'r1.1') {
  height_3p <- 10 # width for 3 panel figures
  width_3p <- 5 # height for 3 panel figures
} else {
  stop('update needed')
}

# prepare dataframes ------------------------------------------------------

ba3 <- df_factor(ba3a) %>% 
  filter_clim_extremes(years = years) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
  mutate(across(matches("^area_"), .fns = \(x) x/area_total*100,
                .names = '{.col}_perc')) 

area_age_group3 <- area_age_group3 %>% 
  filter_clim_extremes(years = years) %>% 
  left_join(area_eco, by = join_by(ecoregion)) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         # not reproducible (yet)
         age_group = factor(age_group, levels = names(age_groups)),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
  df_factor()
  
area_age_group4_pw <- area_age_group3_pw %>% 
  filter_clim_extremes(years = years) %>% 
  left_join(area_eco, by = join_by(ecoregion)) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(age_group = factor(age_group, levels = rev(names(age_groups))),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE),
         area_median_perc = area_median/area*100) %>% 
  df_factor()
  
# *trade-offs ba vs core --------------------------------------------------

# showing area of core vs expected burned area to illustrate the grazing
# trade-off
sei_pcent2 <- sei_pcent1 %>% 
  select(-id, -group, -type) %>% 
  pivot_wider(values_from = c('percent_csa', 'percent_goa'),
              names_from = 'summary') %>% 
  df_factor() %>% 
  # correct run is in file name when read in
  # but run column is truncated (not including the yr month,
  # potentially b/ of a renaming issue in "05_interpolated_summarize_sei_scd-adj.R)
  # for that reason not using run to join
  select(-run)

sei_pcent3 <- ba3 %>% 
  select(-rcp_year) %>% 
  rename_with(.fn = \(x) paste0('ba_', x), .cols = matches('^area')) %>% 
  left_join(rename(area_eco, total_area = area), by = 'ecoregion') %>% 
  # adding percent of total area
  mutate(across(matches('area_'), .fns = \(x) x/total_area*100, 
                .names = '{.col}_perc'),
         region = region_factor(ecoregion)) %>% 
  select(-ecoregion) %>% 
  left_join(sei_pcent2, by = join_by(RCP, years, graze, region)) %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE),
         percent_csagoa_median = percent_goa_median + percent_csa_median,
         percent_csagoa_low = percent_goa_low + percent_csa_low,
         percent_csagoa_high = percent_goa_high + percent_csa_high)
  
# ** gcm level ------------------------------------------------------------

sei_pcent_gcm2 <- sei_pcent_gcm1 %>% 
  mutate(run = str_replace(run, !!runv, !!run)) %>% 
  select(-type, -id, -group)

sei_pcent_gcm3 <- ba_gcm1 %>% 
  rename(ba_area = area) %>% 
  select(-type, -run2) %>% 
  filter_clim_extremes(years = years) %>% 
  left_join(rename(area_eco, total_area = area), by = 'ecoregion') %>% 
  mutate(ba_area_perc = ba_area/total_area*100) %>% 
  rename(region = ecoregion) %>% 
  left_join(sei_pcent_gcm2, 
            by = join_by(region, run, RCP, years, graze, GCM)) %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE),
         rcp_year_GCM = paste(rcp_year, GCM),
         percent_csagoa = percent_goa + percent_csa) %>% 
  df_factor() %>% 
  arrange(rcp_year, GCM, graze) # arranging for geom_path


# ** by SEI class ---------------------------------------------------------

# pixelwise summaries across gcms
c3eco_smry1 <- c3eco_sei2_smry %>% 
  select(-run, -id, -run2, -type, -group) %>% 
  left_join(ba_c3_smry1, 
            join_by(c3, region, RCP, years, graze, summary))  %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE),
         rcp_year_c3 = paste0(rcp_year, c3)) %>% 
  filter_clim_extremes(years = years) %>% 
  df_factor() %>% 
  select(-type) %>% 
  arrange(rcp_year, c3, graze) # order matters for geom_path


test <- c3eco_smry1$area.x - c3eco_smry1$area.y
stopifnot(abs(test) < 1)# areas calculated in the two datasets need to ~ match

c3eco_smry2 <- c3eco_smry1 %>% 
  select(-area.y, -run2) %>% 
  rename(area = area.x) %>% 
  mutate(expected_ba_perc = expected_ba/area*100)

# calculating values for entire region
c3eco_smry3 <- c3eco_smry2 %>% 
  group_by(c3, RCP, years, graze, summary, run, rcp_year, rcp_year_c3) %>% 
  summarize(across(c(SEI_mean, expected_ba, expected_ba_perc),
                   .fns = \(x)weighted.mean(x, w = area)),
            area = sum(area),
            .groups = 'drop') %>% 
  mutate(region = entire) %>% 
  bind_rows(c3eco_smry2) %>% 
  df_factor()


# values by GCM
c3eco_gcm1 <- c3eco_sei2_gcm %>% 
  select(-run, -id, -run2, -type, -group) %>% 
  left_join(ba_c3_gcm1, 
            join_by(c3, region, RCP, years, graze, GCM))  %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE),
         rcp_year_c3_gcm = paste0(rcp_year, c3, GCM)) %>% 
  filter_clim_extremes(years = years) %>% 
  df_factor() %>% 
  select(-type) %>% 
  arrange(rcp_year, c3, graze) # order matters for geom_path

c3eco_gcm2 <- c3eco_gcm1 %>% 
  select(-area.y, -c3eco, -run2, -id_noGCM) %>% 
  rename(area = area.x) %>% 
  mutate(expected_ba_perc = expected_ba/area*100)

test <- c3eco_gcm1$area.x - c3eco_gcm1$area.y
stopifnot(abs(test) < 1)# areas calculated in the two datasets need to ~ matc

# calculate for entire study area
c3eco_gcm3 <- c3eco_gcm2 %>% 
  group_by(c3, RCP, years, graze, GCM, run, rcp_year, rcp_year_c3_gcm) %>% 
  summarize(across(c(SEI_mean, expected_ba, expected_ba_perc),
                   .fns = \(x)weighted.mean(x, w = area)),
            area = sum(area),
            .groups = 'drop') %>% 
  mutate(region = entire) %>% 
  bind_rows(c3eco_gcm2) %>% 
  df_factor()

# * drivers of fire change ------------------------------------------------

dom_drivers2 <- driver_l$dominant_driver1 %>% 
  filter_clim_extremes(years = years) %>% 
  select(-id) %>% 
  rename(dom_driver = dominant_driver) %>% 
  filter(run == !!run) %>% 
  full_join(weights1, by = 'site',
            relationship = 'many-to-many') %>% 
  mutate(area_dom = weight2area(weight),
         dom_driver = driver2factor(dom_driver, include_none = TRUE)) %>% 
  group_by(region, run, RCP, years, graze, dom_driver) %>%
  summarize(area_dom = sum(area_dom),
            .groups = 'drop_last') %>% 
    mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE),
         perc_area = area_dom/sum(area_dom)*100)  %>% 
  ungroup() %>% 
  complete(dom_driver, region, run, nesting(RCP, years, graze, rcp_year),
           fill = list(area_dom = 0, perc_area = 0)) %>% 
  df_factor()

# * fire change vs driver change ------------------------------------------

# change in fire probability due to one predictor variable,
# and the change in that predictor, by site
one_change_smry1 <- driver_l$one_change_smry1 # median across GCMs (site-wise)
one_change_gcm1 <- driver_l$one_change_gcm1 # values by GCM

one_change_prep <- function(df) {
  df %>% 
    left_join(weights1, by = 'site', relationship = 'many-to-many') %>% 
    rename(pred_var = pred_var_cur) %>% 
    filter_clim_extremes(years = years) %>% 
    filter(run == !!run) %>% 
    mutate(pred_var = str_replace(pred_var, 'psp', 'PSP'),
           pred_var = driver2factor(pred_var))
}

one_change_smry2 <- one_change_smry1 %>% 
  one_change_prep()

one_change_gcm2 <- one_change_gcm1 %>% 
  one_change_prep()

# summarized across sites (sas = summarized across sites)
create_sas <- function(df) {
  df %>% 
    summarize(across(c(delta_pred_value, delta_1var),
                   .fns = \(x) weighted.mean(x, w = weight)),
            .groups = 'drop') %>% 
    mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
    arrange(run, region, years, RCP, pred_var, graze)
}


# "pixel level" medians (across gcms) then summarized across sites/pixels
one_change_smry_sas1 <- one_change_smry2 %>% 
  group_by(run, years, RCP, graze, region, pred_var) %>% 
  create_sas() 

one_change_gcm_sas1 <- one_change_gcm2 %>% 
  group_by(run, years, RCP, graze, GCM, region, pred_var) %>% 
  arrange(GCM) %>% 
  create_sas()%>% 
  mutate(rcp_year_gcm = paste(rcp_year, GCM))


# * write files -------------------------------------------------------------

write_csv(c3eco_gcm3, paste0('data_processed/raster_means/', 
                             runv, '_', vr, "_", years, 
       '_sei-mean_ba_by-GCM-region-c3.csv'))

if(create_figs){
# expected burned area figs --------------------------------------------------

plots <- map(ecoregions, function(region) {
  total_area <- area_eco$area[area_eco$ecoregion == region]
  ba3 %>% 
    filter(ecoregion == region) %>% 
    ggplot(aes(x = rcp_year, y = area_median)) +
    geom_errorbar(aes(ymin = area_low, ymax = area_high, group = graze), 
                  position = position_dodge(width = 0.5), width = 0) +
    geom_point(aes(y = area_median, color = graze), 
               position = position_dodge(width = 0.5)) +
    scale_y_continuous(sec.axis = sec_axis(transform = \(x) x/total_area*100,
                                           name = '% of region')) +
    scale_color_manual(values = cols_graze, name = 'Grazing') +
    geom_vline(xintercept = line_loc2, linetype = 2) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 25, hjust = 1, size = rel(0.8))) +
    labs(x = NULL,
         y = lab_ba0,
         subtitle = region)
})

plots <- remove_y_titles(plots, index_keep_y = 3, 
                         index_keep_y_sec = 4)

g <- patchwork::wrap_plots(plots) +
  plot_layout(guides = 'collect', axes = 'collect', ncol = 2)

g2 <- g&theme(legend.position = 'bottom')

jpeg(paste0("figures/fire/area/expected_ba_dotplot_", 
            v, vr_name, yr_lab, '_', run, '.jpg'),
     units = 'in', width = width_3p, height = height_3p, res = 600)
g2
dev.off()


# age_group figs -----------------------------------------------------


# * grouped by climate scenario -------------------------------------------

# figure of expected area in a age class  
fig_area_by_scenario_group <- function(df, total_area) {
  ggplot(df, aes(x = rcp_year))+
    geom_errorbar(aes(ymin = area_low, ymax = area_high, group = graze), 
                  position = position_dodge(width = 0.5), width = 0) +
    geom_point(aes(y = area_median, color = graze), 
               position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = line_loc2, linetype = 2) +
    scale_color_manual(values = cols_graze, name = 'Grazing') +
    facet_wrap(~age_group, nrow = 1) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                     size = rel(0.8)),
          plot.subtitle = element_text(size = rel(0.7),
                                      margin  = margin(4, 0, 2, 0)),
          strip.text = element_text(size = rel(0.7),
                                    margin = margin(0, 0, 0, 0)),
          plot.margin = margin(0, 0, 0, 0)) +
    labs(x = NULL,
         y = 'Expected area in age class (ha/yr)') +
    scale_y_continuous(sec.axis = sec_axis(transform = \(x) x/total_area*100,
                                           name = '% of ecoregion')) +
    expand_limits(y = 0)
}

old_group <- levels(area_age_group3$age_group)[length(age_groups)]


plots0 <- map(ecoregions, function(region) {
  total_area <- area_eco$area[area_eco$ecoregion == region]
  
  df <- area_age_group3  %>% 
    filter(ecoregion == region) 
  
  # only keep axis text in bottom row
  tmp_theme <- if(region == ecoregions[[length(ecoregions)]]) {
    theme()
  } else {
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  }
  g1 <- df %>% 
    filter(age_group != old_group) %>% 
    fig_area_by_scenario_group(total_area = total_area) + 
    labs(subtitle = region) +
    tmp_theme
  
  
  g2 <- df %>% 
    filter(age_group == old_group) %>% 
    fig_area_by_scenario_group(total_area = total_area) +
    tmp_theme
  
  list(g1, g2)
})

plots <- reduce(plots0, c)
g <- patchwork::wrap_plots(plots, ncol = 2, guides = 'collect',
                           tag_level = 'keep',byrow = TRUE, widths = c(3, 1)) +
  plot_layout(axes = 'collect')

g2 <- g&theme(legend.position = 'bottom',
              axis.title.y = element_text(),
              axis.title.y.right = element_text()
)

jpeg(paste0("figures/fire/area/area_age_group_scen-grouping_",
            v, vr_name, yr_lab, "_", run,'.jpg'),
     units = 'in', width = 6, height = 12, res = 600)
g2
dev.off()


# * stacked barchart ------------------------------------------------------

# check--area should sum to 100% (with some rounding error)
check <- area_age_group4_pw %>% 
  group_by(ecoregion, graze, rcp_year) %>% 
  summarize(sum_perc = sum(area_median_perc)) %>% 
  pull(sum_perc)

stopifnot(abs(check - 100) < 0.01)

map(rcps, function(rcp){
  g <- area_age_group4_pw %>% 
    filter(RCP %in% c('Current', rcp)) %>% 
    ggplot(aes(x = graze, y = area_median_perc, fill = age_group)) +
    geom_bar(stat = 'identity') +
    facet_grid(ecoregion~rcp_year, scales = 'free_y') +
    adjust_graze_x_axis() +
    scale_fill_manual(values = rev(cols_agegroup),
                      name = lab_agegroup) +
    labs(y = lab_areaperc0b) +
    theme(strip.text = element_text(size = rel(0.7)))
  
  filename <- paste0("figures/fire/area/area_age_group_by-region_bar_", rcp, 
                     yr_lab, "_", 
                     suffix, '.png')
  ggsave(filename, g, height = 2+6*mr, width = 4, dpi = 600)
})

  
# trade-offs ba vs core ---------------------------------------------------

# *dotplot ----------------------------------------------------------------

# see pre 8/7/2025 commits

# ** gcm level ------------------------------------------------------------

# see pre 8/7/2025 commits

# ** median and GCM level --------------------------------------------------
# figure that shows the median results, and then faint GCM lines
# b/ the low/high error bars in the x and y diretions are misleading
# (they're correlated so direction of effect isn't actually as uncertain
# as they suggest)

args <- list(
  rcp = c(rcps, rcps),
  xvar = c('percent_csa', 'percent_csa', 'percent_csagoa', 'percent_csagoa')
)

xlabs <- c(percent_csa = '% Core Sagebrush Area',
           percent_csagoa = '% CSA + GOA')

pmap(args, function(rcp, xvar) {
  df_med <- sei_pcent3 %>% 
    filter(RCP %in% c('Current', rcp))
  df_gcm <- sei_pcent_gcm3 %>% 
    filter(RCP %in% c('Current', rcp))
    
  x_med <- paste0(xvar, '_median')
  g <- ggplot(df_med, aes(.data[[x_med]], ba_area_median_perc)) +
    geom_path(data = df_gcm, aes(.data[[xvar]], ba_area_perc, 
                                         group = rcp_year_GCM, linetype = rcp_year,
                                         color = rcp_year),
              alpha = 0.5,
              linewidth = 0.2, 
              show.legend = FALSE) + # b/ different linewidth get's double plotted on legend +
    base_tradeoff(linetypes_scen = linetypes_scen,
                  xlab = xlabs[xvar]) +
    facet_manual_region(legend.position.inside = c(.1, 0.15),
                        v = vr) +
    make_legend_small()
  g
  
  ggsave_tradeoff(
    g = g,
    xvar = str_replace(xvar, "percent_", ""),
    prefix = paste0('med-GCM-fix-', rcp, yr_lab),
    height = 4.5*mr,
    width = 7*mr
  )
})


# * by current C3 ---------------------------------------------------------
# showing mean SEI vs expected burned area for each combination of ecoregion 
# and current SEI class

args <- expand_grid(rcp = rcps,
                    # also show lines for each GCM
                    gcm_path = c(TRUE, FALSE))

width_tradeoff <- 8
height_tradeoff <- 7
if(nr == 9) {
  remove_x = 1:6
  remove_y = c(2, 3, 5 , 6, 8, 9)
} else if(nr == 5) {
  remove_x = c(2, 3)
  remove_y = c(2, 3, 5)
  width_tradeoff <- 7
  height_tradeoff <- 5
} else {
  stop('update')
}

pmap(args, function(rcp, gcm_path) {
  
  df_smry <- c3eco_smry3 %>% 
    filter(RCP %in% c('Current', rcp),
           summary == 'median') %>% 
    mutate(c3_cur = relable_c3_current(c3),
           graze_long = relable_graze_long(graze))
  
  df_gcm <- c3eco_gcm3 %>% 
    filter(RCP %in% rcp)%>% 
    mutate(c3_cur = relable_c3_current(c3),
           graze_long = relable_graze_long(graze))
  
  
  # barcharts showing % of total area by SEI class
  # (to be used as insets)
  bar_l <- area_bar_map(df_smry)
  plotsl1 <- tradeoff_lines_map(df_smry = df_smry,
                                df_gcm = df_gcm,
                                v = vr,
                                gcm_path = gcm_path)

  
  plotsl1 <- remove_axis_labels(plotsl1, remove_x = remove_x, 
                                remove_y = remove_y)
  plotsl2 <- map2(plotsl1, bar_l, tradeoff_add_inset) 
  
  comb <- combine_panels_labs(plotsl2,
                                xlab = 'Mean Sagebrush Ecological Integrity',
                                ylab = 'Expected burned area (%/year)')
  
  prefix <- if(gcm_path) '' else 'noGCM_'

  ggsave(paste0("figures/sei/tradeoff/", "sei", "-scd-adj-vs-ba_perc_dotplot_", 
              'c3eco_', prefix, rcp, yr_lab, "_", suffix, ".png"), 
         plot = comb, dpi = 600, width = width_tradeoff, height = height_tradeoff)
})


# * conceptual figure -----------------------------------------------------

dummy <- tibble(
  group = rep(c('a', 'b', 'c'), each = 4),
  graze = rep(graze2factor(c("Very Heavy", "Heavy", "Moderate", "Light")), 
              times = 3),
  y = c(1, 2, 3, 4,
        1, 2, 3, 4,
        1, 1, 1, 1),
  x = c(rep(0.75, 4),
        c(1, 2, 3, 4)*1.2 + 0.5,
        (2:5)*1.2)
)

g <- ggplot(dummy, aes(x, y)) +
  geom_path(aes(group = group), alpha = 0.5) +
  geom_point(aes(color = graze)) +
  scale_color_graze() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Mean Sagebrush Ecological Integrity',
       y = 'Burned Area') +
  expand_limits(x = c(0.2, 6.5),
                y = c(0.5, 4.5))

ggsave("figures/sei/tradeoff/conceptual_fig_v1.png", 
       plot = g, dpi = 600,
       width = 4, height = 3)

# attribution -------------------------------------------------------------


# * burned area--attribution ------------------------------------------------

drivers2 <- drivers1 %>% 
  select(-id)

ba_gcm2 <- ba_gcm1 %>% 
  select(-type)

ba_gcm3 <- drivers2 %>% 
  filter(type == "climate") %>% 
  select(-graze, -type) %>%  # grazing doesn't apply to climate variables
  pivot_wider(values_from = 'mean',
              names_from = variable) %>% 
  right_join(ba_gcm2, by = join_by(ecoregion, RCP, years, GCM))
  
ba_gcm4 <- drivers2 %>% 
  filter(type == "biomass") %>% 
  select(-type) %>% 
  pivot_wider(values_from = 'mean',
              names_from = variable) %>% 
  right_join(ba_gcm3, by = join_by(ecoregion, RCP, years, GCM, graze)) %>% 
  rename(PSP = psp) %>% 
  df_factor() %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE))

driver_vars <- c("MAP", "MAT", "PSP", "Aherb", "Pherb", "Sagebrush")
ba_gcm5 <- ba_gcm4 %>% 
  pivot_longer(cols = all_of(driver_vars),
               values_to = 'mean_driver',
               names_to = 'driver') %>% 
  mutate(GCM = factor(GCM, levels = c('Current', names(cols_GCM1))),
         driver = factor(driver, levels = driver_vars)) %>% 
  df_factor()

ba_current <- ba_gcm5 %>% 
  filter(GCM == 'Current')

rcp_year <- unique(ba_gcm4$rcp_year) %>% 
  str_subset('RCP')


plots <- map(rcp_year, function(x) {
  ba_gcm5 %>% 
    filter(rcp_year == x) %>% 
    ggplot(aes(mean_driver, area)) +
    geom_smooth(aes(linetype = graze), se = FALSE, color = 'gray') +
    geom_point(aes(shape = graze, color = GCM)) +
    geom_point(data = ba_current,
               aes(shape = graze, color = 'Historical')) +
    facet_grid(ecoregion~driver, scales = 'free', switch = 'x') +
    scale_color_manual(values = cols_GCM2,
                       name = 'GCM (or historical)') +
    scale_linetype(name = 'Grazing') +
    scale_shape(name = 'Grazing')+
    labs(x = " ",
         y = lab_ba0,
         subtitle = x,
         caption = 'Mean of fire probability predictor variable was calculated across
         pixels in ecoregion, for a given GCM') +
    theme(strip.placement.x = 'outside')
  
})

pdf(paste0("figures/fire/area/expected_ba_vs_driver_by-GCM_", suffix, '.pdf'),
    width = 13, height = 10*mr)
print(plots)
dev.off()


# * delta fire prob attribution -------------------------------------------


for (rcp in rcps) {
df_gcm <- one_change_gcm2 %>% 
  filter(RCP == rcp) 

df_med <- one_change_smry2 %>% 
  filter(RCP == rcp) 

cap1 <- paste0(unique(df_gcm$rcp_years, 
                      '\nSite level data weighted by area interpolated to'))

g <- ggplot(mapping = aes(delta_pred_value, delta_1var, color = graze)) +
  facet_grid(region~pred_var, scales = 'free', switch = 'x') +
  theme(strip.placement.x = 'outside') +
  #scale_color_manual(values = scales::alpha(cols_graze, alpha = 0.25)) +
  scale_color_graze() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2, linewidth = 0.5) +
  labs(x = 'Change in predictor variable (future - historical)',
       y = 'Change in fire frequency due to change in single predictor (#/100 yrs)',
       caption = cap1)

g_gcm0 <- g+ coord_cartesian(ylim = c(-2, 2)) +
  labs(subtitle = 'Results shown by GCM')
  
g_gcm <- g_gcm0 + 
  geom_smooth(
    data = df_gcm,
    aes(weight = weight, group = interaction(GCM, graze)),
    method.args = list(span = 1, degree = 1), # more smoothing
    method = 'loess',
    se = FALSE,
    linewidth = 0.2
  ) 
ggsave_delta_prob(g_gcm, rcp, 'by-GCM')

g_gcm_lm <- g_gcm0 + 
  geom_smooth(
    data = df_gcm,
    aes(weight = weight, group = interaction(GCM, graze)),
    method = 'lm',
    se = FALSE,
    linewidth = 0.2
  ) 
ggsave_delta_prob(g_gcm_lm, rcp, 'by-GCM-lm')
g_gcm_point <-  g_gcm0 + 
  geom_point(
    data = df_gcm,
    aes(group = interaction(GCM, graze)),
    alpha = 0.5,
    size = 0.3
  ) 

ggsave_delta_prob(g_gcm_point, rcp, 'by-GCM-point')

g_med0 <- g +
  labs(subtitle = 'Median across GCMs (by site) shown')

g_med1 <- g_med0 +
  geom_smooth(data = df_med,
    aes(weight = weight),
    method.args = list(span = 1, degree = 1),
    method = 'loess',
    se = FALSE,
    linewidth = 0.75
  )

ggsave_delta_prob(g_med1, rcp, 'med')

ggsave_delta_prob(plot = g_med0 + geom_point(data = df_med, size = 0.5),
                  rcp, 
                  'med-point')
}

# ** summaries across pixels/sites ----------------------------------------

for (rcp in rcps) {
df_smry <- one_change_smry_sas1 %>% 
  filter(RCP == rcp)
df_gcm <- one_change_gcm_sas1 %>% 
  filter(RCP == rcp)

g <- ggplot(mapping = aes(delta_pred_value, delta_1var)) +
  facet_grid(region~pred_var, scales = 'free', switch = 'x') +
  theme(strip.placement.x = 'outside') +
  geom_point(data = df_gcm, aes(color = graze), size = 0.5) +
  geom_path(data = df_gcm, aes(group = GCM), alpha = 0.5,
            linewidth = 0.5) +
  geom_point(data = df_smry , aes(color = graze)) +
  geom_path(data = df_smry ) +
  scale_color_graze() +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.2, linewidth = 0.5) +
  labs(x = 'Mean change in predictor variable (future - historical)',
       y = 'Mean change in fire frequency due to change in single predictor (#/100 yrs)',
       caption = paste(rcp, 'summarized across sites'))

ggsave_delta_prob(g, rcp, smry = 'sas')

}
# combined figures --------------------------------------------------------
# combined figure for manuscript

# expected  burned darea
args <- expand_grid(rcp = rcps,
                    region = ecoregions)
pmap(args, function(rcp, region) {
  ba3_tmp <- ba3 %>% 
    filter(ecoregion == region, RCP %in% c("Current", rcp)) %>% 
    mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                                add_newline = TRUE))
  
  total_area <- unique(ba3_tmp$area_total)
  
  
  
  g_ba1 <- ggplot(ba3_tmp, aes(x = graze, y = area_median_perc)) +
    geom_errorbar(aes(ymin = area_low_perc, ymax = area_high_perc), 
                  position = position_dodge(width = 0.5), width = 0) +
    geom_point(aes(y = area_median_perc), 
               position = position_dodge(width = 0.5)) +
    # scale_y_continuous(sec.axis = sec_axis(transform = \(x) total_area*x/100,
    #                                        name = lab_ba0)) +
    facet_wrap(~rcp_year) +
    adjust_graze_x_axis(xlab = NULL) +
    labs(y = lab_ba1)
  
  
  # years since fire
  
  g_bar1 <- area_age_group4_pw %>% 
    filter(RCP %in% c('Current', rcp), ecoregion == region) %>% 
    ggplot(aes(x = graze, y = area_median_perc, fill = age_group)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~rcp_year) +
    adjust_graze_x_axis(xlab = NULL) +
    scale_fill_manual(values = rev(cols_agegroup),
                      name = lab_agegroup) +
    labs(y = lab_areaperc0)
  
  g_dom <- dom_drivers2 %>% 
    filter(RCP == rcp, region == !!region) %>% 
    ggplot(aes(x = graze, y = perc_area, fill = dom_driver)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~rcp_year) +
    adjust_graze_x_axis() +
    scale_fill_manual(values = c(cols_pred_vars, 'None' = 'grey'),
                      name = "Primary driver of\nchange in fire\nprobability") +
    labs(y = lab_areaperc0)
  
  comb <- guide_area() +  g_dom + plot_layout(guides = 'collect', 
                                              widths = c(0.8, 1))
  plots <- list(g_ba1, g_bar1, comb)
  
  g <- wrap_plots(plots, ncol = 1)&theme(plot.margin = unit(c(0, 0, 0, 3), "pt"))
  
  if(region != entire) {
    g <- g + plot_annotation(subtitle = region)
  }
  
  ggsave(paste0("figures/fire/area/comb_ba-age-driver_", 
                words2abbrev(region), vr_name, "_", rcp, yr_lab, "_", v, "_", 
                runv, '.png'),
         plot = g,
         height = 9, width = 5, dpi = 600)
})


# explanatory figures -----------------------------------------------------

if (explanatory_figures) {
  p <- seq(from = 0.0001, to = 0.04, length.out = 100)
  
  df1 <- expand_grid(
    p = p,
    age_group = factor(names(age_groups), levels = names(age_groups))
  )
  
  df1$a <- map_dbl(age_groups[df1$age_group], \(x) x['a']) # lower bound of group
  df1$b <- map_dbl(age_groups[df1$age_group], \(x) x['b']) # upper bound
  df1$prob_in_group <- pmap_dbl(df1[c('p', 'a', 'b')], prob_geometric_rule)
  
  # probabilities should sum to 1 
  test <- df1 %>% 
    group_by(p) %>% 
    summarize(tot = sum(prob_in_group)) %>% 
    pull(tot) 
  
  stopifnot(round(test, digits = 5) == 1)
  
  g <- ggplot(df1, aes(p, prob_in_group)) +
    geom_line() +
    facet_wrap(~age_group) +
    theme_bw() +
    labs(x = 'Fire probability (per year)',
         y = 'Probability of a site being in the given years since fire group',
         caption = paste("This is to help understand the relationship between", 
                         "probability of occurrence (per year)",
                         "\nand the probability of no occurrence (no fire) for",
                         "a given interval ('year since fire group'),\n calculated via geometric rule"))
  
  jpeg("figures/fire/area/explanatory_fig_geometric-rule.jpg",
       units = 'in', width = 7, height = 6, res = 600)
  print(g)
  dev.off()
}

} # end create figs