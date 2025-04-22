# Purpose: Create figures of expected burned area and area in different
# fire return inverval categories

# Author: Martin Holdrege

# Started: March 25, 2025

# params

# dependencies

# params ------------------------------------------------------------------

run <- 'fire1_eind1_c4grass1_co20_2503'

v <- 'v1' # version of input files (and for now also used in output file names)

# create additional figures to explore how age groups/fire probability relate
# (not using actual data)
explanatory_figures <- TRUE 

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(patchwork)
source("src/general_functions.R")
source("src/fig_functions.R")
source("src/fig_params.R")
source("src/probability_functions.R")
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

# created in the 05_fire_area.R script
# expected burned area
ba3a <- read_csv(
  paste0("data_processed/area/expected-burn-area_", v, "_", run,".csv"),
  show_col_types = FALSE)

# expected burned area by gcm
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                            v, "_", run, ".csv"),
                    show_col_types = FALSE)

# expected area in different years since fire age groups
area_age_group3 <- read_csv(
  paste0("data_processed/area/area-by-age-group_", v, "_", run, ".csv"),
  show_col_types = FALSE)

# area of our study region in each of 3 ecoregions
area_eco <- read_csv(paste0("data_processed/area/ecoregion-area_", v,".csv"),
                     show_col_types = FALSE)

# created in the 04_summarize_fire_drivers.R script
# means, by ecoregion, of the drivers of fire probability
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, 
                            '_fire-driver-means_by-ecoregion.csv'))

# vectors -----------------------------------------------------------------

ecoregions <- region_factor(area_eco$ecoregion) %>% 
  levels()
age_groups <- create_age_groups()

# fig params --------------------------------------------------------------

line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot
line_loc2 <- 1:4 + 0.5
height_3p <- 6 # width for 3 panel figures
width_3p <- 7# height for 3 panel figures

# prepare dataframes ------------------------------------------------------

ba3 <- df_factor(ba3a) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE))

area_age_group3 <- df_factor(area_age_group3) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         # not reproducible (yet)
         age_group = factor(age_group, levels = names(age_groups)),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE))

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
                                           name = '% of ecoregion')) +
    scale_color_manual(values = cols_graze, name = 'Grazing') +
    geom_vline(xintercept = line_loc2, linetype = 2) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 25, hjust = 1, size = rel(0.8))) +
    labs(x = NULL,
         y = lab_ba0,
         subtitle = region)
})

g <- patchwork::wrap_plots(plots) +
  plot_layout(guides = 'collect', ncol = 2)

g2 <- g&theme(legend.position = 'bottom')

jpeg(paste0("figures/fire/area/expected_ba_dotplot_", v, "_", run, '.jpg'),
     units = 'in', width = width_3p, height = height_3p, res = 600)
g2
dev.off()


# age_group figs -----------------------------------------------------


# *grouped by grazing level -----------------------------------------------

plots <- map(ecoregions, function(region) {
  total_area <- area_eco$area[area_eco$ecoregion == region]
  
  df <- area_age_group3  %>% 
    filter(ecoregion == region) 
  anno_data <- box_anno(
    df = df,
    var = 'area_high',
    id = "id2",
    mult = 0.1,
    group_by = c('graze', 'age_group'),
    anno_same_across_panels = TRUE
  )
  ggplot(df, aes(x = id2, y = area_median)) +
    geom_errorbar(aes(ymin = area_low, ymax = area_high), width = 0) +
    geom_point(aes(color = RCP)) + # not showing outliers as points
    scale_x_discrete(labels = id2year) +
    scale_y_continuous(sec.axis = sec_axis(transform = \(x) x/total_area*100,
                                           name = '% of ecoregion')) +
    scale_color_manual(values = cols_rcp, name = "Scenario") +
    geom_text(data = anno_data,
              aes(x, y, label = graze, fill = NULL),
              size = 2.5) +
    geom_vline(xintercept = line_loc, linetype = 2) +
    facet_wrap(~age_group, nrow = 1) +
    labs(y = 'Expected area in age class (ha/yr)',
         subtitle = region,
         x = lab_yrs)
})

g <- patchwork::wrap_plots(plots) +
  plot_layout(guides = 'collect', ncol = 1, axis_titles = 'collect')

g2 <- g&theme(legend.position = 'bottom')


jpeg(paste0("figures/fire/area/area_age_group_", v, "_", run, '.jpg'),
     units = 'in', width = 14, height = 12, res = 600)
g2
dev.off()


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
                                     size = rel(0.8))) +
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
  
  g1 <- df %>% 
    filter(age_group != old_group) %>% 
    fig_area_by_scenario_group(total_area = total_area) + 
    labs(subtitle = region)
  
  
  g2 <- df %>% 
    filter(age_group == old_group) %>% 
    fig_area_by_scenario_group(total_area = total_area)
  
  list(g1, g2)
})

plots <- reduce(plots0, c)
g <- patchwork::wrap_plots(plots, ncol = 2, guides = 'collect',
                           tag_level = 'keep',byrow = TRUE, widths = c(3, 1)) +
  plot_layout(axis_titles = 'collect')

g2 <- g&theme(legend.position = 'bottom',
              axis.title.y = element_text(),
              axis.title.y.right = element_text()
)

jpeg(paste0("figures/fire/area/area_age_group_scen-grouping_", v, "_", run, '.jpg'),
     units = 'in', width = 10, height = 7, res = 600)
g2
dev.off()


# burned area--attribution ------------------------------------------------

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

driver_vars <- c("MAP", "MAT", "PSP", "Aherb", "Pherb")
ba_gcm5 <- ba_gcm4 %>% 
  pivot_longer(cols = driver_vars,
               values_to = 'mean_driver',
               names_to = 'driver') %>% 
  mutate(GCM = factor(GCM, levels = c('Current', names(cols_GCM1))),
         driver = factor(driver, levels = driver_vars))

ba_current <- ba_gcm5 %>% 
  filter(GCM == 'Current')

rcp_year <- unique(ba_gcm4$rcp_year) %>% 
  str_subset('RCP') %>% 
  rev()


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

pdf(paste0("figures/fire/area/expected_ba_vs_driver_by-GCM_", v, "_", run, '.pdf'),
    width = 12, height = 7)
plots
dev.off()


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
