# Purpose: Create figures of expected burned area and area in different
# fire return inverval categories

# Author: Martin Holdrege

# Started: March 25, 2025

# params

# dependencies

# params ------------------------------------------------------------------

run <- 'fire1_eind1_c4grass1_co20_2503'

v <- 'v1' # version of input files (and for now also used in output file names)

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(patchwork)
theme_set(theme_custom1())
source("src/general_functions.R")
source("src/fig_functions.R")
source("src/fig_params.R")
source("src/probability_functions.R")

# read in data ------------------------------------------------------------

ba3a <- read_csv(
  paste0("data_processed/area/expected-burn-area_", v, "_", run,".csv"),
  show_col_types = FALSE)

area_age_group3 <- read_csv(
  paste0("data_processed/area/area-by-age-group_", v, "_", run, ".csv"),
  show_col_types = FALSE)

area_eco <- read_csv(paste0("data_processed/area/ecoregion-area_", v,".csv"),
                     show_col_types = FALSE)

# vectors -----------------------------------------------------------------

ecoregions <- area_eco$ecoregion

# fig params --------------------------------------------------------------

line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot
height_3p <- 6 # width for 3 panel figures
width_3p <- 7# height for 3 panel figures
# prepare dataframes ------------------------------------------------------

ba3 <- df_factor(ba3a) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)))

area_age_group3 <- df_factor(area_age_group3) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         # not reproducible (yet)
         age_group = factor(age_group, levels = names(create_age_groups())))

# expected burned area figs --------------------------------------------------

plots <- map(ecoregions, function(region) {
  total_area <- area_eco$area[area_eco$ecoregion == region]
  ba3 %>% 
    filter(ecoregion == region) %>% 
    ggplot(aes(x = id2, y = area_median)) +
    geom_errorbar(aes(ymin = area_low, ymax = area_high), width = 0) +
    geom_point(aes(color = RCP)) + # not showing outliers as points
    scale_x_discrete(labels = id2year) +
    scale_y_continuous(sec.axis = sec_axis(transform = \(x) x/total_area*100,
                                           name = '% of ecoregion')) +
    scale_color_manual(values = cols_rcp, name = "Scenario") +
    geom_text(data = ~box_anno(.,var = 'area_high',
                               id = "id2",
                               mult = 0.1,
                               group_by = c('graze')),
              aes(x, y, label = graze, fill = NULL),
              size = 2.5) +
    geom_vline(xintercept = line_loc, linetype = 2) +
    theme(legend.position = legend_pos_box1,
          axis.text = element_text(size = 7)) +
    labs(x = lab_yrs,
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

plots <- map(ecoregions, function(region) {
  total_area <- area_eco$area[area_eco$ecoregion == region]
  
  df <- area_age_group3  %>% 
    filter(ecoregion == region) 
  anno_data <- box_anno(
    df = df,
    var = 'area_low',
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
