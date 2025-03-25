# Purpose: Create figures of expected burned area and area in different
# fire return inverval categories

# Author: Martin Holdrege

# Started: March 25, 2025

# params

# dependencies

# params ------------------------------------------------------------------

run <- 'fire1_eind1_c4grass1_co20_2502'
t
v <- 'v1' # version of input files (and for now also used in output file names)

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
theme_set(theme_custom1())
source("src/general_functions.R")
source("src/fig_functions.R")
source("src/fig_params.R")

# read in data ------------------------------------------------------------

ba2 <- read_csv(
  paste0("data_processed/area/expected-burn-area_", v, "_", run,".csv"),
  show_col_types = FALSE)

fri_area2 <- read_csv(
  paste0("data_processed/area/area-by-fri-class_", v, "_", run, ".csv"),
  show_col_types = FALSE)


# fig params --------------------------------------------------------------

line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot
height_3p <- 4 # width for 3 panel figures
width_3p <- 6# height for 3 panel figures
# prepare dataframes ------------------------------------------------------

ba3 <- df_factor(ba2) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)))

fri_area3 <- df_factor(fri_area2) %>% 
  arrange(graze, RCP, years)  %>% 
  mutate(id2 = paste(RCP, years, graze, sep = '_'),
         id2 = factor(id2, levels = unique(id2)),
         # not reproducible (yet)
         fri_label = factor(fri_label, levels = fri_labels))



# expected burned area figs --------------------------------------------------------

g <- ggplot(ba3, aes(x = id2, y = area)) +
  geom_point(aes(color = RCP)) + # not showing outliers as points
  facet_wrap(~ecoregion, scales = 'free_y', ncol = 2) +
  scale_x_discrete(labels = id2year) +
  scale_color_manual(values = cols_rcp, name = "Scenario") +
  geom_text(data = ~box_anno(.,var = 'area',
                             id = "id2",
                             mult = 0.1,
                             group_by = c("ecoregion", 'graze')),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = legend_pos_box1,
      axis.text = element_text(size = 7)) +
  labs(x = lab_yrs,
       y = lab_ba0)

jpeg(paste0("figures/fire/area/expected_ba_dotplot_", v, "_", run, '.jpg'),
     units = 'in', width = width_3p, height = height_3p, res = 600)
g
dev.off()



# # FRI interval figs -----------------------------------------------------

g <- ggplot(fri_area3, aes(x = id2, y = area)) +
  geom_point(aes(color = RCP)) + # not showing outliers as points
  facet_wrap(fri_label~ecoregion, scales = 'free_y', ncol = 3) +
  scale_x_discrete(labels = id2year) +
  scale_color_manual(values = cols_rcp, name = "Scenario") +
  geom_text(data = ~box_anno(.,var = 'area',
                             id = "id2",
                             mult = 0.1,
                             group_by = c("ecoregion", 'graze', 'fri_label')),
            aes(x, y, label = graze, fill = NULL),
            size = 2.5) +
  geom_vline(xintercept = line_loc, linetype = 2) +
  theme(legend.position = legend_pos_box1,
        axis.text = element_text(size = 7)) +
  labs(x = lab_yrs,
       subtitle = 'Area falling into different FRI categories',
       y = 'Area (ha)')

jpeg(paste0("figures/fire/area/area_by-fri_dotplot_", v, "_", run, '.jpg'),
     units = 'in', width = 8, height = 13, res = 600)
g
dev.off()
