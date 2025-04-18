# Purpose: TBD

# Author: Martin Holdrege

# Started: April 11, 2025

# params ------------------------------------------------------------------

v_interp <- 'v4' # interpolation version
run <- "fire1_eind1_c4grass1_co20_2503"
q_curve_fig <- TRUE # create figure showing the 'biomass' q curves



# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
source("src/fig_params.R")
source("src/fig_functions.R")
source('src/general_functions.R')
source('src/SEI_functions.R')
source('src/mapping_functions.R')
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

path_rds <- 'data_processed/site_means/summarize_bio.RDS'

if(!file.exists(path_rds)) {
  source("scripts/02_summarize_bio.R") # where needed dataframes are created
}
l <- readRDS(path_rds) # created in 02_sommarize_bio.R

# weight (based on area of interpolation) corresponding to each site
# in each ecoregion
# created in 04_interpolation_weights.R
w1 <- read_csv(paste0('data_processed/interpolation_data/interpolation_weights_', 
                       v_interp, '.csv'),
               show_col_types = FALSE)

bio_gcm1 <- l$pft5_bio1 # 
bio2 <- l$pft5_bio2# summarized across GCMs
fire_med1 <- l$fire_med1

# combine site level and weights ------------------------------------------

bio3 <- bio2 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  filter_clim_extremes() %>% 
  filter(run == !!run) %>% 
  df_factor() 


# longer format
bio4 <- bio3 %>% 
  select(-matches('utilization'), -matches('indivs')) %>% 
  pivot_longer(matches('biomass'),
               names_to = 'summary',
               values_to = 'biomass') %>% 
  mutate(summary = str_extract(summary, '(?<=_).*'),
         summary = ifelse(is.na(summary), "median", summary),
         summary = summary2factor(summary),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE))

bio_gcm2 <- bio_gcm1 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  filter_clim_extremes() %>% 
  filter(run == !!run) %>% 
  df_factor() %>% 
  filter(region != levels(region)[1])


# calculate Q's and SEI ---------------------------------------------------

pfts <- c('Sagebrush', 'Pherb', 'Aherb')

# need to calculate for GCM wise q values and SEI, and then summaryies
# b/ SEI is multiplicative median sei need to be from median Q values


q_gcm1 <- bio_gcm2 %>% 
  filter(region != levels(region)[1],
         PFT %in% pfts) %>%  # entire study area
  group_by(PFT, region) %>% 
  mutate(Q = bio2q(biomass, as.character(unique(PFT)), as.character(unique(region))),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE))

sei1 <- q_gcm1 %>% 
  select(-matches('biomass'), - matches('indivs'), -matches('utilization')) %>% 
  pivot_wider(values_from = 'Q',
              names_from = 'PFT',
              names_prefix = 'Q_') %>% 
  mutate(SEI = Q_Aherb*Q_Sagebrush*Q_Pherb) %>% 
  select(-matches('^Q_')) %>% 
  group_by(run, years, RCP, graze, site, weight, region, rcp_year) %>% 
  # summarize across GCMs
  summarize_across_GCMs(var = 'SEI')

sei2 <- sei1 %>% 
  # 'total' region
  mutate(region = levels(region)[1]) %>% 
  # individual regions
  bind_rows(sei1) 

q1 <- q_gcm1 %>% 
  group_by(run, years, RCP, rcp_year, graze, site, weight, region, PFT) %>% 
  summarize_across_GCMs(var = 'Q')

q2 <- q1 %>% 
  # the entire study area is all other regions combined
  # but note that q's are calculated seperately for each ecoregion
  # so a given site may appear multiple times for the 'total'
  mutate(region = levels(region)[1]) %>% 
  bind_rows(q1)

# boxplots ----------------------------------------------------------------

# *biomass ----------------------------------------------------------------

plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(bio4, PFT == pft),
                y_string = 'biomass',
                ylab = lab_bio0,
                subtitle = pft) +
    expand_limits(y = 0)
})

g2 <- combine_grid_panels1(plots)

png(paste0("figures/biomass/bio_weighted_by-region_3pft_boxplot_", run, ".png"),
    width = 8, height = 8, units = 'in', res = 600)
g2
dev.off()

# * Q and SEI -------------------------------------------------------------

plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(q2, PFT == pft),
                y_string = 'Q',
                ylab = lab_q0,
                subtitle = pft) +
    expand_limits(y = c(0, 1))
})
n <- length(plots)
plots[[n]] <- plots[[n]] + theme(strip.text.y = element_blank())

g_sei <- weighted_box1(df = sei2,
                     y_string = 'SEI',
                     ylab = lab_sei0,
                     subtitle = 'SEI')


g2 <- combine_grid_panels1(plots, remove_y_axis_labels = TRUE)
g3 <- g2 + g_sei + plot_layout(guides = 'collect')

png(paste0("figures/sei/q-sei_weighted_by-region_3pft_boxplot_", run, ".png"),
    width = 11, height = 8, units = 'in', res = 600)
g3&theme(legend.position = 'bottom')
dev.off()


# q-curve -----------------------------------------------------------------

if(q_curve_fig) {
  tmp <- bio_gcm1 %>% 
    filter(PFT %in% pfts) %>% 
    group_by(PFT)
  n <- 1000
  
  # Q's for the complete range of biomass in the dataset
  
  q_range <- tmp %>% 
    group_split() %>% 
    map(\(x) range(x$biomass)) %>% 
    set_names(., as.character(group_keys(tmp)[[1]])) %>% 
    map(\(x) seq(x[1], x[2], length.out = n)) %>% 
    bind_cols() %>% 
    mutate(row = 1:n()) %>% 
    pivot_longer(-row, names_to = 'PFT', values_to = 'biomass') %>% 
    right_join(expand_grid(region = unique(q1$region),row = 1:n), 
               by = 'row', 
               relationship = 'many-to-many') %>% 
    group_by(PFT, region) %>% 
    mutate(Q = bio2q(biomass, pft = unique(as.character(PFT)), 
                     region = unique(as.character(region)))) %>% 
    select(-row) %>% 
    df_factor()
  
  g <- ggplot(q_range, aes(biomass, Q, color = region, linetype = region)) +
    geom_line() + 
    facet_wrap(~PFT, scales = 'free_x') +
    scale_color_manual(values = cols_ecoregion) +
    labs(y = lab_q0,
         x = 'Biomass (from stepwat, prior to quantile mapping)',
         caption = 'Q calculated by quantile mapping Pherb and Aherb to 
         RAP scale, then to cover via linear eq., for Sagebrush, straight
         quantile mapping to RCMAP cover') +
    theme(legend.position = 'bottom')
  
  ggsave('figures/sei/q-curves-biomass_v1.png', plot = g, dpi = 600,
         width = 6, height = 4)
}


# comparing to fire -------------------------------------------------------

# 'heat map' type figure showing binned fire on the y axis, grazing
# on the x and SEI (or area) as the color

# * prepare dataframes ----------------------------------------------------


fire_breaks <- c(seq(0, 2.5, by = 0.25), 3, 4)

stopifnot(max(fire_breaks) > max(fire_med1$fire_prob))

sei3 <- sei2 %>% 
  filter(summary == 'median') %>% 
  left_join(fire_med1, by = join_by(run, years, RCP, graze, site))

df_fire_bins1 <- sei3 %>% 
  mutate(fire_bin = cut(fire_prob, breaks = fire_breaks, right = FALSE)) %>% 
  group_by(run, years, RCP, graze, region, rcp_year, fire_bin) %>% 
  summarize(n = n(),
            SEI = weighted.mean(SEI, w = weight),
            weight = sum(weight),
            .groups = 'drop_last') %>% 
  mutate(percent_area = weight/sum(weight) *100)
  

# * create figures ----------------------------------------------------------

base_fire_bins <- function() {
  list(geom_tile(),
       labs(x = lab_graze,
            y =  lab_firep0),
    theme(axis.text.y = element_text(size = rel(0.8)),
          axis.text.x = element_text(angle = 45, hjust = 1)),
    facet_grid(region~rcp_year))
}


g <- ggplot(df_fire_bins1, aes(x = graze, y = fire_bin, fill = SEI)) +
  base_fire_bins() +
  scale_fill_gradientn(
    colors = c("#eee1ba", "#a6611a", "#78c679", "#006837", "#2166ac", "#053061"),
    values = scales::rescale(c(0, 0.173, 0.174, 0.431, 0.432, 1)),
    guide = "colorbar",
    name = "SEI"
  ) 

ggsave(paste0("figures/sei/sei_by-fire-graze_", run, ".png"), 
       plot = g, dpi = 600,
       width = 8, height = 9)

g <- ggplot(df_fire_bins1, aes(x = graze, y = fire_bin, fill = percent_area)) +
  base_fire_bins() +
  scale_fill_viridis_c(option = "C", name = "% of area", trans = 'sqrt')
  

ggsave(paste0("figures/sei/area-perc_by-fire-graze_", run, ".png"), 
       plot = g, dpi = 600,
       width = 8, height = 9)
