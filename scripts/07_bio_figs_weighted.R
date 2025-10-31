# Purpose: Biomass etc figures, weighted by interpolation area

# Author: Martin Holdrege

# Started: April 11, 2025

# dependencies ------------------------------------------------------------

library(tidyverse)
library(patchwork)
source("src/fig_params.R")
source("src/fig_functions.R")
source('src/general_functions.R')
source('src/SEI_functions.R')
source('src/mapping_functions.R')
theme_set(theme_custom1())

# params ------------------------------------------------------------------

source('src/params.R')
yr_lab <- opt$yr_lab
years <- opt$years
run <- opt$run
runv <- paste0(run, v_interp)
q_curve_fig <- FALSE # create figure showing the 'biomass' q curves
v <- '_v2'
suffix <- paste0(run, v, vr_name, yr_lab) # for figures
pfts <- c("Sagebrush", 'Pherb', 'Aherb') # code won't work with any other pfts

regions <- region_factor(v = vr, return_levels = TRUE)
entire <- regions[1]
ref_graze <- opt$ref_graze
# read in data ------------------------------------------------------------

path_rds <- 'data_processed/site_means/summarize_bio.RDS'

if(!file.exists(path_rds)) {
  source("scripts/02_summarize_bio.R") # where needed dataframes are created
}
l <- readRDS(path_rds) # created in 02_sommarize_bio.R
l <- filter_list_run(l, run)

# weight (based on area of interpolation) corresponding to each site
# in each ecoregion
# created in 04_interpolation_weights.R
w1 <- read_csv(paste0('data_processed/interpolation_data/interpolation_weights_', 
                       v_interp, vr_name, '.csv'),
               show_col_types = FALSE)

bio_gcm1 <- l$pft5_bio1 # 
bio2 <- l$pft5_bio2# summarized across GCMs
fire_med1 <- l$fire_med1
pft5_bio_d2 <- l$pft5_bio_d2
pft5_d_wgcm <- l$pft5_d_wgcm

# weighted means and percentiles of change biomass and climate variables
# created in "scripts/06_summarize_fire_drivers.R"
driver1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                  '_delta-fire-driver_by-ecoregion.csv'))

# combine site level and weights ------------------------------------------
w1 <- df_factor(w1)

bio3 <- bio2 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  filter_clim_extremes(years = years) %>% 
  df_factor() 


# longer format
bio4 <- bio3 %>% 
  pivot_longer(
    cols = matches("^(utilization|biomass|indivs)"),
    names_to = c(".value", "summary"),
    names_pattern = "^(utilization|biomass|indivs)_?(.*)$"
  ) %>% 
  mutate(summary = ifelse(summary == "", "median", summary),
         summary = summary2factor(summary),
         rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                              add_newline = TRUE))

bio_gcm2 <- bio_gcm1 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  filter_clim_extremes(years = years) %>% 
  df_factor() %>% 
  filter(region != levels(region)[1])

# % change (climate effect)

prepare_bio_diff <- function(df) {
  df %>% 
    filter_clim_extremes(years = years) %>% 
    filter(run == !!run) %>% 
    rename(bio_diff_median = bio_diff) %>% 
    select(run, site, years, RCP, PFT, graze, matches('bio_diff')) %>% 
    pivot_longer(cols = matches('bio_diff_'),
                 values_to = 'bio_diff',
                 names_to = 'summary',
                 names_pattern = '([[:alpha:]]+$)') %>% 
    left_join(w1[w1$region == entire, ], 
              by = 'site') %>% 
    df_factor(v = vr) %>% 
    mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                                add_newline = TRUE))
}

pft5_bio_d3 <- pft5_bio_d2 %>% 
  prepare_bio_diff()


# % change (grazing effect)
pft5_d_wgcm2 <- pft5_d_wgcm %>% 
  prepare_bio_diff()
  
# other summaries ---------------------------------------------------------

# % of total utilization
bio_pcent1 <- bio4 %>% 
  filter(summary == 'median',
         PFT %in% pfts) %>% 
  group_by(region, run, years, RCP, rcp_year, graze, PFT) %>% 
  summarize(util_median = spatstat.univar::weighted.median(utilization, w = weight),
            .groups = 'drop_last') %>% 
  mutate(util_pcent = util_median/sum(util_median)*100)

# setup so can plot climate variables against biomass variables
driver_bio1 <- driver1 %>% 
  filter(variable %in% pfts)

driver2 <- driver1 %>% 
  filter(!variable %in% pfts) %>% 
  select(-graze, -weight) %>% 
  right_join(driver_bio1, by = c('years', 'RCP', 'GCM', 'region'),
             suffix = c('_clim', '_bio'),
             relationship = "many-to-many") %>% 
  df_factor() %>% 
  mutate(variable_clim = driver2factor(variable_clim),
         variable_bio = driver2factor(variable_bio, include_sagebrush = TRUE)) %>% 
  droplevels() %>% 
  arrange(region, RCP, years)

# calculate Q's and SEI ---------------------------------------------------



# need to calculate for GCM wise q values and SEI, and then summaryies
# b/ SEI is multiplicative median sei need to be from median Q values


q_gcm1 <- bio_gcm2 %>% 
  filter(region != levels(region)[1],
         PFT %in% pfts) %>%  # entire study area
  group_by(PFT, region) %>% 
  mutate(Q = bio2q(biomass, as.character(unique(PFT)), as.character(unique(region)),
                   v = vr),
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
  bind_rows(sei1) %>% 
  df_factor()

# for use in the 08_fire_area_figs.R script
saveRDS(sei2, paste0('data_processed/temp_rds/sei_df', vr_name, yr_lab, '.rds'))

q1 <- q_gcm1 %>% 
  group_by(run, years, RCP, rcp_year, graze, site, weight, region, PFT) %>% 
  summarize_across_GCMs(var = 'Q')

q2 <- q1 %>% 
  # the entire study area is all other regions combined
  # but note that q's are calculated seperately for each ecoregion
  # so a given site may appear multiple times for the 'total'
  mutate(region = levels(region)[1]) %>% 
  bind_rows(q1) %>% 
  df_factor()

# boxplots ----------------------------------------------------------------

# *biomass ----------------------------------------------------------------

plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(bio4, PFT == pft),
                y_string = 'biomass',
                ylab = lab_bio0,
                subtitle = pft) +
    expand_limits(y = 0)
})

g2 <- combine_grid_panels1(plots, remove_y_axis_labels = FALSE)

png(paste0("figures/biomass/bio_weighted_by-region_3pft_boxplot_", suffix, ".png"),
    width = 8, height = 11, units = 'in', res = 600)
g2
dev.off()

# absolute biomass for current and %change for future
pft5 <- c(pfts, 'C3Pgrass', 'C4Pgrass')
args <- list(
  pfts = list(pfts, pft5, pfts, pft5),
  RCP = c('RCP45', 'RCP45', 'RCP85', 'RCP85')
)

ylab_diff_gref = paste0('\u0394Biomass (%), relative to ', str_to_lower(ref_graze),
                        ' grazing')

ylab_diff_cref = '\u0394Biomass (%), relative to historical climate'

# figure for entire region
pmap(args, function(pfts, RCP) {
  df_abs <- bio4 %>% 
    filter(RCP %in% c('Current', !!RCP), 
           region == entire, PFT %in% pfts)# %>% 
    #filter(summary == 'median')
  
  # % change in biomass due to climate
  df_diff_cref <- pft5_bio_d3 %>% 
    filter(region == entire) %>% 
    filter(RCP == !!RCP, PFT %in% pfts)
  
  df_diff_gref <- pft5_d_wgcm2 %>% 
    filter(region == entire) %>% 
    filter(RCP %in% c('Current', !!RCP), 
           PFT %in% pfts)
  
  g1 <- box_abs_diff(df_abs = df_abs,
               df_diff_cref = df_diff_cref,
               df_diff_gref = df_diff_gref,
               y_abs = 'biomass',
               y_diff = 'bio_diff',
               ylab_abs = lab_bio0,
               ylab_diff_gref = ylab_diff_gref,
               ylab_diff_cref = ylab_diff_cref,
               scales_cref = 'free_y',
               scales_gref = 'fixed')
  
  # Pub qual for the main manuscript
  n_pft <- length(pfts)
  filename <- paste0("bio-abs-diff_weighted_entire_",
                     n_pft, "pft_boxplot_", RCP, '_', years, "_", 
                     run, v, ".png")
  
  ggsave(file.path("figures/biomass", filename),plot = g1,
      width = 7, height = 5*n_pft/3, dpi = 600)

})




# *testing -----------------------------------------------------------------

# this code can be deleted (used to confirm that these
# numbers correspond to the Core area figure)
# df <- sei2 %>% 
#   filter(region == 'Southern Great Basin',
#          years == 'Current', summary == 'median')
# 
# ggplot(df, aes(x = graze, y = SEI, weight = weight, fill = graze)) +
#   geom_violin() +
#   scale_fill_manual(values = cols_graze) +
#   geom_hline(yintercept = 0.431)
# 
# 
# l <- split(df, df$graze)[c(1, 3)]
# map(l, \(df) Hmisc::wtd.quantile(df$SEI, weights = df$weight,
#                                  probs = seq(0.5, 0.75, by = 0.05)))


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
                     subtitle = 'SEI',
                     strip.text.y = element_text(size = rel(0.7)))


g2 <- combine_grid_panels1(plots, remove_y_axis_labels = TRUE)
g3 <- g2 + g_sei + plot_layout(guides = 'collect')

png(paste0("figures/sei/q-sei_weighted_by-region_3pft_boxplot_", suffix, ".png"),
    width = 11, height = 10, units = 'in', res = 600)
g3&theme(legend.position = 'bottom')
dev.off()

# *utilization ----------------------------------------------------------------

plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(bio4, PFT == pft),
                y_string = 'utilization',
                ylab = lab_util0,
                subtitle = pft) +
    expand_limits(y = 0)
})

g2 <- combine_grid_panels1(plots, remove_y_axis_labels = FALSE)

png(paste0("figures/util/util_weighted_by-region_3pft_boxplot_", suffix, ".png"),
    width = 8, height = 10, units = 'in', res = 600)
g2
dev.off()


# *individuals ----------------------------------------------------------------

plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(bio4, PFT == pft),
                y_string = 'indivs',
                ylab = lab_indivs0,
                subtitle = pft) +
    expand_limits(y = 0)
})

g2 <- combine_grid_panels1(plots, remove_y_axis_labels = FALSE)

png(paste0("figures/indivs/indivs_weighted_by-region_3pft_boxplot_", suffix, ".png"),
    width = 8, height = 11, units = 'in', res = 600)
g2
dev.off()


# stacked bar charts ------------------------------------------------------

# * utilization -----------------------------------------------------------

g <- ggplot(bio_pcent1, aes(graze, util_pcent, fill = PFT)) +
  geom_col(position = 'stack') +
  facet_grid(region~rcp_year) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = 'bottom') +
  scale_fill_manual(values = cols_pft3) +
  labs(x = lab_graze,
       y = lab_util0b)


ggsave(
  paste0('figures/util/util-perc_bar_pft3_', suffix, '.png'),
  g, width = 4, height = 9
)


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
    df_factor() %>% 
    mutate(region = region2wafwa(region)) %>% 
    group_by(PFT, region) %>% 
    mutate(Q = bio2q(biomass, pft = unique(as.character(PFT)), 
                     region = unique(as.character(region)))) %>% 
    select(-row)
  
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


# clim vs bio by GCM ------------------------------------------------------
# show the 25th to 75 percentiles of biomass (y axis) and climate variables 
# (x axis)

linewidth <- 0.3
rcps <- unique(driver1$RCP)

for (rcp in rcps) { 
  for (reg in regions) {

    df <- driver2 %>% 
      filter(years == !!years,
             RCP == rcp,
             graze == 'Moderate',
             region == reg) %>% 
      mutate(GCM = factor(GCM, levels = names(cols_GCM1)))
    

    subtitle <- rcp_label(rcp, years, include_parenth = FALSE)
    if(reg != entire) {
      subtitle <- paste0(subtitle, ', ', reg)
    }
    
    # Get per-facet x-ranges (so geom_vline only plotted in som panels)
    ranges <- df |>
      dplyr::group_by(variable_clim, variable_bio) |>
      dplyr::summarise(xmin = min(c(p25_clim, median_clim, p75_clim), na.rm = TRUE),
                       xmax = max(c(p25_clim, median_clim, p75_clim), na.rm = TRUE),
                       .groups = "drop") |>
      dplyr::filter(xmin <= 0 & xmax >= 0)  # keep only panels that include 0
  
    g <- ggplot(df) +
      geom_hline(yintercept = 0, alpha = 0.3, linetype = 2, 
                 linewidth = linewidth) +
      # add geom_vline only for panels that include 0
      geom_vline(
        data = ranges,
        mapping = aes(xintercept = 0),
        alpha = 0.3, linetype = 2, linewidth = linewidth
      ) +
      geom_segment(aes(x = median_clim, y = p25_bio, yend = p75_bio,
                       color = GCM), linewidth = linewidth) +
      geom_segment(aes(x = p25_clim, xend = p75_clim, y = median_bio,
                       color = GCM), linewidth = linewidth) +
      geom_point(aes(x = median_clim, y = median_bio, color = GCM, shape = GCM)) +
      #geom_smooth(aes(x = median_clim, y = median_bio), se = FALSE, method = lm) +
      facet_grid(variable_bio ~ variable_clim, scales = 'free', 
                 switch = 'x', 
                 labeller = labeller(variable_clim = driver_labeller(delta = TRUE))) +
      scale_color_manual(values = cols_GCM1) +
      scale_shape_manual(values = shapes_GCM1) +
      theme(strip.placement.x = 'outside') +
      labs(x = NULL, y = lab_bio1,
           subtitle = subtitle)

    filename <- paste0('figures/biomass/by_gcm/bio_vs_clim_crossplot_',
                        vr, '_', words2abbrev(reg), '_', rcp, '_', 
                       years, '_', runv, '.png')
    ggsave(filename, plot = g,
           width = 7, height = 5.5, dpi = 600)
  }
}



# clim vs clim by GCM ------------------------------------------------------
# crossplots comparing interquartile ranges of delta climate
# variables, separately by GCM


vars <- tibble(var1 = c('MAT', 'MAT', 'MAP'),
               var2 = c('MAP', 'PSP', 'PSP'))
for (rcp in rcps) { 
  for (reg in regions) {
    
    df <- driver2 %>% 
      filter(years == !!years,
             RCP == rcp,
             graze == 'Moderate',
             region == reg) %>% 
      mutate(GCM = factor(GCM, levels = names(cols_GCM1)))
    
    
    subtitle <- rcp_label(rcp, years, include_parenth = FALSE)
    if(reg != entire) {
      subtitle <- paste0(subtitle, ', ', reg)
    }
    
    df_wide <- df %>% 
      select(-matches('_bio')) %>% 
      filter(!duplicated(.)) %>% 
      pivot_wider(
        names_from = variable_clim,
        values_from = c(median_clim, p25_clim, p75_clim, mean_clim),
        names_glue = "{variable_clim}_{.value}"
      ) %>% 
      rename_with(.fn = \(x) str_replace(x, '_clim$', ''))
  
    plots <- pmap(vars, function(var1, var2) {
      crossplot_1panel(df_wide = df_wide, var1 = var1, var2 = var2)
    })
    
    p <- wrap_plots(c(plots, list(patchwork::guide_area())),
                 nrow = 2, guides = 'collect')
    
    p2 <- p & guides(color=guide_legend(ncol = 2),
             linetype= guide_legend(ncol = 2))
    
    p3 <- p2 + plot_annotation(subtitle = subtitle)
    filename <- paste0('figures/climate/by_gcm/clim_vs_clim_crossplot_',
                       vr, '_', words2abbrev(reg), '_', rcp, '_',
                       years, '_', runv, '.png')
    ggsave(filename, plot = p3,
           width = 6, height = 5.5, dpi = 600)
  }
}

