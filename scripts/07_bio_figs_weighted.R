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
# this script runs for all RCP/years
source('src/params.R')
vr <- opt$vr
run <- opt$run

runv <- opt$runv
ref_graze <- opt$ref_graze
q_curve_fig <- FALSE # create figure showing the 'biomass' q curves
v <- '_v2'

suffix <- paste0(runv, '_', vr)
pfts <- c("Sagebrush", 'Pherb', 'Aherb') # code won't work with any other pfts

regions <- region_factor(v = vr, return_levels = TRUE)
entire <- regions[1]
ref_graze <- opt$ref_graze
# read in data ------------------------------------------------------------
# created in "scripts/02_summarize_bio.R"
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
pft5_d_wgcm <- l$pft5_d_wgcm # summarized across gcms
pft5_d_wgcm_gcm <- l$pft5_d_wgcm_gcm # values by GCM
fire_d_wgcm_gcm <- l$fire_d_wgcm_gcm

# weighted means and percentiles of change biomass and climate variables
# created in "scripts/06_summarize_fire_drivers.R"
driver1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                  '_delta-fire-driver_by-ecoregion.csv'))

# combine site level and weights ------------------------------------------
w1 <- df_factor(w1)

bio3 <- bio2 %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  df_factor() 

yrs <- levels(bio2$years)
rcps <- levels(bio2$RCP)
yrs <- yrs[yrs != 'Current']
rcps <- rcps[rcps != 'Current']
# for looping over
df_rcp_yr <- expand_grid(RCP = rcps,
                         years = yrs)

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
  df_factor() %>% 
  filter(region != levels(region)[1])

# % change (climate effect)

prepare_bio_diff <- function(df) {
  df %>% 
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


join_vars <- c('run', 'years', 'RCP', 'graze',
               'site', 'GCM')

# GCM level delta biomass and delta fire responses to grazing
gref_diff_gcm1 <- pft5_d_wgcm_gcm %>% 
  filter(PFT %in% pfts) %>% 
  select(all_of(join_vars), PFT, bio_diff) %>% 
  left_join(
    select(fire_d_wgcm_gcm, all_of(join_vars), fire_prob_diff),
    by = join_vars
  ) %>% 
  rename(fire = fire_prob_diff,
         bio = bio_diff) %>% 
  left_join(w1, by = 'site',relationship = "many-to-many") %>% 
  group_by(run, years, RCP, graze, GCM, region, PFT) %>% 
  # weighted percentiles
  summarize_weighted(varname = c('bio', 'fire'),
                     .names = '{.col}_{.fn}') %>% 
  df_factor()
  
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

# boxplots ----------------------------------------------------------------

# *biomass ----------------------------------------------------------------

# biomass boxplots for all scenarios
plots <- map(pfts, function(pft) {
  weighted_box1(df = filter(bio4, PFT == pft),
                y_string = 'biomass',
                ylab = lab_bio0,
                subtitle = pft) +
    expand_limits(y = 0)
})

g2 <- combine_grid_panels1(plots, remove_y_axis_labels = FALSE)

png(paste0("figures/biomass/bio_weighted_by-region_3pft_boxplot_", suffix, 
           ".png"),
    width = 12, height = 11, units = 'in', res = 600)
g2
dev.off()


# absolute biomass for current and %change 
# pub qual, first biomass figure
pft5 <- c(pfts, 'C3Pgrass', 'C4Pgrass')
args <- expand_grid(
  pfts = list(pfts, pft5),
  RCP = c('RCP45', 'RCP85'),
  years = yrs
)

ylab_diff_gref = paste0('\u0394Biomass (%), relative to ', str_to_lower(ref_graze),
                        ' grazing')

# figure for entire region
pmap(args, function(pfts, RCP, years) {
  df_abs <- bio4 %>% 
    filter(RCP %in% c('Current', !!RCP), 
           years %in% c('Current', !!years),
           region == entire, PFT %in% pfts)# %>% 
    #filter(summary == 'median')
  
  # % change in biomass due to climate
  df_diff_cref <- pft5_bio_d3 %>% 
    filter(region == entire) %>% 
    filter(RCP == !!RCP, PFT %in% pfts,
           years == !!years)
  
  df_diff_gref <- pft5_d_wgcm2 %>% 
    filter(region == entire) %>% 
    filter(RCP %in% c('Current', !!RCP), 
           years %in% c('Current', !!years),
           PFT %in% pfts)
  
  g1 <- box_abs_diff(df_abs = df_abs,
               df_diff_cref = df_diff_cref,
               df_diff_gref = df_diff_gref,
               y_abs = 'biomass',
               y_diff = 'bio_diff',
               ylab_abs = lab_bio0,
               ylab_diff_gref = lab_bio2_gref,
               ylab_diff_cref = lab_bio2_cref,
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
    width = 12, height = 10, units = 'in', res = 600)
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
    width = 12, height = 11, units = 'in', res = 600)
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
  g, width = 6, height = 9
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
    mutate(region = region2wafwa(region, v = vr)) %>% 
    group_by(PFT, region) %>% 
    mutate(Q = bio2q(biomass, pft = unique(as.character(PFT)), 
                     region = unique(as.character(region)),
                     v = vr)) %>% 
    select(-row)
  
  g <- ggplot(q_range, aes(biomass, Q, color = region, linetype = region)) +
    geom_line() + 
    facet_wrap(~PFT, scales = 'free_x') +
    scale_color_manual(values = cols_wafwa_region) +
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



for (i in 1:nrow(df_rcp_yr)) { 
  rcp <- df_rcp_yr$RCP[i]
  years <- df_rcp_yr$years[i]
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
for (i in 1:nrow(df_rcp_yr)) { 
  rcp <- df_rcp_yr$RCP[i]
  years <- df_rcp_yr$years[i]
  for (reg in regions) {
    
    df <- driver2 %>% 
      filter(years == !!years,
             RCP == rcp,
             graze == ref_graze,
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

# gref bio vs fire by GCM ------------------------------------------------------
# crossplots showing delta fire against delta biomass, within a GCM
# comparing to moderate grazing

# i <- 1
# yr <- df_rcp_yr$years[i]
# rcp <- df_rcp_yr$RCP[i]
# pft <- pfts[1]
for(pft in pfts) {
  pmap(df_rcp_yr, function(RCP, years) {
    rcp <- RCP
    yr <- years
    df <- gref_diff_gcm1 %>% 
      filter(years %in% c('Current', yr),
             RCP %in% c('Current', rcp),
             PFT == pft) %>% 
      mutate(GCM = factor(GCM, 
                          levels = c('Current', names(cols_GCM1)),
                          labels = names(cols_GCM2)))
    
    subtitle <- paste0(pft, ', ', rcp_label(rcp, yr, include_parenth = FALSE))
    
    g <- crossplot_1panel(df, var1 = 'fire', var2 = 'bio',
                          colors = cols_GCM2,
                          shapes = shapes_GCM2) +
      # geom_smooth(aes(x = fire_median, y = bio_median), method = 'lm',
      #             se = FALSE, formula = ' y ~ x') +
      facet_grid(graze ~ region , scales = 'free_y') +
      labs(y = lab_bio2_gref,
           x = lab_firep1_gref,
           subtitle = subtitle)
    g <- tag(g, fig_letters)
    filename <- paste0('figures/biomass/by_gcm/bio_vs_fire_gref_crossplot_',
                       pft, "_", vr, '_by-region_', rcp, '_', 
                       yr, '_', runv, '.png')
    ggsave(filename, plot = g,
           width = 10, height = 6, dpi = 600)
  })
}
