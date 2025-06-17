# Martin Holdrege

# Started: March 21, 2025

# Purpose: Maps that examine what variables are most responsible for 
# the changes in projected fire probability. 


# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(stars)
library(patchwork)
source("src/general_functions.R")
source("src/fig_params.R")
source("src/mapping_functions.R")
source('src/fig_functions.R')

# params ------------------------------------------------------------------

source('src/params.R')
v_out <- "v1" # version appended to output
graze_levels <- c("grazL" = "Light", "grazVH" = "Very Heavy")

test_run <- FALSE
pfts <- c("Pherb", "Aherb")
# Read in data ------------------------------------------------------------

# selecting which rasters to load
# interpolated rasters of stepwat data
path_r <- "data_processed/interpolated_rasters"

# * median fire probability (across GCMs) --------------------------------------------

# *median delta fire-prob ---------------------------------------------------

# raw difference relative to historical climate conditions
path_rdiff <- file.path(path_r,
                        paste0(run, '_fire-prob-rdiff-cref_median.tif'))
rdiff1 <- rast(path_rdiff)

into <- c("type", "RCP", "years", 
          "graze")
info_rdiff0 <- create_rast_info(rdiff1, into = into)%>% 
  filter_clim_extremes()

info_rdiff1 <- info_rdiff0 %>% 
  filter(graze %in% graze_levels) 

# * drivers of fire change ----------------------------------------------

# dominant driver of change
paths_dom <- list.files(
  path = file.path(path_r, "fire_attrib", v_interp),
  pattern = paste0(run, ".*dom-driver"),
  full.names = TRUE
)

r_dom1 <- rast(paths_dom)

# change in fire probability due to a given predictor variable
paths_driver <- list.files(
  path = file.path(path_r, "fire_attrib", v_interp),
  pattern = paste0(run, ".*delta-"),
  full.names = TRUE
)
r_delta_driver1 <- rast(paths_driver)

# created in 02b_fire_attribution.R
pred_vars <- readRDS('data_processed/site_means/fire_dominant_drivers.RDS')$pred_vars

# * predictor variables ---------------------------------------------------

# **climate ---------------------------------------------------------------

# raw difference relative to historical 
# created in "scripts/04_interpolated_summarize_clim.R"
clim_diff1 <- rast(paste0(path_r, '/climate_rdiff-cref_median_', v_interp, '.tif'))

into_pred <- c("pred_var", "type", "RCP", "years")
info_clim_diff1 <- create_rast_info(clim_diff1, into = into_pred, run_regex = '^')

# **vegetation --------------------------------------------------------------

# raw difference relative to historical climate conditions
path_bio_diff <- paste0(path_r, "/", run, "_bio-rdiff-cref_median.tif")

bio_diff1 <- rast(path_bio_diff)

info_bio_diff1 <- create_rast_info(bio_diff1, 
                                    into =  c(into_pred, 'graze')) %>% 
  filter(pred_var %in% pfts)
bio_diff1 <- bio_diff1[[info_bio_diff1$id]]

if(test_run) {
  rdiff1 <- downsample(rdiff1)
  r_dom1 <- downsample(r_dom1)
  r_delta_driver1 <- downsample(r_delta_driver1)
  clim_diff1 <- downsample(clim_diff1)
  bio_diff1 <- downsample(bio_diff1)
}

# prepare rasters ---------------------------------------------------------

# adding categorical levels
pred_vars2 <- str_replace(pred_vars, "psp", "PSP")
tmp <- data.frame(ID = 1:length(pred_vars),
                  category = pred_vars2)
r_dom2 <- r_dom1

for(i in 1:nlyr(r_dom1)) {
  levels(r_dom2[[i]]) <- tmp
}
names(r_dom2) <- names(r_dom1)

# raster info -------------------------------------------------------------

into_attrib <- c("type", "RCP", "years", "graze", "summary_stat")

info_dom1 <- create_rast_info(r_dom2, into = into_attrib) %>% 
  filter_clim_extremes() %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
  arrange(RCP, graze)

info_delta1 <- create_rast_info(r_delta_driver1, into = into_attrib) %>% 
  # incorrect file naming using, so manually parsing
  mutate(pred_var = str_extract(type, "(?<=fire-delta-)[[:alpha:]]+"),
         type = str_replace(type, pred_var, "pred")) %>% 
  dplyr::select(-summary_stat)

# actual change in fire probability 
info_rdiff1 <- create_rast_info(rdiff1,
                               into = c("type", "RCP", "years", "graze")) %>% 
  mutate(pred_var = 'total')

info_delta2 <- bind_rows(info_delta1, info_rdiff1)%>% 
  filter_clim_extremes() %>% 
  filter(graze %in% graze_levels) %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE),
         pred_var = factor(pred_var, levels = c('total', pred_vars),
                           labels = c("total", pred_vars2)))

r_delta2 <- c(r_delta_driver1, rdiff1)[[info_delta2$id]]


# * delta pred vars -------------------------------------------------------

info_pred_diff1 <- bind_rows(info_bio_diff1, info_clim_diff1) %>% 
  filter(pred_var %in% pred_vars,
         graze %in% graze_levels | is.na(graze)) %>% 
  mutate(pred_var = driver2factor(pred_var)) %>% 
  filter_clim_extremes()

pred_delta1 <- c(clim_diff1, bio_diff1)[[info_pred_diff1$id]]

# Maps of delta fire prob -------------------------------------------------

# top right map shows absolute difference in fire prob,
# the remaining maps show change in probability due to a given 

m <- minmax(r_delta2) %>% 
  unlist() %>% 
  abs() %>% 
  max()
range_delta <- c(-m, m)

info_delta3 <- info_delta2 %>% 
  mutate(plot = vector(mode = "list", length = n())) %>% 
  arrange(RCP, years, graze, pred_var) %>% 
  group_by(rcp_year, graze) 
  
info_delta_l <- group_split(info_delta3)

pdf(paste0("figures/fire_attribution/maps/fire-delta-driver_", v_out, "_", run, ".pdf"),
          width = 8, height = 7)  

for(df in info_delta_l){

  # plots of biomass
  maps_delta <- pmap(df[c('id', 'pred_var')], function(id, pred_var) {
    
    if(pred_var == 'total') {
      tag <- "Total delta fire probability"
      limits <- range_delta
      scale_name <- expression(Delta ~ "# fires/century (Total)")
      colors <- rev(cols_map_bio_d)
    } else {
      tag <- paste("Change due to", pred_var)
      limits <- range_delta*0.4
      scale_name <- expression(atop(Delta ~ "# fires/century", "(due to change in predictor)"))
      colors <- rev(cols_map_bio_d2)
    }
    plot_map_inset(r = r_delta2[[id]],
                   colors = colors,
                   tag_label = tag,
                   limits = limits,
                   scale_name = scale_name)
    
  })
  
  # putting into list for use in other plots
  info_delta3$plot[info_delta3$id %in% df$id] <- maps_delta
  
  p <- wrap_plots(maps_delta, ncol = 3) +
    plot_layout(guides = 'collect') +
    patchwork::plot_annotation(
      subtitle = paste0('Attribution of changes in fire probability ',
                        '(',df$graze[1], " grazing, ",
                       df$rcp_year[1], ")"),
      caption = paste(run, "\n Change in fire probability due to a predictor",
                      "variable calculated by holding that variable at historical",
                      "levels", 
                      "\nand calculating difference from future"))
  
  p2 <- p&theme(legend.position = 'bottom')
  
  print(p2)
  }
dev.off()

# maps of primary driver of change ----------------------------------------

info_dom_l <- info_dom1 %>% 
  group_by(rcp_year) %>% 
  group_split()

pdf(paste0("figures/fire_attribution/maps/fire-dom-driver_", v_out, "_", run, ".pdf"),
    width = 6, height = 7)  
for (df in info_dom_l) {
  maps <- pmap(df[c('id', 'graze')], function(id, graze) {
    
    r <- r_dom2[[id]]
        lvl <- levels(r)[[1]]
    s <- st_as_stars(r)

    s[[1]] <- factor(s[[1]], levels = c(lvl$ID), labels = c(lvl$category))
    
    plot_map2(
      r = s,
      add_coords = TRUE
    ) +
      newRR3::add_tag_as_label(paste(graze, 'grazing')) +
      scale_fill_manual(
        name = "Primary driver of change",
        values = cols_pred_vars,
        na.value = 'transparent',
        na.translate = FALSE,
        drop = FALSE
      )
  })
  
  p <- wrap_plots(maps, ncol = 2) +
    plot_layout(guides = 'collect') +
    patchwork::plot_annotation(
      subtitle = paste0('Primary driver of changes in fire probability',
                        '\nunder a fixed grazing level ',
                        '(', df$rcp_year[1], ")"),
      caption = run)
  
  p2 <- p&theme(legend.position = 'bottom')
  
  print(p2)
}
dev.off()


# delta fire and delta pred var maps --------------------------------------
# side by side, change in fire due to a variable and change in the variable

info_pred_diff2 <- info_pred_diff1 %>% 
  group_by(pred_var, RCP) %>% 
  nest() %>% 
  mutate(
    lims = map(data, function(x) {
      r <- pred_delta1[[x$id]]
      range_raster(r, absolute = TRUE)
    })
  ) %>% 
  select(-data) %>% 
  right_join(info_pred_diff1, join_by(pred_var, RCP))

# create the plots
info_pred_diff2$plot_pred <- pmap(info_pred_diff2[c("pred_var", "id", "lims")], 
                              function(pred_var, id, lims) {
                                
                              tag <- paste("Change in", pred_var)
                              scale_name <-  labs_pred1[[as.character(pred_var)]]
                              p <- plot_map_inset(r = pred_delta1[[id]],
                                             colors = rev(cols_map_bio_d2),
                                             tag_label = tag,
                                             limits = lims,
                                             scale_name = scale_name
                                             ) 
                              p&theme(legend.position = 'left')
                                
                              })

# repeating climate plots for each grazing level
info_pred_diff3 <- info_pred_diff2 %>% 
  filter(!is.na(graze))
info_delta4 <-  info_delta3 %>% 
  filter(!is.na(graze))

tmp1 <- info_pred_diff2 %>% 
  filter(is.na(graze))
tmp2 <- info_delta3 %>% 
  filter(is.na(graze))

for(lvl in graze_levels) {
  tmp1$graze <- lvl
  tmp2$graze <- lvl
  info_pred_diff3 <- bind_rows(info_pred_diff3, tmp1)
  info_delta4 <- bind_rows(info_delta4, tmp2)
}


info_plots_l <- info_delta4 %>%
  ungroup() %>% 
  select(RCP, years, graze, pred_var, plot) %>% 
  right_join(info_pred_diff3, by = join_by(RCP, years, graze, pred_var)) %>% 
  arrange(RCP, years, graze, pred_var) %>% 
  group_by(RCP, years, graze) %>% 
  group_split()

map(info_plots_l, function(df) {
  left <- wrap_plots(df$plot_pred) +
    patchwork::plot_layout(ncol = 1) +
    plot_annotation(subtitle = "Change in driver of wildfire")

  right <- wrap_plots(df$plot) +
    patchwork::plot_layout(ncol = 1, guides = 'collect') +
    plot_annotation(subtitle = "Change in fire due\nto change in driver")
  right <- right & theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title.position = "top"
    ) 
  g <- wrap_plots(left, right, ncol = 2, widths = c(1, 1))
  graz <- names(graze_levels[graze_levels == df$graze[1]])
  filename <- paste0("figures/fire_attribution/maps/delta-driver-delta-fire_",
                     df$RCP[1], "_", graz, "_", v_out, "_", run, ".png")
  
  ggsave(filename, plot = g,
          width = 6, height = 13, dpi = 600)
})

