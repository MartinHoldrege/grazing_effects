# Martin Holdrege

# Script started Jan 7, 2022


# Purpose: Create maps of biomass and change in biomass across the west
# based on STEPWAT2 output for 200 sites (for each plant functional
# type and climate and grazing scenario) that has been upscaled across
# the west. The rasters of interpolated data are created
# in the 03_interpolate.R script

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

run <- 'fire1_eind1_c4grass1_co20_2502'
graze_levels <- c("grazL" = "Light")

# PFTs for which to plot 
PFTs2plot <- c("Sagebrush", "Aherb", "Pherb", 'C4Pgrass', "C3Pgrass", "Pforb", 
               "Cheatgrass", "Aforb")

v <- 'v1' # for appending to output file names

# Read in data ------------------------------------------------------------

# selecting which rasters to load
# interpolated rasters of stepwat data
path_r <- "data_processed/interpolated_rasters"

# * median biomass (across GCMs) --------------------------------------------
list.files(path_r, "bio_future_median_across_GCMs")
path <- paste0(path_r, "/", run, "_bio_future_median_across_GCMs.tif")
r1 <- rast(path)

into <- c("PFT", "type", "RCP", "years", 
          "graze")
info0 <- create_rast_info(r1, into = into)

info1 <- info0 %>% 
  filter(graze %in% graze_levels,
         PFT %in% PFTs2plot)

r3 <- r1[[info1$id]]


# *median delta biomass ---------------------------------------------------

# raw difference relative to historical climate conditions
path_rdiff <- paste0(path_r, "/", run, "_bio-rdiff-cref_median.tif")

rdiff1 <- rast(path_rdiff)

info_rdiff1 <- create_rast_info(rdiff1, into = into)%>% 
  filter(graze %in% graze_levels,
         PFT %in% PFTs2plot) %>% 
  # for exploratory reasons just looking at the most and least
  # extreme scenarios
  filter_clim_extremes()

rdiff2 <- rdiff1


# Figures -----------------------------------------------------------------

# combining difference and absolute biomass
info_c1 <- bind_rows(info1, info_rdiff1) %>% 
  filter_clim_extremes() %>% 
  filter(PFT %in% PFTs2plot) %>% 
  mutate(type = fct_rev(factor(type))) %>% 
  arrange(run, PFT, RCP, type, years)

r_c1 <- c(r3, rdiff2) # combined raster

# maps of biomass and raw difference --------------------------------------

# One big map on the left of historical biomass, and 4 panels on the left
# showing 2 future time periods and for each time period 1 map shows
# biomass the other shows delta biomass
info_c_l <- info_c1 %>% 
  group_by(PFT, run2) %>% 
  group_split() %>%  # split into list
  map(\(x) arrange(x, RCP, years))

pdf(paste0("figures/biomass_maps/bio-rdiff-cref_", v, "_", run, ".pdf"),
          width = 11, height = 7)
for(df in info_c_l){
  print(df$id[1])
  bio_id <- df$id[df$type == 'biomass']
  diff_id <- df$id[df$type != 'biomass']
  stopifnot(length(bio_id) == 3)
  stopifnot(length(diff_id) == 2)
  
  # ids for all runs for this PFT
  bio_id_all <- info_c1 %>% 
    filter(df$PFT[[1]] == PFT,
           type == "biomass") %>% 
    pull(id)
  
  stopifnot(is.factor(df$RCP)) # check so that RCP correctly sorted
  
  diff_id_all <- info_c1 %>% 
    filter(df$PFT[[1]] == PFT,
           type != "biomass") %>% 
    pull(id)
  
  # want to get range across run types so that figures will
  # have comparable colors across runs for given pft
  range_b <- range(as.numeric(minmax(r_c1[[bio_id_all]])))
  range_d <- range(as.numeric(minmax(r_c1[[diff_id_all]])))
  
  m <- max(abs(range_d)) # for colour gradient b/ can't sent midpoint
  title_diff <- "Delta biomass" # delta
  
  # plots of biomass
  maps_bio1 <- map(bio_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = r_c1[[id]],
                   colors = cols_map_bio(10),
                   tag_label = paste("Biomass", rcp_label(d$RCP, d$years)),
                   limits = range_b,
                   scale_name = lab_bio0)
    
  })
  
  # maps of biomass difference (for each time period)
  maps_diff1 <- map(diff_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = r_c1[[id]],
                   colors = cols_map_bio_d,
                   tag_label = paste(title_diff, rcp_label(d$RCP, d$years)),
                   limits = c(-m, m),
                   scale_name = lab_bio1)
    
  })
  
  # combining the plots
  p <- maps_bio1[[1]] + ((maps_bio1[[2]] + maps_diff1[[1]])/(maps_bio1[[3]] + maps_diff1[[2]])) 
  
  p2 <- (p + plot_layout(guides = 'collect'))&
    theme(legend.position = 'bottom')
  
  p3 <- p2+
    patchwork::plot_annotation(df$PFT[1],
                               caption = df$run2[1])
  print(p3)
  
  # create same plot but with right figures showing % change instead 
  # of absolute difference
  
  perc_diff <-(r_c1[[diff_id]]/r_c1[[bio_id[1]]])*100
  tmp <- max(abs(as.vector(minmax(perc_diff))))
  if(tmp>100) {
    tmp <- 100
  }
  range_perc <- c(-tmp, tmp)
  maps_perc1 <- map(diff_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = perc_diff[[id]],
                   colors = cols_map_bio_d,
                   tag_label = paste('% change', rcp_label(d$RCP, d$years)),
                   limits = range_perc,
                   scale_name = '% change')
    
  })
  # combining the plots
  p <- maps_bio1[[1]] + ((maps_bio1[[2]] + maps_perc1[[1]])/(maps_bio1[[3]] + maps_perc1[[2]])) 
  
  p2 <- (p + plot_layout(guides = 'collect'))&
    theme(legend.position = 'bottom')
  
  p3 <- p2+
    patchwork::plot_annotation(df$PFT[1],
                               caption = df$run2[1])
  print(p3)
  
}

dev.off()
