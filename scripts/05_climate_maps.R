# Martin Holdrege

# Started: March 15, 2024

# Purpose: Map and examine interpolated climate data from STEPWAT2
# (i.e. maps showing what climate is being represented across space,
# by stepwat results, maps also show change in climate)


# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(stars)
library(patchwork)
library(RColorBrewer)
source("src/general_functions.R")
source("src/fig_functions.R")
source('src/mapping_functions.R')

# params ------------------------------------------------------------------

# date in input file names
date_in <- "20230919"
date <- '20240515' # appended to output file names

# Read in data ------------------------------------------------------------

# selecting which rasters to load
# interpolated rasters of stepwat data
path_r <- "data_processed/interpolated_rasters"


# * median climate vars (across GCMs) -------------------------------------

r1 <- rast(paste0(path_r, '/climate_median_across_GCMs_', date_in, '.tif'))

into <- c("type", "description", "RCP", "years")
info1 <- create_rast_info(r1, into = into, run_regex = "^")

r2 <- r1[[info1$id]]

# *median delta climate ---------------------------------------------------

# raw difference relative to historical 
rdiff1 <- rast(paste0(path_r, '/climate_rdiff-cref_median_', date_in, '.tif'))

info_rdiff1 <- create_rast_info(rdiff1, into = into, run_regex = '^')

# Figures -----------------------------------------------------------------

# combining difference and absolute biomass
info_c1 <- bind_rows(info1, info_rdiff1) %>% 
  # for exploration just looking at subset of scenarios
  filter_clim_extremes() %>% 
  mutate(type = fct_rev(factor(type))) %>% 
  arrange(run, RCP, type, years)

info_c_l <- split(info_c1, info_c1$type)

r_c1 <- c(r2, rdiff1) # combined raster


# maps of climate and raw difference --------------------------------------
# One big map on the left of historical climate, and 4 panels on the right
# showing 2 future time periods and for each time period 1 map shows
# climate the other shows delta climat


title_diff <- "\u0394 # fires/century" # delta

labels <- list(
  'MAT' = c('title' = 'MAT', 
            'unit' = 'MAT (C)'),
  'MAP' = c('title' = 'MAP', 
            'unit' = 'MAP (mm)'),
  'psp' = c('title' = 'Proportion summer precipitation (Jun-Aug)', 
            'unit' = 'PSP'),
  'PTcor' = c('title' = 'Seasonality (monthly Precip and Temp correlation)', 
            'unit' = 'PTcor')
)

clim_cols <- list(
  'MAT' = rev(brewer.pal(9, 'RdYlBu')),
  'MAP' = brewer.pal(9, 'YlGnBu'),
  'psp' = brewer.pal(9, 'RdYlBu'),
  'PTcor' = brewer.pal(9, 'RdYlBu')
)
c <- cols_map_bio_d
diff_cols <- list(MAT = rev(c), MAP = c, psp = c, PTcor = c)

# using cair_pdf so 'delta' symbol printed
cairo_pdf(paste0("figures/climate/climate-rdiff-cref_", date, ".pdf"),
          width = 12, height = 7, onefile = TRUE)
for(df in info_c_l){
  type <- as.character(unique(df$type))
  stopifnot(length(type) == 1)
  print(type)
  
  # fig labels
  unit <- labels[[type]]['unit']
  title_diff <- paste("\u0394", unit)

  # raster layer id's
  clim_id <- df$id[df$description == 'climate']
  diff_id <- df$id[df$description == 'rdiff']
  stopifnot(length(clim_id) == 3)
  stopifnot(length(diff_id) == 2)
  
  # determining axis limits
  
  range_c <- range(as.numeric(minmax(r_c1[[clim_id]])))
  # expanding range so don't have gray patches at the limit
  range_c <- range_c + c(-0.0001, 0.0001) 
  
  m <- max(abs(range(as.numeric(minmax(r_c1[[diff_id]])))))
  range_d <- c(-m, m) # for colour gradient b/ can't set midpoint
  
  # plots of climat variable
  maps_clim1 <- map(clim_id, function(id) {
    
    d <- create_rast_info(id, into = into, run_regex = '^')
    
    plot_map_inset(r = r_c1[[id]],
                   colors = clim_cols[[type]],
                   tag_label = paste(unit, rcp_label(d$RCP, d$years)),
                   limits = range_c,
                   scale_name = unit)
    
  })
  
  # maps of biomass difference (for each time period)
  maps_diff1 <- map(diff_id, function(id) {
    
    d <- create_rast_info(id, into = into, run_regex = '^')
    
    plot_map_inset(r = r_c1[[id]],
                   colors = diff_cols[[type]],
                   tag_label = paste(title_diff, rcp_label(d$RCP, d$years)),
                   limits = range_d,
                   scale_name = title_diff)
    
  })
  
  # combining the plots
  p <- maps_clim1[[1]] + ((maps_clim1[[2]] + maps_diff1[[1]])/(maps_clim1[[3]] + maps_diff1[[2]])) 
  
  p2 <- (p + plot_layout(guides = 'collect'))&
    theme(legend.position = 'bottom')
  
  p3 <- p2+
    patchwork::plot_annotation(
      caption = 'CMIP5 climate data interpolated from 200 STEPWAT2 sites \n median across GCMs shown',
      subtitle = labels[[type]][['title']])
  print(p3)
}
dev.off()

