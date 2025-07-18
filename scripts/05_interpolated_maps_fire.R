# Martin Holdrege

# Started: March 3, 2025

# Purpose: Map and examine interpolated fire probability
# code adapted from SEI/scripts/sw_results/01_examine_sw_results_fire.R


# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
library(stars)
library(patchwork)
source("src/general_functions.R")
source("src/fig_params.R")
source("src/mapping_functions.R")
source('src/fig_functions.R')
source('src/params.R')
# params ------------------------------------------------------------------

runs <- opt$run
names(runs) <- runs
v_out <- "v1" # version appended to output
years <- opt$years
yr_lab <- opt$years
graze_levels <- c("grazL" = "Light")
ref_graze <- 'Light' # reference graze level for within gcm comparisons

# Read in data ------------------------------------------------------------

# selecting which rasters to load
# interpolated rasters of stepwat data
path_r <- "data_processed/interpolated_rasters"


# * median fire probability (across GCMs) --------------------------------------------

paths <- list.files(path_r, full.names = TRUE,
                    pattern = 'fire-prob_future_summary_across_GCMs.tif') %>% 
  str_subset(pattern = paste("(", runs, ")", collapse = "|", sep = "")) 

r1 <- rast(paths)

into0 <- c("type", "RCP", "years", 
          "graze", 'summary') 
into <- into0[-length(into0)]
info0 <- create_rast_info(r1, into = into0)

info1a <- info0 %>% 
  filter(summary == 'median') 
  

r2 <- r1[[info1a$id]]
names(r2) <- str_replace(names(r2), "_median$", "")

info1 <- info1a %>% 
  mutate(id = str_replace(id, "_median$", "")) %>% 
  filter(graze %in% graze_levels,)

# *median delta fire-prob ---------------------------------------------------

# raw difference relative to historical climate conditions
paths_rdiff <- list.files(path_r, full.names = TRUE,
                          pattern = 'fire-prob-rdiff-cref_median') %>% 
  str_subset(pattern = paste("(", runs, ")", collapse = "|", sep = "")) 

rdiff1 <- rast(paths_rdiff)

info_rdiff0 <- create_rast_info(rdiff1, into = into)%>% 
  filter_clim_extremes(years = years)

info_rdiff1 <- info_rdiff0 %>% 
  filter(graze %in% graze_levels) 
  
# * median delta fire prob wgcm ------------------------------------------

paths_wgcm <- map(runs, function(run) {
  p <- list.files(path_r, full.names = TRUE,
             pattern = paste0(run, "_fire-prob-wgcmDiff-", ref_graze, "_median"))
  stopifnot(length(p) == 1)
  p
})

wgcm_l1 <- map(paths_wgcm, rast)


# figures--20 panel -------------------------------------------------------

# 20 panel figures ---------------------------------------------------------

# preparing args
r <- c(r1, rdiff1)
# r <- spatSample(r, 100, method = 'regular', as.raster = TRUE) # for testing

info_20panel <- info1a %>% 
  bind_rows(info_rdiff0) %>% 
  filter(years %in% c('Current', !!years)) 

g <- plot_map_20panel(
  r = r,
  info = info_20panel,
  type_absolute = 'fire-prob',
  type_diff =  'fire-prob-rdiff-cref',
  title = 'Fire probability',
  name4absolute = lab_firep0,
  name4diff = NULL,
  legend_title_absolute = lab_firep0,
  legend_title_diff = lab_firep1,
  palette_absolute = cols_firep,
  palette_diff = rev(cols_map_bio_d2)
  )


filename <- paste0('figures/fire/maps/20panel_fire-prob_', run, yr_lab, ".png")
png20panel(filename)
print(g)
dev.off()
  

# * percent change --------------------------------------------------------

info_perc <- info1a %>% 
  filter(RCP == 'Current') %>% 
  select(-years, -RCP, -type) %>% 
  right_join(info_rdiff0, by = c('run', 'graze', 'run2'),
             suffix = c('_Current', '_diff'))

# percent change
r_diff_perc <- r[[info_perc$id_diff]]/r[[info_perc$id_Current]]*100

g <- plot_map_20panel(
  r = c(r1, r_diff_perc),
  info = info_20panel,
  type_absolute = 'fire-prob',
  type_diff =  'fire-prob-rdiff-cref',
  title = 'Fire probability',
  name4absolute = lab_firep0,
  name4diff = '% Change',
  legend_title_absolute = lab_firep0,
  legend_title_diff = '% Change in # fires/century',
  palette_absolute = cols_firep,
  palette_diff = rev(cols_map_bio_d2),
  lims_diff = c(-100, 400),
  midpoint_diff = 0
)


filename <- paste0('figures/fire/maps/20panel_fire-prob-perc_', run, yr_lab, ".png")
png20panel(filename)
print(g)
dev.off()

# Figures -----------------------------------------------------------------

# combining difference and absolute biomass
info_c1 <- bind_rows(info1, info_rdiff1) %>% 
  filter_clim_extremes(years = years) %>% 
  mutate(type = fct_rev(factor(type))) %>% 
  arrange(run, RCP, type, years)

r_c1 <- c(r2, rdiff1) # combined raster


# maps of biomass and raw difference --------------------------------------
# One big map on the left of historical biomass, and 4 panels on the left
# showing 2 future time periods and for each time period 1 map shows
# biomass the other shows delta biomass
info_c_l <- info_c1 %>% 
  group_by(run2) %>% 
  group_split() # split into list

# ranges for axes
# ids for all runs for this PFT
fire_id_all <- info_c1 %>% 
  filter(type == 'fire-prob') %>% 
  pull(id)

diff_id_all <- info_c1 %>% 
  filter(type != 'fire-prob') %>% 
  pull(id)

# want to get range across run types so that figures will
# have comparable colors across runs for given pft
max_f <- max(as.numeric(minmax(r_c1[[fire_id_all]])))
# adding small amount b/ I think rounding is causing this max to be too low (and some values not display)
range_f <- c(0, max_f + 0.1)
range_d0 <- range(as.numeric(minmax(r_c1[[diff_id_all]])))
m <- max(abs(range_d0)) # for colour gradient b/ can't sent midpoint
range_d <- c(-m, m)

title_diff <- "\u0394 # fires/century" # delta

# using cair_pdf so 'delta' symbol printed
cairo_pdf(paste0("figures/fire/maps/fire-prob-rdiff-cref_", v_out, yr_lab,
                 ".pdf"),
          width = 12, height = 7, onefile = TRUE)
for(df in info_c_l){
  print(df$id[1])
  fire_id <- df$id[df$type == 'fire-prob']
  diff_id <- df$id[df$type != 'fire-prob']
  stopifnot(length(fire_id) == 3)
  stopifnot(length(diff_id) == 2)
  
  # plots of biomass
  maps_fire1 <- map(fire_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = r_c1[[id]],
                   colors = cols_firep,
                   tag_label = paste("Fire probability", rcp_label(d$RCP, d$years)),
                   limits = range_f,
                   scale_name = lab_firep0)
    
  })
  
  # maps of biomass difference (for each time period)
  maps_diff1 <- map(diff_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = r_c1[[id]],
                   colors = rev(cols_map_bio_d),
                   tag_label = paste(title_diff, rcp_label(d$RCP, d$years)),
                   limits = range_d,
                   scale_name = lab_firep1)
    
  })
  
  # combining the plots
  p <- maps_fire1[[1]] + ((maps_fire1[[2]] + maps_diff1[[1]])/(maps_fire1[[3]] + maps_diff1[[2]])) 
  
  p2 <- (p + plot_layout(guides = 'collect'))&
    theme(legend.position = 'bottom')
  
  p3 <- p2+
    patchwork::plot_annotation(caption = paste('simulation settings:',  df$run2))
  print(p3)
  
  # create same plot but with right figures showing % change instead 
  # of absolute difference
  
  perc_diff <-(r_c1[[diff_id]]/r_c1[[fire_id[1]]])*100
  tmp <- max(abs(as.vector(minmax(perc_diff))))
  if(tmp>100) {
    range_perc <- c(-100, min(300, tmp))
    cap_perc <- '(% Change range truncated)'
  } else {
    range_perc <- c(-tmp, tmp)
    cap_perc <- ""
  }
  
  maps_perc1 <- map(diff_id, function(id) {
    
    d <- create_rast_info(id, into = into)
    
    plot_map_inset(r = perc_diff[[id]],
                   colors = rev(cols_map_bio_d),
                   tag_label = paste('% change', rcp_label(d$RCP, d$years)),
                   limits = range_perc,
                   midpoint = 0,
                   scale_name = '% change')
    
  })
  # combining the plots
  p <- maps_fire1[[1]] + ((maps_fire1[[2]] + maps_perc1[[1]])/(maps_fire1[[3]] + maps_perc1[[2]])) 
  
  p2 <- (p + plot_layout(guides = 'collect'))&
    theme(legend.position = 'bottom')
  
  p3 <- p2+
    patchwork::plot_annotation(caption = paste(df$run2[1], cap_perc))
  print(p3)
}


dev.off()



# comparing grazing levels ------------------------------------------------

for (run in runs) {
  r <- wgcm_l1[[run]]
  m <- max(abs(unlist(minmax(r))))
  range_wgcm <- c(-m, m)
  info_wgcm1 <- create_rast_info(r,
                                 into = c("type", 'RCP', 'years', 'graze')) %>% 
    mutate(id_noGraze = str_replace(id, '_[[:alpha:]]+$', '')) %>% 
    filter_clim_extremes() %>% 
    arrange(RCP, years, graze)
  
  info_wgcm_l <- split(info_wgcm1, info_wgcm1$id_noGraze)
  # iterating across scenarios (1 page per scenario)
  pages <- map(info_wgcm_l, function(info) {

    # iterating across grazin levels
    maps_delta <- map(info$id, function(id) {
      tag_label <- paste0("Delta fire probability (", 
                          info$graze[info$id == id], "\n relative to ", 
                          ref_graze, " grazing)")
      plot_map_inset(r = r[[id]],
                     colors = rev(cols_map_bio_d),
                     tag_label = tag_label,
                     limits = range_wgcm,
                     scale_name = lab_firep1)
    })
    
    id <- paste0(unique(info$id_noGraze), '_', ref_graze)
    map_ref <-plot_map_inset(r = r2[[id]],
                             colors = cols_firep,
                             tag_label = paste0("Fire probability (", ref_graze, 
                                                " grazing)"),
                             limits = range_f,
                             scale_name = lab_firep0)
    
    
    comb1 <- patchwork::wrap_plots(c(list(map_ref), maps_delta),
                                   ncol = 2, guides = 'collect')
    
    comb2 <- comb1 +
      patchwork::plot_annotation(
        subtitle = paste0('Change in fire probability due to grazing', 
                          ' (comparisons within climate scenario)\n',
                          rcp_label(unique(info$RCP), unique(info$years))))
    comb2
  })
  
  pdf(paste0("figures/fire/maps/", run, yr_lab, "_delta-prob_wgcm-", ref_graze, '_', 
             v_out, '.pdf'),
      width = 10, height = 10)
    map(pages, print)
  dev.off()
  
}
