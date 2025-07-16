# Purpose: create map of study area showing ecoregions, these are the regions
# some the data analyses use for data summarization


# params ------------------------------------------------------------------

source('src/params.R')
v <- v_interp
# version of the ecoregions used (i.e. how the study area is devided up)
# vr <- 'r1.1' # 'r1.0' is normal 4 regions # this object created in params.R
wafwa_only <- FALSE # only plot the wafwa ecoregions

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source('src/mapping_functions.R')
source('src/fig_params.R')


# read in data ------------------------------------------------------------

r_eco <- load_wafwa_ecoregions_raster(wafwa_only = wafwa_only,
                                      v = vr)

# map ---------------------------------------------------------------------
n <- nrow(levels(r_eco)[[1]])
rows <- if(n <=4) 2 else 3

g <- plot_map2(r_eco) +
  scale_fill_manual(values = get_cols_ecoregion(vr),
                    na.value = 'transparent',
                    na.translate = FALSE,
                    name = NULL) +
  theme(legend.position = 'bottom',
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )+
  guides(fill = guide_legend(nrow = rows))

suffix <- if(wafwa_only) {
  ""
} else {
  paste0('_', vr)
}

name <- paste0("figures/ecoregions_", v, suffix)

ggsave(filename = paste0(name, '.png'), plot = g, width = 6, height = 7, dpi = 900)

saveRDS(g, paste0(name, '.rds'))
