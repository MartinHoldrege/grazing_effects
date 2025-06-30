# Purpose: create map of study area showing ecoregions, these are the regions
# some the data analyses use for data summarization


# params ------------------------------------------------------------------

v <- 'v4' # interpolation version
wafwa_only <- FALSE # only plot the wafwa ecoregions

# dependencies ------------------------------------------------------------

library(terra)
library(tidyverse)
source('src/mapping_functions.R')
source('src/fig_params.R')


# read in data ------------------------------------------------------------

r_eco <- load_wafwa_ecoregions_raster(wafwa_only = wafwa_only)
eco1 <- load_wafwa_ecoregions()
levels(eco1$ecoregion)

levels(r_eco)
# map ---------------------------------------------------------------------

g <- plot_map2(r_eco) +
  scale_fill_manual(values = cols_ecoregion,
                    na.value = 'transparent',
                    na.translate = FALSE,
                    name = NULL) +
  theme(legend.position = 'bottom')+
  guides(fill = guide_legend(nrow = 2))

suffix <- if(wafwa_only) {
  ""
} else {
  '_all'
}

name <- paste0("figures/ecoregions_", v, suffix)
png(paste0(name,  '.png'),
    width = 6, height = 7, units = 'in', res = 900)
g
dev.off()
saveRDS(g, paste0(name, '.rds'))
