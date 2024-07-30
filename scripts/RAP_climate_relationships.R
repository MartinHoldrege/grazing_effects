# Purpose: Examine relationship between RAP/RCMAP cover and climate variables
# to then compare to these same relationships based on STEPWAT2 output (i.e
# to see how STEPWAT simulations are performing). Here using smoothed 
# RAP/RCMAP data to get a average 'potential' cover in a given area. 

# Author: martin holdrege
# Started: July 29, 2024

# Dependencies --------------------------------------------------------------
library(terra)
library(tidyverse)
source('src/fig_functions.R')
source('../cheatgrass_fire/src/general_functions.R') # for making quantile plots
source('src/mapping_functions.R')
theme_set(theme_custom1())

# Read in data ------------------------------------------------------------

# cell numbers that correspond to the following csv
cellnumbers <- rast("data_processed/interpolation_data/cellnumbers.tif")

# csv the contains daymet normals, as well as data to construct multiple study areas
clim1 <- read_csv("data_processed/interpolation_data/clim_for_interp_daymet-v4_1991-2020.csv",
                  show_col_types = FALSE); # parsing fails for site_id col, which is fine

# contains average proportion summer precipitation (the average of the psp
# predictor variable used for the fire model)
psp1 <- rast("../cheatgrass_fire/data_processed/pred-fire-clim-veg_avg-across-yrs_annf3_A-P_entire.tif")$prcpPropSum

# Cover of sage, perennials and annuals, smoothed within 707 m. with multiple summaries taken
# the first number is the temporal percentile, and the second is the spatial percentile
cov1 <- rast("../SEI/data_processed/cover/cover_rap-rcmap_1986_2021_1000m_707msmooth_20230905.tif")

# create mask -------------------------------------------------------------

clim2 <- clim1 %>% 
  # 1 km grid cells containing > 50% of 30 km cells that have at least 1%
  # sage
  filter(prop_sage_gt_1 > 0.5) %>% 
  mutate(study_area = 1) 

# "study area"
study_area <- fill_raster(clim2[ , c('cellnumber', 'study_area')], cellnumbers)

# align layers ------------------------------------------------------------

psp2 <- crop(psp1, study_area) %>% 
  mask(study_area)

cov2 <- cov1 %>% 
  crop(study_area) %>% 
  mask(study_area)

stopifnot(compareGeom(psp2, study_area),
          compareGeom(cov2, study_area))

# extract and create dataframe --------------------------------------------

psp_df <- as.data.frame(psp2)
psp_df$cellnumber <- as.numeric(row.names(psp_df))

cov_df <- as.data.frame(cov2)
cov_df$cellnumber <- as.numeric(row.names(cov_df))

comb0 <- clim2 %>% 
  select(cellnumber, bio1, bio12) %>% 
  left_join(psp_df, by = c("cellnumber")) %>% 
  left_join(cov_df, by = c("cellnumber")) %>% 
  rename(MAT = bio1, MAP = bio12, psp = prcpPropSum)

comb1 <- comb0 %>% 
  # only keep cells for which we have all datasets (many rows don't have
  # psp data)
  drop_na()

# calculate 'quantiles'  --------------------------------------------------

# create dataframes of quantile summaries across climate variables,
# for the different cover variables (as per Holdrege et al. 2024, Fire Ecology)
# e.g. avg percent sagebrush cover across all pixels that are in the first 
# percentile of MAP
  
pred_vars <- c('MAT', 'MAP', 'psp')
  
response_vars <- c("sagebrush_p50_p95", 'PFG_p50_p95', 'AFG_p50_p50')

quant1 <- predvars2deciles(comb1, response_vars = response_vars,
                           pred_vars = pred_vars)

quant_filt1 <- predvars2deciles(comb1, response_vars = response_vars,
                                pred_vars = pred_vars,
                                filter_var = TRUE,
                                filter_vars = pred_vars)

cap1 <- 'Study area includes 1km cells with >50% of 30km cells w/ >1% sage'
cap_95 <- paste0(cap1, '\n RAP/RCMAP cover calculated as temporal median and spatial (707m) 95th percentile')
cap_50 <- paste0(cap1, '\n RAP cover calculated as temporal and spatial (707m) median')

pdf('figures/RAP/RAP_quantile_plots_v1.pdf',
    height = 7, width = 7)
decile_dotplot('sagebrush_p50_p95', quant1,
               ylab = 'Sagebrush cover')  +
  labs(caption = cap_95)
decile_dotplot('PFG_p50_p95', quant1,
               ylab = 'PFG cover') +
  labs(caption = cap_95)
decile_dotplot('AFG_p50_p50', quant1,
               ylab = 'AFG cover') +
  labs(caption = cap_50)

decile_dotplot_filtered('sagebrush_p50_p95', df = quant_filt1,
                        ylab = 'Sagebrush cover')

decile_dotplot_filtered('PFG_p50_p95', df = quant_filt1,
                        ylab = 'Perennials cover')

decile_dotplot_filtered('AFG_p50_p50', df = quant_filt1,
                        ylab = 'Annuals cover')
dev.off()

