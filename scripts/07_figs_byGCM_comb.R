# Purpose:
# create figure showing how responses of biomass, SEI, and fire
# differ between GCMs

# dependencies ------------------------------------------------------------
library(tidyverse)
library(patchwork)
source('src/fig_params.R')
source('src/fig_functions.R')
source('src/general_functions.R')
# params ------------------------------------------------------------------
source('src/params.R')
runv <- opt$runv
yr_lab <- opt$yr_lab
vr_name <- opt$vr_name
ref_graze <- opt$ref_graze
vr <- opt$vr
options(readr.show_col_types = FALSE)

theme_set(theme_custom1())
# load files --------------------------------------------------------------

# file created in 05_interpolated_summarize_sei_scd-adj.R
# mean SEI and % core by ecoregion and GCM
p <- paste0('data_processed/raster_means/', 
            runv, yr_lab, vr_name,
            '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv')

sei_byGCM1 <- read_csv(p)

# average of climatic drivers, per gcm and region
# file output by 05_summarize_fire_drivers.R
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                            '_fire-driver-means_by-ecoregion.csv'))

area1 <-  read_csv(paste0('data_processed/area/ecoregion-area_v2', vr_name,
                          '.csv'))

# expected burned area by gcm
# created in 06_fire_area.R
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_v2",
                           "_", run, ".csv"))

# combining dataframes ----------------------------------------------------

ba_gcm2 <- ba_gcm1 %>% 
  mutate(variable = 'ba') %>% 
  # burned area isn't really a median (it is more equivelant to 
  # weighted average of fire probability), but so that
  # this df has the same structure
  rename(p50 = area, region = ecoregion) %>% 
  select(-run, -run2, -type) %>% 
  left_join(area1, by = c('region' = 'ecoregion')) %>% 
  mutate(p50 = p50/area*100) %>%  # converting to %/year (from ha/year)
  select(-area)

sei_byGCM2 <- sei_byGCM1 %>% 
  select(region, RCP, years, graze, GCM, matches('SEI_p')) %>% 
  mutate(variable = 'SEI') %>% 
  rename_with(.fn = \(x) str_replace(x, 'SEI_', ''),
              .cols = matches('SEI_p\\d+'))

drivers2 <- drivers1 %>% 
  select(-weight, -mean) %>% 
  rename(p50 = median) 

comb1 <- bind_rows(drivers2, ba_gcm2, sei_byGCM2) %>% 
  mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
  df_factor()


# crossplot figure --------------------------------------------------------


# * variable setup ------------------------------------------------------------

vars_clim <- c('MAT', 'MAP', 'psp')
vars_veg <- c('Sagebrush', 'Pherb', 'Aherb')
vars_pred <- c(vars_clim, rev(vars_veg))

# which x and y variables are plotted against each other
vars_l <- list(
  veg = expand_grid(y = vars_veg, x = vars_clim),
  SEI = expand_grid(y = "SEI", x = vars_pred),
  # sagebrush does not affect fire
  ba = expand_grid(y = 'ba', x = vars_pred[vars_pred != 'Sagebrush'])
)
  
vars_df <- bind_rows(vars_l)


# * df prep ---------------------------------------------------------------

yr <- '2070-2100'
rcp <- 'RCP45'
region <- levels(comb1$region)[1]
df <- comb1 %>% 
  # treating grazing of climate variables as moderate 
  # (climate doesn't differ with grazing )
  mutate(graze = ifelse(is.na(graze), ref_graze, as.character(graze)),
         graze = graze2factor(graze),
         GCM = gcm2factor(GCM, include_current = TRUE)) %>% 
  filter(region == !!region,
         graze == ref_graze,
         years %in% c('Current', yr),
         RCP %in% c('Current', rcp))

make_wide_4crossplot <- function(df, var1, var2) {
  df_wide <- df %>% 
    filter(.data$variable %in% c(var1, var2)) %>% 
    rename(median = p50) %>% 
    pivot_wider(values_from = c('p25', 'median', 'p75'), names_from = variable,
                names_glue = "{variable}_{.value}")
  df_wide
}

x <- 'MAT'
y <- 'ba'


plots <- pmap(vars_df, function(x, y) {
  df_wide <- make_wide_4crossplot(df, var1 = x, var2 = y)
  crossplot_1panel(df_wide, x, y, colors = cols_GCM2,
                   shapes = shapes_GCM2,
                   labeller = driver_labeller(delta = FALSE))
})


design <- '
abcvvv
defvvv
ghivvv
jklmno
pqrst#
'
plots2 <- c(plots, list(guide_area()))

g <- patchwork::wrap_plots(plots2,
                      design = design) + 
  plot_layout(guides = 'collect', axis_titles = 'collect',
              axes = 'collect')

filename <- paste0('figures/by_gcm_comb/response_v_driver_20panel_',
                   vr, "_", words2abbrev(region), '_', rcp, '_', yr, '_',
                   runv, '.png')
ggsave(
  filename = filename,
  plot = g,
  width = 11,
  height = 10,
  dpi = 600
)
