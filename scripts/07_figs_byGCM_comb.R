

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
# load files --------------------------------------------------------------

# file created in 05_interpolated_summarize_sei_scd-adj.R
# mean SEI and % core by ecoregion and GCM
p <- paste0('data_processed/raster_means/', 
            runv, yr_lab, vr_name,
            '_sei-mean_pcent-csa_scd-adj_by-GCM-region.csv')

if(file.exists(p)) {
  sei_byGCM1 <- read_csv(p)
} else {
  # delete this later
  message('loading test data for code development')
  sei_byGCM1 <- readRDS('test.RDS')
}

# average of climatic drivers, per gcm and region
# file output by 05_summarize_fire_drivers.R
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                            '_fire-driver-means_by-ecoregion.csv'))

# expected burned area by gcm
# created in 06_fire_area.R
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_", 
                           v_input, "_", run, ".csv"))



# combining dataframes ----------------------------------------------------

ba_gcm2 <- ba_gcm1 %>% 
  mutate(variable = 'ba') %>% 
  # burned area isn't really a median (it is more equivelant to 
  # weighted average of fire probability), but so that
  # this df has the same structure
  rename(p50 = area, region = ecoregion) %>% 
  select(-run, -run2, -type)


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

vars_clim <- c('MAT', 'MAP', 'PSP')
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

df <- comb1 %>% 
  # treating grazing of climate variables as moderate 
  # (climate doesn't differ with grazing )
  mutate(graze = ifelse(is.na(graze), ref_graze, as.character(graze)),
         graze = graze2factor(graze)) %>% 
  filter(region == levels(region)[1],
         graze == ref_graze)

make_wide_4crossplot <- function(df, var1, var2) {
  df_wide <- df %>% 
    filter(.data$variable %in% c(var1, var2)) %>% 
    rename(median = p50) %>% 
    pivot_wider(values_from = c('p25', 'median', 'p75'), names_from = variable,
                names_glue = "{variable}_{.value}")
  df_wide
  
}

# CONTINUE HERE