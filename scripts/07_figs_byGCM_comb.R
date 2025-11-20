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
years <- opt$years
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

# change in drivers, per gcm and region
# file output by 06_summarize_fire_drivers.R
drivers1 <- read_csv(paste0('data_processed/raster_means/', run, vr_name,
                            '_delta-fire-driver_by-ecoregion.csv'))

area1 <-  read_csv(paste0('data_processed/area/ecoregion-area_v2', vr_name,
                          '.csv'))

# expected burned area by gcm
# created in 06_fire_area.R
ba_gcm1 <- read_csv(paste0("data_processed/area/expected-burn-area_by-GCM_v2",
                           "_", run, ".csv"))


# calculate change --------------------------------------------------------
# changes relative to historical conditions

sei_byGCM2 <- sei_byGCM1 %>% 
  # combined goa & CSA area
  #mutate(class_area = percent_goa/100*area + percent_csa/100*area) %>% 
  # just CSA
  mutate(class_area = percent_csa/100*area) %>% 
  select(-id, -matches('SEI'), -area, -group, -type, -matches('percent_'),
         -run)

sei_cref_gcm <- join4cref(sei_byGCM2, by = c('region', 'graze')) %>% 
  # % change in area
  mutate(median = (class_area_fut - class_area_cur)/class_area_cur *100,
         variable = 'csa') %>% 
  select(-matches('class_area'))

# * burned area -----------------------------------------------------------

ba_gcm2 <- ba_gcm1 %>% 
  rename(ba = area, region = ecoregion) %>% 
  select(region, RCP, years, GCM, graze, ba)

# change in burned area (%) relative to historical climate
ba_cref_gcm <- join4cref(ba_gcm2, by = c('region', 'graze')) %>% 
  # just to make naming consistent (this is actually a median)
  mutate(median = (ba_fut - ba_cur)/ba_cur*100,
         variable = 'ba_delta_perc') %>% 
  select(-ba_cur, -ba_fut)


# combining dataframes ----------------------------------------------------


drivers2 <- drivers1 %>% 
  select(-weight, -mean) 

comb1 <- bind_rows(drivers2, ba_cref_gcm, sei_cref_gcm) %>% 
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
  SEI = expand_grid(y = "csa", x = vars_pred),
  ba = expand_grid(y = 'ba_delta_perc', x = vars_pred)
)
  
vars_df <- bind_rows(vars_l)
vars_df$tag <- fig_letters[1:nrow(vars_df)]
# * df prep ---------------------------------------------------------------

comb2 <- comb1 %>% 
  # treating grazing of climate variables as moderate 
  # (climate doesn't differ with grazing )
  mutate(graze = ifelse(is.na(graze), ref_graze, as.character(graze)),
         graze = graze2factor(graze),
         GCM = gcm2factor(GCM, include_current = TRUE))


# * create fig ------------------------------------------------------------


yr <- years

region <- levels(comb1$region)[1]

rcps <- unique(comb1$RCP)
for(rcp in rcps) {
  df <-  comb2%>% 
    filter(region == !!region,
           graze == ref_graze,
           years == yr,
           RCP == rcp)
  
  legend_title <- paste0('GCM\n', rcp_label(rcp, yr, include_parenth = TRUE))
  
  g <- crossplot_multipanel(df = df,
                            vars_df = vars_df,
                            legend_title = legend_title)

  filename <- paste0('figures/by_gcm_comb/response_v_driver_', nrow(vars_df), 
                   'panel_cref_csa_',
                   vr, "_", words2abbrev(region), '_', rcp, '_', yr, '_',
                   runv, '.png')
  ggsave(
    filename = filename,
    plot = g,
    width = 11,
    height = 10,
    dpi = 600
  )
}

