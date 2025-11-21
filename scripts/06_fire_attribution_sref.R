# Purpose: Do calculations that help explain the 
# spatial pattern of wildfire under historical conditions
# Do this by calculating how different fire in a location
# due to a single variable at a time (attribution). Where
# the comparison is the variable being at it's median. 

# the 'sref' in the filename stands for 'spatial reference' 
# i.e. were' determine what's causing fire to be different from
# a median (across space) location, not a different simulation scenario

# dependencies ------------------------------------------------------------
library('tidyverse')
source('src/general_functions.R')
source('src/fire_functions.R')
source('src/qm_functions.R')
# params ------------------------------------------------------------------

source('src/params.R')
v_interp <- opt$v_interp
vr_name <- opt$vr_name
ref_graze <- opt$ref_graze
entire <- opt$entire
run <- opt$run
runv <- opt$runv
# load data ---------------------------------------------------------------

# create in "scripts/02b_fire_attribution.R"
d_l1 <- readRDS('data_processed/site_means/fire_dominant_drivers.RDS')
d_l1 <- filter_list_run(d_l1, run = run)
# created in "scripts/05_interpolation_weights.R"
wt <- read_csv(paste0('data_processed/interpolation_data/interpolation_weights_', 
                      v_interp, vr_name, '.csv'))

bio_wide0 <- d_l1$bio_wide

pred_vars <- d_l1$pred_vars

predict_fire_with <- d_l1$predict_fire_with

# calculating medians of predictor variable -------------------------------

bio_wide1 <- wt %>% 
  df_factor() %>% 
  filter(region == entire) %>% 
  right_join(bio_wide0, by = c('site')) %>% 
  filter(graze == ref_graze,
         RCP == 'Current')
  
stopifnot(nrow(bio_wide1) == 200)

med <- bio_wide1 %>% 
  summarise(
    across(all_of(pred_vars), .fns = \(x) Hmisc::wtd.quantile(x, 
                                                              weights = weight,
                                                              probs = 0.5))
  )


# predict fire ------------------------------------------------------------

# predict fire probability with one variable at a
# time adjusted to it's median
preds_1change <- map_dfr(pred_vars, function(var) {
  df <- bio_wide1
  df[[var]] <- med[[var]]
  df$fire_pred <- predict_fire_with(df)
  df$pred_var <- var
  df
})

# fire probability under reference conditions
bio_wide2 <- bio_wide1
bio_wide2$fire_pred_ref <- predict_fire_with(bio_wide1)

preds_1change2 <- preds_1change %>% 
  left_join(select(bio_wide2, -all_of(pred_vars)),
            join_by(site, region, weight, run, years, RCP, graze, id, GCM)) %>% 
  # change in fire due single variable
  mutate(delta_1var = fire_pred_ref - fire_pred)

# this is the fire probability if every variable were at it's median
predict_fire_with(med)
# 0.299838 


preds_1change3 <- preds_1change2 %>% 
  dplyr::select(site, years,RCP, graze, GCM, id, pred_var, delta_1var)


# save --------------------------------------------------------------------

filename <- paste0('data_processed/site_means/fire_attribution_spatial-reference_',
                   runv, '.csv')

write_csv(preds_1change3, filename)

