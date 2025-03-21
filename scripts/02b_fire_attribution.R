# Purpose: Determine which predictor variables of the fire equation
# are most responsible the change in fire probability at a given site
# under a given climate scenario

# Author: Martin Holdrege

# Started: March 19, 2025


# params ------------------------------------------------------------------


# dependencies ------------------------------------------------------------

library(tidyverse)
source("src/fire_functions.R")
source("src/qm_functions.R")
source("src/general_functions.R")
source("src/fig_functions.R")

# read in data ------------------------------------------------------------

# list of dataframes created in "scripts/02_summarize_bio.R" 
l <- readRDS('data_processed/site_means/summarize_bio.RDS') 

clim_all2 <- l$clim_all2

bio1 <- l$pft5_bio1

# quantiles used in STEPWAT2 for quantile mapping in fire equation
# created in scripts/bio_matching/02_quantile_matching_simple.R
# NOTE--as the quantiles used are updated in STEPWAT, update the file that is read in
quants <- read_csv("data_processed/qm/quantiles_for_qm_0.5match_fire0_eind1_c4grass1_co20.csv")


# create qm functions -----------------------------------------------------
# functions for quantile mapping

qm_Aherb <- qm_quant_factory(
  from = quants$biomass_stepwat[quants$PFT == 'afgAGB'], 
  to = quants$biomass_rap[quants$PFT == 'afgAGB'])

qm_Pherb <- qm_quant_factory(
  from = quants$biomass_stepwat[quants$PFT == 'pfgAGB'], 
  to = quants$biomass_rap[quants$PFT == 'pfgAGB'])


# prepare dataframes  --------------------------------------------------------

bio2 <- bio1 %>% 
  filter(str_detect(run, 'fire1')) %>% 
  calc_aherb() %>% 
  filter(PFT %in% c('Aherb', 'Pherb')) %>% 
  select(run, years, RCP, graze, id, site, GCM, biomass, PFT) %>% 
  left_join(clim_all2, by = join_by(years, RCP, site, GCM))

bio3a <- bio2 %>% 
  pivot_wider(values_from = biomass, names_from = PFT) %>% 
  df_factor()

# calculate the quantile mapped biomass,
# so that predictions made with fire equation are appropriate
bio3 <- bio3a %>% 
  mutate(Aherb = qm_Aherb(Aherb),
         Pherb = qm_Pherb(Pherb))

bio_cur1 <- bio3 %>% 
  filter(years == 'Current')

bio_fut1 <- bio3 %>% 
  filter(years != 'Current')

pred_vars <- c("MAP", "MAT", "psp", "Aherb", "Pherb")

bio_cur2 <- bio_cur1 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_cur',
               values_to = "value_cur")

# columns labeled with the pred vars, are the 'current' values
bio_cur3 <- bio_cur2 %>% 
  select(-years, -RCP, -GCM, -id) %>% 
  full_join(bio_fut1, by = c('run', 'graze', 'site'),
             relationship = "many-to-many")

test <- bio_cur3 %>% 
  filter(site == 1,
         years == levels(years)[3],
         RCP == levels(RCP)[2],
         GCM == unique(GCM)[1],
         run == unique(run)[1],
         graze == "Light")

stopifnot(nrow(test) == length(pred_vars))

# all pred vars (columns) are equal to future values, 
# except the one named in pred_var_cur column, where the current (historical)
# value is instead used
one_change <- bio_cur3 %>% 
  pivot_longer(cols = all_of(pred_vars),
               names_to = 'pred_var_fut',
               values_to = 'value_fut') %>% 
  # replacing one pred var at a time, with the current value
  mutate(value = ifelse(pred_var_fut == pred_var_cur, value_cur, value_fut)) %>% 
  select(-value_cur, -value_fut) %>% 
  pivot_wider(values_from = value, names_from = pred_var_fut) 


# predict fire probability ------------------------------------------------

predict_fire_with <- function(df) {
  p <- with(df, predict_fire(mat = MAT, map = MAP, psp = psp, afg = Aherb,
                        pfg = Pherb,run_checks = TRUE))
  p*100 # convert to %
}

# predicted fire probability when 1 variable changed to the current value
one_change$fire_prob_1var <- predict_fire_with(one_change)

# predicted fire probability under future conditions
bio_fut2 <- bio_fut1
bio_fut2$fire_prob_fut <- predict_fire_with(bio_fut2)

# predicted fire probability under current conditions
bio_cur1b <- bio_cur1
bio_cur1b$fire_prob_cur <- predict_fire_with(bio_cur1b)

# combine
one_change2 <- bio_cur1b %>% 
  select(run, site, graze, fire_prob_cur) %>% 
  right_join(one_change) %>% 
  select(-all_of(pred_vars)) %>% 
  left_join(select(bio_fut2, run, years, RCP, GCM, graze, id, site, fire_prob_fut)) %>% 
  # How much different is fire probability, when 1 variable isn't changed
  # (i.e. 1 variable changed to current value)
  mutate(delta_1var = fire_prob_fut - fire_prob_1var ,
         rcp_years = rcp_label(RCP, years, include_parenth = FALSE)) 


# summarize results -------------------------------------------------------

dominant_driver0 <- one_change2 %>% 
  group_by(run, site, graze, RCP, years, id, GCM) %>% 
  filter(abs(delta_1var) == max(abs(delta_1var))) %>% 
  mutate(n = n()) %>% 
  rename(dominant_driver = pred_var_cur)


# if n>1 that means there were ties for the dominant driver
# and these would need to be dealt with
stopifnot(all(dominant_driver0$n == 1))

# summary across GCMs
set.seed(1234)
dominant_driver1 <- dominant_driver0 %>% 
  group_by(run, site, graze, RCP, years, id) %>% 
  summarize(dominant_driver = mode_rand(dominant_driver),
            .groups = 'drop') %>% 
  mutate(dominant_driver = factor(dominant_driver,
                                  levels = !!pred_vars))

one_change3 <- one_change2 %>% 
  group_by(run, site, graze, years, RCP, id, pred_var_cur) %>% 
  summarize(across(.cols = c('fire_prob_fut', "fire_prob_cur", "delta_1var"),
                   .fns = median),
            delta_1var_abs = median(abs(delta_1var)),
            .groups = 'drop') %>% 
  mutate(pred_var_cur = factor(pred_var_cur, levels = !!pred_vars))


# checks ------------------------------------------------------------------


test <- dominant_driver1 %>% 
  group_by(run, id) %>% 
  summarize(n = n()) %>% 
  pull(n)
stopifnot(all(test == 200))  

# output ------------------------------------------------------------------

out <- list(
  one_change3 = one_change3,
  dominant_driver1 = dominant_driver1,
  pred_vars = pred_vars
)

saveRDS(out, 'data_processed/site_means/fire_dominant_drivers.RDS')
