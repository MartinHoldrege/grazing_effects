# Purpose: See what soil depths and months water is being transpired by
# different plant functional types. Point is to see if there is a reason
# simulation results don't show a strong trend in sagebrush biomass with
# precipitation seasonality. 


# dependencies ------------------------------------------------------------

library(tidyverse)
library(DBI)
source('src/paths.R') 
source('src/general_functions.R')
source("src/fig_functions.R") # box_ann function defined here
theme_set(theme_custom1())

# read in data ------------------------------------------------------------

# KP compiled this database for current conditions for years 101-150 for moderate grazing, C4grass expansion off, and wildfire on.

path <- file.path(path_sw, "WildfireGrazingC4offFeb2024Runs/Output_Compiled_transp.sqlite")

connect <- dbConnect(RSQLite::SQLite(), path)

dbListTables(connect)

q <- "SELECT * FROM sw2_monthly_slyrs"
dat1 <- dbGetQuery(connect,
                   statement = q)

clim_db0 <- read_csv("data_processed/site_means/dbWeather_200sites.csv",
                    show_col_types = FALSE)

# make tidy ---------------------------------------------------------------

dat2 <- dat1 %>% 
  as_tibble() %>% 
  #head() %>% 
  pivot_longer(cols = matches('TRANSP_transp'),
               names_to = c('PFT', 'layer', 'summary'),
               names_pattern = "TRANSP_transp_(.+)_Lyr_(\\d)_(.+)") %>% 
  pivot_wider(names_from = 'summary') %>% 
  rename(transp = Mean,
         transp_sd = SD)

# discarding columns that have just 1 unique value
dat3 <- keep(dat2, \(x) !(is.character(x) & length(unique(x)) == 1))

# mean across years
dat4 <- dat3 %>% 
  group_by(site, Month, PFT, layer) %>% 
  summarize(transp = mean(transp)*10, # convert to mm
            .groups = 'drop') %>% 
  mutate(layer = as.numeric(layer),
         depth = lyr2depth(layer))

clim_db1 <- clim_db0 %>% 
  rename(site = Site_id,
         MAT = MAT_C,
         MAP = MAP_mm) %>% 
  select(-CorrTP, -Latitude, -Longitude, -Label)

# summarize ---------------------------------------------------------------

# calculate proportion of 

tot_t <- dat4 %>% 
  group_by(site) %>% 
  summarise(transp = sum(transp)
            ) %>% 
  left_join(clim_db1, by = 'site') %>% 
  # ecosystem water use efficiency
  mutate(eWUE = transp/MAP)

# sum across depths by PFT
t_by_pft0 <- dat4 %>% 
  group_by(site, PFT) %>% 
  summarize(transp = sum(transp)) %>% 
  left_join(tot_t, by = 'site',
            suffix = c('', "_tot")) %>% 
  mutate(transp_perc = transp/transp_tot*100 # percent of total transpiration (by PFT)
         )

t_by_depth0 <- dat4 %>% 
  group_by(site, depth, layer) %>% 
  summarize(transp = sum(transp))

dat5 <- dat4 %>% 
  left_join(t_by_depth0, by = c('site', 'depth', 'layer'),
            suffix = c('', '_tot_depth')) %>% 
  left_join(t_by_pft0, by = c('site', 'PFT'),
            suffix = c('', '_tot_pft')) %>% 
  mutate(transp_perc_depth = transp/transp_tot_depth*100,
         transp_perc_pft = transp/transp_tot_pft*100,
        transp2 = transp_per_cm(transp, layer))

t_by_depth <- left_join(t_by_depth0, clim_db1, by = 'site') %>% 
  # defining 'transp2' as transpiration as mm water per cm soil
  mutate(transp2 = transp_per_cm(transp, layer))

t_by_pft <- left_join(t_by_pft0, clim_db1, by = 'site')

# figures -----------------------------------------------------------------

# total transpiration

ggplot(t_by_depth, aes(CorrTP2, transp2)) +
  geom_point() +
  facet_wrap(~depth, scales = 'free_y') +
  geom_smooth(method = 'loess', se = FALSE) +
  expand_limits(y = range(t_by_depth$transp2[t_by_depth$depth != 5]))

ggplot(dat5, aes(CorrTP2, transp2)) +
  geom_point() +
  facet_grid(depth~PFT, scales = 'free_y') +
  geom_smooth(method = 'loess', se = FALSE)
  expand_limits(y = range(t_by_depth$transp2[t_by_depth$depth != 5]))