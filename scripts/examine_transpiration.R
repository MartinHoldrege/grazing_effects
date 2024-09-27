# Purpose: See what soil depths and months water is being transpired by
# different plant functional types. Point is to see if there is a reason
# simulation results don't show a strong trend in sagebrush biomass with
# precipitation seasonality. 


# dependencies ------------------------------------------------------------

library(tidyverse)
library(DBI)
library(patchwork)
source('src/paths.R') 
source('src/general_functions.R')
source("src/fig_functions.R") # box_ann function defined here
source('src/fig_params.R')
theme_set(theme_custom1())


# params ------------------------------------------------------------------

suffix <- 'fire1eind1c4grass0co20grazeMCurrent' # for file names

# read in data ------------------------------------------------------------

# KP compiled this database for current conditions for years 101-150 for moderate grazing, C4grass expansion off, and wildfire on.

path <- file.path(path_sw, "WildfireGrazingC4offFeb2024Runs/Output_Compiled_transp.sqlite")

connect <- dbConnect(RSQLite::SQLite(), path)

dbListTables(connect)

q <- "SELECT * FROM sw2_monthly_slyrs"
dat1 <- dbGetQuery(connect,
                   statement = q)

# weater data compiled from stepwat weather database 00_query_weather_db.R
clim_db0 <- read_csv("data_processed/site_means/dbWeather_200sites.csv",
                    show_col_types = FALSE)
# mean monthly ppt and temp
clim_m_db0 <- read_csv("data_processed/site_means/dbWeather_200sites_monthly.csv",
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
         transp_sd = SD,
         month = Month)

# discarding columns that have just 1 unique value
dat3 <- keep(dat2, \(x) !(is.character(x) & length(unique(x)) == 1))

# mean across years
dat4 <- dat3 %>% 
  group_by(site, month, PFT, layer) %>% 
  summarize(transp = mean(transp)*10, # convert to mm
            .groups = 'drop') %>% 
  mutate(layer = as.numeric(layer),
         depth = lyr2depth(layer),
         depth_group = cut_depth(depth),
         month = factor(month, levels = 1:12,
                        labels = month.abb))

datg4 <- dat4 %>% # same as dat4 but with depth groups (object name has 'g')
  group_by(site, month, PFT, depth_group) %>% 
  summarize(transp = sum(transp));

clim_db1 <- clim_db0 %>% 
  rename(site = Site_id,
         MAT = MAT_C,
         MAP = MAP_mm) %>% 
  select(-CorrTP, -Latitude, -Longitude, -Label)
clim_vars <- c("corrTP2", 'MAP', 'MAT')

clim_m_db1 <- clim_m_db0 %>% 
  mutate(month = month.abb[month]) %>% 
  left_join(clim_db1, by = 'site') %>% 
  # percent of mean annual precip received during that month
  mutate(PPT_m_perc = PPT/MAP*100)

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
  summarize(transp = sum(transp)) 

t_by_depth0 <- dat4 %>% 
  group_by(site, depth, layer) %>% 
  summarize(transp = sum(transp))

t_by_depthg0 <- dat4  %>% 
  group_by(site, depth_group) %>% 
  summarize(transp = sum(transp), .groups = 'drop') 

# total by month
months2keep <- month.abb[3:8] # focus on a subset of months
t_by_m_depthg <- dat4  %>% 
  group_by(site, depth_group, month) %>% 
  filter(month %in% months2keep) %>% 
  summarize(transp = sum(transp), .groups = 'drop') %>% 
  mutate(transp2 = transp_per_cm(transp, depth_group)) %>% 
  left_join(clim_m_db1, by = c('site', 'month'))

datg5 <- datg4 %>% 
  group_by(site, month, depth_group) %>% 
  # percent of that depth's transpiration by that PFT
  mutate(transp_perc_depth = transp/sum(transp)*100) %>% 
  ungroup() %>% 
  group_by(site, month, PFT) %>% 
  # percent of total t for that pft that comes from the given
  # depth (during the given month)
  mutate(transp_perc_pft_m = transp/sum(transp)*100,
         transp2 = transp_per_cm(transp, depth_group)) %>% 
  ungroup() %>% 
  group_by(site, PFT) %>% 
  # proportion of a PFTs total annual transpiration, occuring in a given
  # month and depth
  mutate(transp_perc_pft = transp/sum(transp)*100) %>% 
  left_join(clim_m_db1, by = c('site', 'month'))

datg6 <- datg5 %>% 
  filter(month %in% months2keep)

t_by_depthg <- left_join(t_by_depthg0, clim_db1, by = 'site') %>% 
  # defining 'transp2' as transpiration as mm water per cm soil
  mutate(transp2 = transp_per_cm(transp, depth_group))

t_by_pft <- t_by_pft0 %>% 
  left_join(tot_t, by = 'site',
            suffix = c('', "_tot")) %>% 
  mutate(transp_perc = transp/transp_tot*100 # percent of total transpiration (by PFT)
  )

t_by_pft_depthg <- datg5 %>% 
  group_by(site, PFT, depth_group) %>% 
  summarize(transp = sum(transp)) %>% 
  left_join(t_by_depthg, by = c('site', 'depth_group'),
            suffix = c('', '_tot_depth')) %>% 
  # percent of transpiration at a given depth that comes from a particular
  # PFT
  mutate(transp_perc_depth = transp/transp_tot_depth*100)


# save key dataframes -----------------------------------------------------
# for use in the stepwat_params/scripts/02_evaluate-test-sites.R script

tmp <- list(
  tot_t = tot_t,
  t_by_depthg = t_by_depthg,
  t_by_pft = t_by_pft,
  t_by_pft_depthg,
  datg6 = datg6
)

saveRDS(tmp, paste0('data_processed/site_means/transpiration_dfs_', suffix))

# figures -----------------------------------------------------------------


# * totals (across PFT and depth) -------------------------------------------

# base (b) of the next figures
b0 <- function() {
  list(geom_point(),
       geom_smooth(method = 'loess', se = FALSE))
}

b1 <- function() {
  list(b0(),
       labs(y = lab_transp0))
}

pdf(paste0('figures/water/transpiration_total_', suffix, '.pdf'),
    width = 8, height = 8)
# climate variables

GGally::ggpairs(clim_db1[, c("MAT", "MAP", "CorrTP2")])

# total transpiration


g1 <- ggplot(tot_t, aes(CorrTP2, transp)) +
  b1() +
  labs(x = lab_corrtp)

g2 <- ggplot(tot_t, aes(MAT, transp)) +
  b1()

g3 <- ggplot(tot_t, aes(MAP, transp)) +
  b1()

(g1 + g2)/(g3 + patchwork::plot_spacer())

# ecosystem water use efficiency
g1 <- ggplot(tot_t, aes(CorrTP2, eWUE)) +
  b0() +
  labs(x = lab_corrtp)

g2 <- ggplot(tot_t, aes(MAT, eWUE)) +
  b0()

g3 <- ggplot(tot_t, aes(MAP, eWUE)) +
  b0()

(g1 + g2)/(g3 + patchwork::plot_spacer())



# * totals by depth -------------------------------------------------------

# total transpiration by depth
ggplot(t_by_depthg, aes(CorrTP2, transp)) +
  b0() +
  facet_wrap(~depth_group, ncol = 2) +
  labs(x = lab_corrtp,
       y = lab_transp0)

# total by depth and month
b2 <- \() list(b0(), 
               facet_grid(month~depth_group), 
               labs(y = lab_transp0))

ggplot(t_by_m_depthg, aes(CorrTP2, transp)) +
  b2() +
  labs(x = lab_corrtp)

ggplot(t_by_m_depthg, aes(MAP, transp)) +
  b2()

ggplot(t_by_m_depthg, aes(MAT, transp)) +
  b2()

dev.off()
# * PFT level (across depth and month) ------------------------------------------

pdf(paste0('figures/water/transpiration_by-pft_', suffix, '.pdf'),
    width = 9, height = 8)

ggplot(t_by_pft, aes(CorrTP2, transp_perc)) +
  b0() +
  facet_wrap(~PFT, ncol = 2) +
  labs(x = lab_corrtp,
       y = "% of total transpiration")

ggplot(t_by_pft, aes(CorrTP2, transp)) +
  b0() +
  facet_wrap(~PFT, ncol = 2, scales = "free_y") +
  labs(x = lab_corrtp,
       y = lab_transp0)

ggplot(t_by_pft_depthg, aes(CorrTP2, transp_perc_depth)) +
  b0() +
  facet_grid(depth_group~PFT)+
  labs(x = lab_corrtp,
       y = "% of total transpiration of that depth",
       subtitle = 'rows sum to 100%')

ggplot(t_by_pft_depthg, aes(CorrTP2, transp)) +
  b0() +
  facet_grid(depth_group~PFT, scales = 'free_y')+
  labs(x = lab_corrtp,
       y = lab_transp0)

for (m in months2keep) {
  g <- datg6 %>% 
    filter(month == m) %>% 
    ggplot(aes(CorrTP2, transp_perc_depth) # percent of total Transp for that pft
           # coming from that month and depth
           ) +
    geom_point(aes(color = PPT_m_perc, size = transp_perc_pft))+
    geom_smooth(method = 'loess', se = FALSE)+
    scale_color_gradient2(low = "#a50026", midpoint = 1/12*100, 
                        mid = "lightgrey", high = "#313695",
                        name = '% MAP this month',
                        limits = range(datg6$PPT_m_perc, na.rm = TRUE))+
    scale_size(range = c(0.2, 3),
               name = '% of PFTs total T',
               limits = range(datg6$transp_perc_pft, na.rm = TRUE)) +
    facet_grid(depth_group~PFT) +
    labs(x = lab_corrtp,
       y = "% of total transpiration of that depth & month",
       title = m,
       subtitle = 'rows sum to 100%') +
    expand_limits(y = c(0, 100)) +
    theme(legend.title = element_text(size = rel(0.5)))
  print(g)
}


for (m in months2keep) {
  g <- datg6 %>% 
    filter(month == m) %>% 
    ggplot(aes(CorrTP2, transp_perc_pft_m)) +
    b0() +
    facet_grid(depth_group~PFT)+
    labs(x = lab_corrtp,
         y = "% of PFT total transpiration for the month",
         title = m,
         subtitle = 'columns sum to 100%') +
    expand_limits(y = c(0, 100))
  print(g)
}

for (m in months2keep) {
  g <- datg6 %>% 
    filter(month == m) %>% 
    ggplot(aes(PPT, transp)) +
    geom_point(aes(color = PPT_m_perc)) +
    geom_smooth(method = 'loess', se = FALSE) +
    scale_color_gradient2(low = "#a50026", midpoint = 1/12*100, 
                        mid = "lightgrey", high = "#313695",
                        name = '% MAP this month',
                        limits = range(datg6$PPT_m_perc, na.rm = TRUE)) +
    facet_grid(depth_group~PFT) +
    labs(x= paste(m, 'precipitation (mm)'),
         y = lab_transp0,
         title = m) +
    expand_limits(y = c(0, 30))
  print(g)
}

dev.off()

