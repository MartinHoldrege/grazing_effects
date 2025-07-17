# Martin Holdrege

# code started 5/1/2023

# Purpose: Compare stepwat results of a few test sites done with the old
# implementation of stepwat, versus done with the new fire model (holdrege et al.)
# implented (for now this is the only change to the model for these test runs)

# dependencies ------------------------------------------------------------

library(tidyverse)
library(DBI)
library(sf)
source("src/general_functions.R")
source("src/fig_params.R")
source("../cheatgrass_fire/src/basemaps.R") # for basemaps_g
theme_set(theme_classic())


# params ------------------------------------------------------------------

path <- "D:/USGS/large_files/stepwat/"


# site locations ----------------------------------------------------------

sites1 <- read_csv(file.path(path, "WildfireMay2023TestRuns/STEPWAT2.testing.sites.csv"))

# connect to databases ----------------------------------------------------

# path to the databases from the original data runs (i.e. only cheatgrass
# fire implementation)
path_orig <- "data_raw/Output.biomass.200sites.grazing.cheatgrass.Nov2021.sqlite"

sites <- c(1, 4, 21,
           25, 84, 87, 89, 178) # sites for which test runs were done (only some done for 178)

# new wildfire module implemented
path_test_runs1 <- list.files(file.path(path, "WildfireMay2023TestRuns"),
                             pattern = ".sqlite",
                             full.names = TRUE)

# run where not wildfire was implemented
path_test_runs2 <- list.files(file.path(path, "WildfireMay2023TestRuns/NoWildfireDatabases"),
                              pattern = ".sqlite",
                              full.names = TRUE)

path_test_runs3 <- file.path(path, "WildfireEindJuly2023TestRuns/Output_Compiled.sqlite")

path_test_runs4 <- file.path(path, "WildfireNov2023TestRuns/Output_Compiled_Nov2023test.sqlite")

db_orig <- dbConnect(RSQLite::SQLite(), path_orig)

db_test1 <- map(path_test_runs1, function(x) {
  dbConnect(RSQLite::SQLite(), x)
})

db_test2 <- map(path_test_runs2, function(x) {
  dbConnect(RSQLite::SQLite(), x)
})

db_test3 <- dbConnect(RSQLite::SQLite(), path_test_runs3)
db_test4 <- dbConnect(RSQLite::SQLite(), path_test_runs4)

# db queries --------------------------------------------------------------

# table is small enough to fit in memory, so just loading the whole table
# then summarizing with R below
q1 <- paste("SELECT *",
            "FROM BIOMASS",
            "WHERE Year > 100;")

q2 <- paste("SELECT *",
            "FROM BIOMASS",
            "WHERE Year > 100 AND (site == ", 
            paste(sites, collapse = ' OR site == '), ");")

# separately querying the two tables
bio1 <- map_dfr(db_test1, dbGetQuery, statement = q1) %>% 
  # run where no fire probability was measured
  mutate(run = '2023Fire1')

bio1a <- map_dfr(db_test2, dbGetQuery, statement = q1) %>% 
  # fire 2023 implementation where only the fire module has been updated
  # (flexible eind etc. not changed here)
  mutate(run = '2023NoFire')

bio1b <- dbGetQuery(db_test3, q1) %>% 
  mutate(run = '2023Fire1Eind')

bio1c <- dbGetQuery(db_test4, q1) %>% 
  mutate(run = '2023Fire2Eind')

bio_orig1 <- dbGetQuery(db_orig, q2)  %>% 
  # 2022 implimentation, with cheatgrass fire equation
  mutate(run = '2022CheatgrassFire')

bio2 <- bind_rows(bio1, bio1a, bio1b, bio1c, bio_orig1) %>% 
  as_tibble() %>% 
  mutate(RCP = rcp2factor(RCP),
         years = years2factor(years),
         id = paste(RCP, years, sep = "_"),
         # sequential ordering of site numbers
         site_simple = as.factor(as.numeric(factor(site))),
         run = factor(run),
         run = fct_relevel(run, '2023NoFire')) 

# summarize ---------------------------------------------------------------

# * fire ------------------------------------------------------------------
# summarize just the results on fire frequency

fire1 <- bio2 %>% 
  group_by(id, site, site_simple, run, RCP, years, GCM) %>% 
  # average across years of a given simulation
  summarize(WildFire = mean(WildFire),
            .groups = 'drop_last') %>% 
  # median across GCMs
  summarize(WildFire = median(WildFire)) %>% 
  # % wildfire probability (calculating from # fires out of 200 years)
  mutate(prob = WildFire/200*100)# 


# *biomass ----------------------------------------------------------------
PFTs <- c("sagebrush", "a.cool.forb", "a.warm.forb", "p.cool.forb", "p.warm.forb",
          "a.cool.grass", "p.cool.grass", "p.warm.grass")

# median biomass across GCMs
bio_med1 <- bio2 %>% 
  select(run, site, site_simple, GCM, years, RCP, Year, all_of(PFTs), id) %>% 
  pivot_longer(cols = all_of(PFTs),
               values_to = "biomass",
               names_to = "PFT") %>% 
  group_by(id, site, site_simple, run, RCP, years, PFT, GCM) %>% 
  summarise(biomass = mean(biomass), .groups = "drop_last") %>% 
  summarize(biomass = median(biomass))



# figures -----------------------------------------------------------------


shapes <- scale_shape_manual(values = seq(0, length(sites)))
base <- function() {
  list(theme(legend.position = 'top'),
       shapes)
}

cap1 <- paste("Simulations run for", length(sites), "sites"
              ,"\nFire2 refers to the newest fire equation")
pdf("figures/fire/compare-fire-implementations_v3.pdf",
    width = 10, height = 8)

# *map of site locations ------------------------------------------------------

# kyle added a site to the great basin, adding coords here
sites2 <- tibble(X_WGS84 = -117.82083, Y_WGS84 = 39.3125,
                 `Site(new)` = 178) %>% 
  bind_rows(sites1) %>% 
  mutate(site_simple = as.numeric(as.factor(`Site(new)`))) %>% 
  sf::st_as_sf(coords = c("X_WGS84", "Y_WGS84"),
               crs = sf::st_crs("EPSG:4326")) %>% 
  sf::st_transform(crs = crs)


ggplot() +
  geom_sf_text(data = sites2, aes(label = site_simple)) +
  #geom_text()
  basemap_g(bbox) +
  labs(title = 'locations of test sites')

# *fire -------------------------------------------------------------------

ggplot(fire1, aes(id, prob, color = run, shape = run, fill = run, group = run)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Scenario",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1)

ggplot(fire1, aes(id, prob, color = run, shape = run, fill = run, group = run)) +
  geom_text(aes(label = site_simple), position = position_dodge(width = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Scenario",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1)


ggplot(fire1, aes(run, prob, group = site_simple, color = site_simple, shape = site_simple)) +
  geom_point() +
  geom_line() +
  facet_wrap(~id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Model run",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1) +
  base()


ggplot(fire1, aes(id, prob, color = site_simple, group = site_simple)) + 
  geom_line()+
  geom_point(aes(shape = site_simple)) +
  facet_wrap(~run) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Scenario",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1) +
  base()

# * biomass ---------------------------------------------------------------
g <- ggplot(bio_med1, aes(id, biomass, color = run, shape = run, fill = run, group = run)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~PFT, scales = 'free_y') +
  labs(x = "Scenario",
       y = lab_bio0,
       caption = cap1)

g + geom_point(position = position_dodge(width = 0.5)) 
g + geom_text(aes(label = site_simple), position = position_dodge(width = 0.5)) 


for(x in PFTs) {
  g <- ggplot(bio_med1[bio_med1$PFT == x, ], aes(run, biomass, group = site,
                                                 color = site_simple,
                                                 shape =site_simple)) +
    geom_point() +
    geom_line() +
    facet_wrap(~id) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    labs(x = "Model run",
         y = lab_bio0,
         subtitle = x,
         caption = cap1) +
    base()
  print(g)
}




dev.off()

