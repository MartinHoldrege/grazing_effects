# Martin Holdrege

# code started 5/1/2023

# Purpose: Compare stepwat results of a few test sites done with the old
# implementation of stepwat, versus done with the new fire model (holdrege et al.)
# implented (for now this is the only change to the model for these test runs)

# dependencies ------------------------------------------------------------

library(tidyverse)
library(DBI)
source("src/general_functions.R")
theme_set(theme_classic())

# params ------------------------------------------------------------------

path <- "D:/USGS/large_files/stepwat/"

# connect to databases ----------------------------------------------------

# path to the databases from the original data runs (i.e. only cheatgrass
# fire implementation)
path_orig <- "data_raw/Output.biomass.200sites.grazing.cheatgrass.Nov2021.sqlite"

sites <- c(1, 4, 21,
           25, 84, 87, 89) # sites for which test runs were done
# where test dbs are located
db_loc <- "D:/USGS/large_files/stepwat/WildfireMay2023TestRuns"
test_dbs <- paste0("Output_site_", sites, ".sqlite")
names(test_dbs) <- paste0("site_", sites)
path_test_runs <- file.path(db_loc, paste())

# new wildfire module implemented
path_test_runs1 <- list.files(file.path(path, "WildfireMay2023TestRuns"),
                             pattern = ".sqlite",
                             full.names = TRUE)

# run where not wildfire was implemented
path_test_runs2 <- list.files(file.path(path, "WildfireMay2023TestRuns/NoWildfireDatabases"),
                              pattern = ".sqlite",
                              full.names = TRUE)

db_orig <- dbConnect(RSQLite::SQLite(), path_orig)

db_test1 <- map(path_test_runs1, function(x) {
  dbConnect(RSQLite::SQLite(), x)
})

db_test2 <- map(path_test_runs2, function(x) {
  dbConnect(RSQLite::SQLite(), x)
})

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
  mutate(run = '2023Fire')

bio1a <- map_dfr(db_test2, dbGetQuery, statement = q1) %>% 
  # 2023 implementation where only the fire module has been updated
  # (flexible eind etc. not changed here)
  mutate(run = '2023NoFire')

bio_orig1 <- dbGetQuery(db_orig, q2)  %>% 
  # 2022 implimentation, with cheatgrass fire equation
  mutate(run = '2022CheatgrassFire')

bio2 <- bind_rows(bio1, bio1a, bio_orig1) %>% 
  as_tibble() %>% 
  mutate(RCP = rcp2factor(RCP),
         years = years2factor(years),
         id = paste(RCP, years, sep = "_"))


# summarize ---------------------------------------------------------------

# * fire ------------------------------------------------------------------
# summarize just the results on fire frequency

fire1 <- bio2 %>% 
  group_by(id, site, run, RCP, years, GCM) %>% 
  # average across years of a given simulation
  summarize(WildFire = mean(WildFire),
            .groups = 'drop_last') %>% 
  # median across GCMs
  summarize(WildFire = median(WildFire)) %>% 
  # % wildfire probability (calculating from # fires out of 200 years)
  mutate(prob = WildFire/200*100)# 


# figures -----------------------------------------------------------------

# *fire -------------------------------------------------------------------
fire1

cap1 <- paste("Simulations run for", length(sites), "sites")
pdf("figures/fire/compare-fire-implementations_v1.pdf",
    width = 6, height = 5)
ggplot(fire1, aes(id, prob, color = run, shape = run, fill = run)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Scenario",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1)


ggplot(fire1, aes(run, prob, group = site)) +
  geom_point() +
  geom_line() +
  facet_wrap(~id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Scenario",
       y = "Observed fire probability (%) in STEPWAT2",
       caption = cap1)
dev.off()

