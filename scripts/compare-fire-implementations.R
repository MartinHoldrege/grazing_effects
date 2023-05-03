# Martin Holdrege

# code started 5/1/2023

# Purpose: Compare stepwat results of a few test sites done with the old
# implementation of stepwat, versus done with the new fire model (holdrege et al.)
# implented (for now this is the only change to the model for these test runs)

# dependencies ------------------------------------------------------------

library(DBI)

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

db_orig <- dbConnect(RSQLite::SQLite(), path_orig)

db_test <- map(test_dbs, function(x) {
  path <- file.path(db_loc, x)
  dbConnect(RSQLite::SQLite(), path)
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
bio1 <- map_dfr(db_test, dbGetQuery, statement = q1) %>% 
  # 2023 implementation where only the fire module has been updated
  # (flexible eind etc. not changed here)
  mutate(run = '23FireOnly')

bio_orig1 <- dbGetQuery(db_orig, q2)  %>% 
  # 2022 implimentation, with cheatgrass fire equation
  mutate(run = '22CheatgrassFire')

bio2 <- bind_rows(bio_orig1) %>% 
  as_tibble() %>% 
  mutate(RCP = ifelse(GCM == "Current", "Current", RCP),
         years = ifelse(GCM == "Current", "Current", years))


# summarize ---------------------------------------------------------------


# * fire ------------------------------------------------------------------
# summarize just the results on fire frequency

fire1 <- bio2 %>% 
  group_by(site, run, RCP, years, GCM) %>% 
  # average across years of a given simulation
  summarize(WildFire = mean(WildFire),
            .groups = 'drop_last') %>% 
  # median across GCMs
  summarize(WildFire = median(WildFire)) %>% 
  # % wildfire probability (calculating from # fires out of 200 years)
  mutate(prob = WildFire/200*100)# 
fire1
