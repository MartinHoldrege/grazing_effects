# Martin Holdrege

# Script started 1/4/2022

# Purpose: to access the database of STEPWAT2 biomass output, and
# create summary df's for later use. 
# Simulations were done with 4 grazing levels, 200 sites, and 52 GCM*RCP*time 
# period combinations, with the cheatgrass/fire cycle turned on. 
# simulations run for 150 years, last 50 yrs used.


# dependencies ------------------------------------------------------------

library(tidyverse)
library(DBI)

# connect to db -----------------------------------------------------------

path_sw <- "D:/USGS/large_files/stepwat/" # db files stored on hard-drive

# naming convention:
# Wildfire – 0/1
# Grazing intensity –  L, M, H, VH
# Site-specific eind – 0/1
# C4 grass expansion under future conditions – 0/1
# CO2 – 0/1

db_paths <- c(
  # first full run done by Kyle in 2023 (new dynamic eind implementation, 
  # and new fire equation, co2 water use efficiency adjustment is turned off)
  "fire1_eind1_c4grass1_co20" =
    file.path(path_sw, "WildfireJuly2023Runs/Output_Compiled.sqlite"),
  # co2 functionality turned on
  "fire1_eind1_c4grass1_co21" =
    file.path(path_sw, "Wildfire.CO2.August2023Runs/Output_fire1_grazL_eind1_c4grass1_co21.sqlite"),
  # no co2, no co2
  "fire0_eind1_c4grass1_co20" =
    file.path(path_sw, "NoWildfireJuly2023Runs/Output_fire0_grazL_eind1_c4grass1_co20.sqlite"))

db_connects <- map(db_paths, function(x) dbConnect(RSQLite::SQLite(), x))

# understand db structure -------------------------------------------------

# db only contains one table
tables <- map(db_connects, dbListTables)
tables

# code needs to be adjusted if multiple tables in db
stopifnot(map_dbl(tables, length) == 1) 

# list columns
cols <- map2(db_connects, tables, dbListFields)

cols

#stopifnot(cols[[1]] == cols[[2]]) # the two tables should have exactly the same
# structure

# db queries --------------------------------------------------------------

# table is small enough to fit in memory, so just loading the whole table
# then summarizing with R below
q1 <- paste("SELECT *",
            "FROM BIOMASS",
            "WHERE Year > 100;")

# separately querying the two tables
bio1 <- map(db_connects, dbGetQuery, statement = q1)

# summarize biomass -------------------------------------------------------

# add in sites that were not included in the c4off run because they contain
# c4 grasses under current conditions
#sites <- unique(bio1$c4on$site)

#stopifnot(sites %in% 1:200) # check

# c4sites <- sites[!sites %in% bio1$c4off$site]
# stopifnot(length(c4sites) == 102)

bio1a <- bio1
#bio1a$c4off <- bind_rows(bio1$c4off, bio1$c4on[bio1$c4on$site %in% c4sites, ])

# now should have complete set of sites
#stopifnot(sort(unique(bio1a$c4off$site)) == sites) 


# now creating single df from the different runs
bio2 <- bind_rows(bio1a, .id = "run")%>% 
  as_tibble()

names(bio2)
# character cols with > 1 unique value (i.e. to help 
# understand what will need to group by)
bio2 %>% 
  summarise(across(.cols = where(is.character), ~length(unique(.)))) %>% 
  select(which(colSums(.) > 1)) %>% # columns w/ > 1 unique character
  names()

# different sites have different species parameters (i.e. cheatgrass),
# so don't need to group by it (no site has 2 levels of SpeciesTreatment)
# additionally RgroupTreatment col is redundant to the [grazing]] intensity
# column

group_cols <- c("run", 'GCM', 'years', 'RCP', 
                'intensity', # the grazing intensity
                'site')
map(bio2[, group_cols], unique)

# * avg across years --------------------------------------------------------

# avg across years and make into long form
bio3 <- bio2 %>% 
  select(
    # numeric columns not interested in
    # only including biomass and indivs (# of individuals), columns
    # doesn't look like indivs has standard deviation in this column
    -matches("_(std)|(Pfire)|(PRstd)|(PR)|(RSize)|(graz)$"),
    -matches("StdDev"),
    -Year, # will avg across years
    -SpeciesTreatment, # sites have different spp. treatments (i.e. different
    # cheatgrass params based on location)
    -RGroupTreatment, # this variable is redundant with [grazing] intensity column
    -c(SoilTreatment, grazing) # cols with only 1 unique variable
  ) %>% 
  # avg across years
  group_by(across(all_of(group_cols))) %>% 
  summarise(across(.cols = where(is.numeric),
                   .fns = mean),
            n = n(),
            .groups = "drop") %>% 
  # long format
  pivot_longer(# cols selection needs to be updated as different columns are 
               # selected above
               cols = sagebrush:oppo_Indivs) %>% 
  # determine whether name value belongs to an indivs data or biomass data
  mutate(PFT = str_replace(name, "_Indivs", ""),
         # create 
         type = ifelse(str_detect(name, "Indivs"), "indivs", "biomass")) %>% 
  select(-name) %>% 
  # create separate number of individuals (indivs) and biomass columns
  pivot_wider(names_from = "type")


# Check: if grouping above missed a variable,
# then would expect value other than 50 years, per set of grouping variables
stopifnot(bio3$n == 50) 

bio3$n <- NULL

# save files --------------------------------------------------------------

# mean biomass across years for each site/scenario/treatment combination
# original version of this file (bio_mean_by_site-PFT.csv) was created with
# data from the 2021/2022 implementationof stepwat (old cheatgrass fire, no C02,
# and no dynamic eind implementation)
write_csv(bio3, "data_processed/site_means/bio_mean_by_site-PFT_v2.csv")


# dbDisconnect()
