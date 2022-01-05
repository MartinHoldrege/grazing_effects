# Martin Holdrege

# Script started 1/4/2022

# Purpose: to access the database of STEPWAT2 biomass output, and
# create summary df's for later use. 
# Simulations were done with 4 grazing levels, 200 sites, and 52 GCM*RCP*time 
# period combinations, with the cheatgrass/fire cycle turned on. 
# simulations run for 150 years, last 50 yrs used.


# dependencies ------------------------------------------------------------

library(tidyverse)
library(readr)
library(DBI)

# connect to db -----------------------------------------------------------

db_path <- "data_raw/Output.biomass.200sites.grazing.cheatgrass.Nov2021.sqlite"

db_connects <- dbConnect(RSQLite::SQLite(), db_path)

# understand db structure -------------------------------------------------

# db only contains one table
table <- dbListTables(db_connects)
table

# code needs to be adjusted if multiple tables in db
stopifnot(length(table) == 1) 

# list columns
cols <- dbListFields(db_connects, table)

cols

# db queries --------------------------------------------------------------

# table is small enough to fit in memory, so just loading the whole table
# then summarizing with R below
q1 <- paste("SELECT *",
            "FROM BIOMASS",
            "WHERE Year > 100;")

bio1 <- dbGetQuery(db_connects, q1)

# summarize biomass -------------------------------------------------------

names(bio1)
str(bio1)

bio2 <- bio1 %>% 
  as_tibble()

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

group_cols <- c('GCM', 'years', 'RCP', 'intensity', 'site')
map(bio2[, group_cols], unique)

# * avg across years --------------------------------------------------------

# avg across years and make into long form
bio3 <- bio2 %>% 
  select(
    # numeric columns not interested in
    -matches("_(std)|(Indivs)|(Pfire)|(PRstd)|(PR)|(RSize)$"),
    -matches("StdDev"),
    -Year, # will avg across years
    -SpeciesTreatment, # sites have different spp. treatments (i.e. different
    # cheatgrass params based on location)
    -RGroupTreatment, # this variable is redundant with [grazing] intensity column
    -c(SoilTreatment, dst, grazing) # cols with only 1 unique variable
  ) %>% 
  # avg across years
  group_by(across(all_of(group_cols))) %>% 
  summarise(across(.cols = where(is.numeric),
                   .fns = mean),
            n = n(),
            .groups = "drop") %>% 
  # long format
  pivot_longer(cols = sagebrush:oppo,
               names_to = "PFT",
               values_to = "biomass")

# Check: if grouping above missed a variable,
# then would expect value other than 50 years, per set of grouping variables
stopifnot(bio3$n == 50) 

bio3$n <- NULL

# save files --------------------------------------------------------------

# mean biomass across years for each site/scenario/treatment combination
write_csv(bio3, "data_processed/site_means/bio_mean_by_site-PFT.csv")


