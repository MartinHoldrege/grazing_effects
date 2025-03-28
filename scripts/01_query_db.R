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
source('src/paths.R') # db files stored on hard-drive
# connect to db -----------------------------------------------------------

# naming convention:
# Wildfire – 0/1
# Grazing intensity –  L, M, H, VH
# Site-specific eind – 0/1
# C4 grass expansion under future conditions – 0/1
# CO2 – 0/1
# runs done in 2025--quantile mapping for fire equation has been 
# implemented
# 2502--first runs that used quantile mapping for fire probability inputs
# 2503 runs have corrected the quantile mapping values used for fire probability


# enter details here (for constructing file paths)
run_details <- list(
  c('run' = 'fire1_eind1_c4grass1_co20_2503', 'graze' = 'grazLMHVH'),
   c('run' = 'fire1_eind1_c4grass0_co20_2503', 'graze' = 'grazLMHVH')
)

# construct file paths
db_paths <- map(run_details, function(x) {
  run2 <- str_replace(x['run'], '_\\d{4}$', 
                      paste0('_', x['graze'], "\\0"))
  
  # adding in the grazing levels to the filename
  file.path(path_sw, run2, 
            paste0('Output_', run2, '.sqlite'))
})
names(db_paths) <- map_chr(run_details, \(x) x['run'])


stopifnot(map_lgl(db_paths, file.exists))

db_connects <- map(db_paths, function(x) dbConnect(RSQLite::SQLite(), x))

# understand db structure -------------------------------------------------

# db only contains one table
tables <- map(db_connects, dbListTables)
tables

# code needs to be adjusted if multiple tables in db
stopifnot(map_dbl(tables, length) == 1) 

# list columns
cols <- map2(db_connects, tables, dbListFields)

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
bio2a <- bio2 %>% 
  # filter(SoilTreatment == 'soils_fixed1') %>% 
  select(
    # numeric columns not interested in
    # only including biomass and indivs (# of individuals), columns
    # doesn't look like indivs has standard deviation in this column
    -matches("_(std)|(Pfire)|(PRstd)|(PR)|(RSize)$"),
    -matches("StdDev"),
    -Year, # will avg across years
    -SpeciesTreatment, # sites have different spp. treatments (i.e. different
    # cheatgrass params based on location)
    -RGroupTreatment, # this variable is redundant with [grazing] intensity column
    -c(SoilTreatment, grazing) # cols with only 1 unique variable
  ) %>% 
  # avg across years
  group_by(across(all_of(group_cols))) %>% 
  summarise(across(.cols = where(is.numeric) & !matches("WildFire"),
                   .fns = mean),
            WildFire = sum(WildFire),
            n = n(),
            .groups = "drop") %>% 
  #select(-(artr:oppo)) %>% 
  # long format
  pivot_longer(# cols selection needs to be updated as different columns are 
               # selected above
               cols = sagebrush:oppo_Indivs) 


bio3a <- bio2a %>% 
  # determine whether name value belongs to an indivs data or biomass data
  mutate(PFT = str_replace(name, "(_Indivs)|(_graz)", ""),
         # create 
         type = case_when(
           str_detect(name, "Indivs") ~ "indivs",
           str_detect(name, 'graz') ~ 'utilization',
           TRUE ~ 'biomass')) %>% 
  select(-name) %>% 
  # create separate number of individuals (indivs) and biomass columns
  pivot_wider(names_from = "type")

PFT_lookup <- c(artr = "sagebrush", cryp = "a.cool.forb", chen = "a.warm.forb", 
                phho = "p.cool.forb", arfr = "p.warm.forb", brte = "a.cool.grass", 
                pssp = "p.cool.grass", bogr = "p.warm.grass", chvi = "shrub", 
                oppo = "succulents")

stopifnot(unique(bio3a$PFT) %in% c(PFT_lookup, names(PFT_lookup)))

# rows with spp codes
bio_spp <- bio3a %>% 
  filter(PFT %in% names(PFT_lookup)) %>% 
  mutate(PFT4join = PFT_lookup[PFT]) %>% 
  select(run, GCM, years, RCP, intensity, site, PFT4join, indivs)

# joining in # of indivs column (so don't have duplicated rows, one showing
# spp code with indivs, and the others show pft code and no indivs)
bio3 <- bio3a %>% 
  filter(PFT %in% PFT_lookup) %>% 
  select(-indivs) %>% 
  left_join(bio_spp, by = c("run", "GCM", "years", "RCP", "intensity", "site", 
                            PFT = "PFT4join"))
  
# Check: if grouping above missed a variable,
# then would expect value other than 50 years, per set of grouping variables
# also note--if this isn't true (e.g. if simulations run for 300 not 150 years,
# then fire return interval calculation (next script) needs to be adusted)
stopifnot(bio3$n == 50) 

bio3$n <- NULL

# save files --------------------------------------------------------------

# mean biomass across years for each site/scenario/treatment combination
# original version of this file (bio_mean_by_site-PFT.csv) was created with
# data from the 2021/2022 implementationof stepwat (old cheatgrass fire, no C02,
# and no dynamic eind implementation)
write_csv(bio3, "data_processed/site_means/bio_mean_by_site-PFT_v6.csv")

map(db_connects, dbDisconnect) # disconnect

