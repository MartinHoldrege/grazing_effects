# Martin Holdrege

# Script started 1/4/2022

# Purpose: to access the database of STEPWAT2 biomass output, and
# compare fire return interval and fire probability for whole time period of
# simulation (150 years) vs last 50 years


# dependencies ------------------------------------------------------------

library(DBI)
theme_set(theme_classic())
source("scripts/02_summarize_bio.R")
# params ------------------------------------------------------------------

n_iter <- 200 # number of iterations run


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
    file.path(path_sw, "Wildfire.CO2.August2023Runs/Output_fire1_grazL_eind1_c4grass1_co21.sqlite"))

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
            "FROM BIOMASS")

# separately querying the two tables
bio1 <- map(db_connects, dbGetQuery, statement = q1)

# summarize fire -------------------------------------------------------

bio1a <- bio1

# now creating single df from the different runs
bio2 <- bind_rows(bio1a, .id = "run")%>% 
  as_tibble()

# * avg across years --------------------------------------------------------

n_years <- max(bio2$Year) # this would be a problem if simulation lengths
# differ between runs

# avg fire return interval and fire probability
fire_allyrs <- bio2 %>% 
  select(run, site, GCM, years, RCP, intensity, WildFire) %>%
  mutate(graze = graze2factor(intensity),
         RCP = rcp2factor(RCP),
         years = years2factor(years),
         # the unique combination of treatments
         id = paste(RCP, years, graze, sep = "_"), 
         id = str_replace(id, " ", ""),
         scenario_label = ifelse(RCP == "Current", "(Historical)", paste0("(",RCP,", ",years, ")"))) %>% 
  dplyr::select(-intensity) %>% 
    group_by(run, years, RCP, id,
             graze, scenario_label,
             site, GCM) %>% # GCM last (so drop last works)
    summarize(n_fires = sum(WildFire),
              .groups = "drop_last") %>% 
    # 200 iterations occured and we're using the last 50 years of data
    # if no fires occurred then set fire return interval to NA
    mutate(fire_return = (n_years*n_iter)/n_fires,
           fire_return = ifelse(is.infinite(fire_return), NA, fire_return),
           fire_prob = n_fires/(n_years*n_iter)*100)
  
# median across GCMs
fire_allyrs_med1 <- fire_allyrs %>% 
  summarize(fire_prob = median(fire_prob),
            n_fires = median(n_fires),
            .groups = "drop")  %>% 
  mutate(fire_return = (n_years*n_iter)/n_fires) 


# combine 50 and 150 yr fire data -----------------------------------------

med2 <- fire_med1 %>% # object from 02_summarize_bio.R where fire was calculated from 50 years
  inner_join(fire_allyrs_med1, by = join_by(run, years, RCP, graze, id, site),
             suffix = c("_50", "_150"))

# figures -----------------------------------------------------------------

# fire probability in the first 
med2_l <- med2 %>% 
  group_by(run) %>% 
  group_split()

fig_l <- map(med2_l, function(df) {
  g <- ggplot(df) +
    labs(caption = paste('simulation settings:', df$run[1]),
         subtitle = 'Comparing fire data from the entire simulation (150 yrs) vs last 50 yrs') +
    geom_abline(slope = 1, intercept = 0) +
    facet_wrap(~scenario_label)
  
  g +
    geom_point(aes(fire_prob_50, fire_return_50))
  g_fri <- g + 
    geom_point(aes(fire_return_50, fire_return_150)) +
    labs(x = 'FRI (years) of last 50 years',
         y = 'FRI (years) of entire 150 years')
  
  g_fri2 <- g_fri +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 200))
  
  g_prob <- g + 
    geom_point(aes(fire_prob_50, fire_prob_150)) +
    labs(x = 'Fire probability (%) of last 50 years',
         y = 'Fire probability (%) of entire 150 years')
  
  return(list(g_prob, g_fri, g_fri2))
})

pdf("figures/fire/fire_50_vs_150_yrs_v1.pdf")
print(fig_l)
dev.off()

