# Martin Holdrege

# started 1/6/2022

# Purpose: parameters (colors, axis titles etc) that are used in other scripts.
# By putting them here colors etc. can be consistent across scripts, and 
# easily changed.

# also putting some functions here that are useful when making figures

# colors ------------------------------------------------------------------

# colors for RCP scenarios
cols_rcp <- c("Current" = "cornflower blue",
              "RCP4.5" = "darkgreen",
              "RCP8.5" = "green3")

# colors for grazing levels (from 5-class OrRd on colorbrewer2)
cols_graze <- c("Light" = "#fdcc8a", 
                "Moderate" = "#fc8d59",
                "Heavy" = "#e34a33",
                "Very Heavy" = "#b30000")

# color for whether c4 present at a given point location
cols_c4present <- c("absent" = "black", "present" = "#377eb8")

# * for maps --------------------------------------------------------------

# Create color ramp that represents the proportional changes in biomass...
cols_map_bio_d <- c(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c",
                "#bd0026")),
          "#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6",
          "#08519c","#08306b")

# Standard color ramp for current biomass:
cols_map_bio <- colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d",
                           "#238443","#005a32"), interpolate = "spline")

# axis labels -------------------------------------------------------------
# meaning of endings: 0= raw value, 1 = raw change, 2 = % change

lab_yrs <- "Time period"


# biomass
lab_bio0 <- expression("Biomass ("*gm^-2*")")
lab_bio2 <- expression(~Delta*" Biomass (scaled %)")
lab_c4on0 <- expression("Biomass ("*gm^-2*"); C4Pgrass on")
lab_c4off0 <- expression("Biomass ("*gm^-2*"); C4Pgrass off")

# number of individuals
lab_indivs0 <- expression("# of individuals ("*m^-2*")")
lab_indivs2 <- expression(~Delta*" # of individuals (scaled %)")
lab_c4on0_indivs <- expression("# of individuals ("*m^-2*"); C4Pgrass on")
lab_c4off0_indivs <- expression("# of individuals ("*m^-2*"); C4Pgrass off")

# C3Pgrass/Pgrass
lab_C3_Pgrass_ratio <- expression(frac(C3Pgrass~biomass, Pgrass~biomass))

# biomass threshold
lab_below0 <- "% of sites above threshold"

# RCP
lab_rcp <- "Scenario"

#fire
lab_fire0 <- "Fire return interval (years)"
lab_fire1 <- "Change in fire return interval (years)"

# climate
lab_map <- "MAP (mm)"
lab_mat <- expression("MAT ("*degree*C*")")

# effect size
lab_es0 <- expression("ln("*frac(biomass, reference~group~biomass)*")")
lab_es0_indivs <- expression("ln("*frac(individuals, reference~group~individuals)*")")
lab_es_on_off <- expression("ln("*frac(biomass~C4off, biomass~C4on)*")")
lab_es_on_off_indivs <- expression("ln("*frac(individuals~C4off, individuals~C4on)*")")

# mapping figure dimensions ------------------------------------------------

# parameters used for layout
mar <- c(1,1,2,1)
mgp = c(3,0.3,0)

# parameters specific to 6 paneled figures
layout.matrix6 <- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, 
                         byrow = TRUE) # for 6 paneled figs

# leaving right column empty, but still same dimensions as for 6 columns
layout.matrix4 <- matrix(c(1,2, 0,3,4,0), nrow = 2, ncol = 3, 
                         byrow = TRUE) # for 4 paneled figs
widths6 <- rep(1, 3)
heights6 = rep(1, 2)

wfig6 <- 9 # width of 9 figure panels (inches)
hfig6 <- 6 # height of 9 figure panels
