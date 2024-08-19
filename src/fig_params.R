# Martin Holdrege

# started 1/6/2022

# Purpose: parameters (colors, axis titles etc) that are used in other scripts.
# By putting them here colors etc. can be consistent across scripts, and 
# easily changed.

# also putting some functions here that are useful when making figures

# colors ------------------------------------------------------------------

# colors for RCP scenarios
cols_rcp <- c("Current" = "cornflower blue",
              "RCP45" = "darkgreen",
              "RCP85" = "green3")

# colors for grazing levels (from 5-class OrRd on colorbrewer2)
cols_graze <- c("Light" = "#fdd49e", 
                "Moderate" = "#fc8d59",
                "Heavy" = "#d7301f",
                "Very Heavy" = "#7f0000")

# color for whether c4 present at a given point location
cols_c4present <- c("absent" = "black", "present" = "#377eb8")

# colors for the main 5 PFts, and an 'other' category
cols_pft5_other <- c("Other" = "darkgrey",
                     "Cheatgrass" = "#f6e8c3",
                     "Pforb" = "#d8b365",
                     "C4Pgrass" = "#c7eae5",
                     "C3Pgrass" = "#5ab4ac",
                     "Sagebrush" = "#01665e")

cols_herb <- cols_pft5_other[c("Cheatgrass","Pforb", "C4Pgrass", "C3Pgrass")]

# fire probability
cols_firep <- RColorBrewer::brewer.pal(11, "RdYlBu")[9:2]

# * for maps --------------------------------------------------------------

# Create color ramp that represents the proportional changes in biomass...
cols_map_bio_d <- c(rev(c("#ffeda0","#fed976","#feb24c","#fc4e2a","#e31a1c",
                "#bd0026")), 'grey', 'grey',
          "#deebf7","#c6dbef","#9ecae1","#4292c6",
          "#08519c","#08306b")

cols_map_bio_dr <- colorRampPalette(cols_map_bio_d, interpolate = 'spline')

# Standard color ramp for current biomass:
cols_map_bio <- colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d",
                           "#238443","#005a32"), interpolate = "spline")

# axis labels -------------------------------------------------------------
# meaning of endings: 0= raw value, 1 = raw change, 2 = % change

lab_yrs <- "Time period"

# change
lab_change0 <- "% Change relative to light grazing"

# biomass
lab_bio0 <- expression("Biomass ("*gm^-2*")")
lab_bio1 <- expression(~Delta*" Biomass ("*gm^-2*")")
lab_bio1b <- "\u0394 Biomass" # delta biomass
lab_bio2 <- expression(~Delta*" Biomass (scaled %)")
lab_c4on0 <- expression("Biomass ("*gm^-2*"); C4Pgrass on")
lab_c4off0 <- expression("Biomass ("*gm^-2*"); C4Pgrass off")
lab_perc_bio0 <- "% of total biomass"

# utilization
lab_util0 <- expression("Utilization ("*gm^-2*")")
lab_util2 <- expression(~Delta*" Utilization (scaled %)")

# cover
lab_cov0 <-  "Cover (%)"
lab_cov1 <- expression(~Delta*" Cover (percentage point)")
lab_cov2 <- expression("% Change in cover")

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

# grazing
lab_graze <- "Grazing"

#fire (return interval)
lab_fire0 <- "Fire return interval (years)"
lab_fire1 <- "Change in fire return interval (years)"

# fire probability
lab_firep0 <- "# fires/century"
lab_firep1 <- expression(Delta~"# fires/century")

# climate
lab_map <- "MAP (mm)"
lab_mat <- expression("MAT ("*degree*C*")")
lab_corrtp <- 'Seasonality (T-P correlation)'
lab_psp <- 'Prop summer ppt'
lab_map0 <- lab_map
lab_mat0 <- lab_mat
lab_map1 <- expression(Delta~"MAP (mm)")
lab_mat1 <- expression(Delta~"MAT ("*degree*C*")")

# transpiration
lab_transp0 <- 'Transpiration (mm)'
lab_transp02 <- 'Transpiration (mm/cm soil)'

# effect size
lab_es0 <- expression("ln("*frac(biomass, reference~group~biomass)*")")
lab_es1 <- expression("ln("*frac(biomass, light~grazing~biomass)*")")
lab_es1_heavy <- expression("ln("*frac(biomass, heavy~grazing~biomass)*")")
lab_es0_indivs <- expression("ln("*frac(individuals, reference~group~individuals)*")")
lab_es0_util <- expression("ln("*frac(utilization, reference~group~utilization)*")")
lab_es1_util_heavy <- expression("ln("*frac(utilization, heavy~grazing~utilization)*")")
lab_es_on_off <- expression("ln("*frac(biomass~C4off, biomass~C4on)*")")
lab_es_on_off_indivs <- expression("ln("*frac(individuals~C4off, individuals~C4on)*")")

# SEI
lab_sei1 <- expression(Delta~"SEI")

# Q
lab_q1 <- expression(Delta~"Q")

# area
lab_areakm0 <- expression("Area ("*km^2*")")
lab_areaha0 <- "Area (million ha)"
lab_areaperc0 <- "Area (%)"

# mapping figure dimensions ------------------------------------------------

# parameters used for layout
mar <- c(1,1,2,1)
mgp = c(3,0.3,0)

# parameters specific to 6 paneled figures
layout.matrix6 <- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, 
                         byrow = TRUE) # for 6 paneled figs

# parameters specific to 8 paneled figures
layout.matrix9 <- matrix(1:9, nrow = 3, ncol = 3, 
                         byrow = TRUE) # for 6 paneled figs

# leaving right column empty, but still same dimensions as for 6 columns
layout.matrix4 <- matrix(c(1,2, 0,3,4,0), nrow = 2, ncol = 3, 
                         byrow = TRUE) # for 4 paneled figs
widths6 <- rep(1, 3)
heights6 = rep(1, 2)
widths9 <- rep(1, 3)
heights9 <- rep(1, 3)

wfig6 <- 9 # width of 9 figure panels (inches)
hfig6 <- 6 # height of 9 figure panels
