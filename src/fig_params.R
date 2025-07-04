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

cols_pft3 <- c("Aherb" = "#f6e8c3",
               "Pherb" = "#c7eae5",
               "Sagebrush" = "#01665e")


# fire probability
cols_firep <- RColorBrewer::brewer.pal(11, "RdYlBu")[9:2]


# colors from the SEI repository, same colors as used in appendix
# of Holdrege et al 2024 REM, cool colors are GCMs with lower temperatures
# across 200 stepwat sites, and warmer colors are the hotter GCMs 
cols_GCM1 <- c(inmcm4 = "#313695", `GISS-E2-R` = "#416AAE", `MRI-CGCM3` = "#649AC7", 
               `FGOALS-s2` = "#8FC3DC", `FGOALS-g2` = "#BCE1ED", `CSIRO-Mk3-6-0` = "#E5F4EE", 
               `CESM1-CAM5` = "#FFFFBF", MIROC5 = "#FEE597", `IPSL-CM5A-MR` = "#FDBE70", 
               `HadGEM2-CC` = "#F88D52", CanESM2 = "#EA5839", `HadGEM2-ES` = "#CE2826", 
               `MIROC-ESM` = "#A50026")

cols_GCM2 <-  c('Historical' = 'black', cols_GCM1)


cols_ecoregion <- c("Great Plains" = "#ECA130",
                    "Intermountain West" = "#33A37A",
                    "Southern Great Basin" = "#1673B2",
                    'Western Intermountain' = "#41ae76",
                    "Eastern Intermountain" = "#005824")

# continuous SEI
cols_seicont = c("#eee1ba", "#a6611a", "#78c679", "#006837", "#2166ac", "#053061")
# these are the SEI category break points
vals_seicont = c(0, 0.173, 0.174, 0.431, 0.432, 1)

# colors for summary across GCMs
cols_smry <-   c(low = "#ffffcc", 'median' = "#a1dab4",
                 'high' = "#2c7fb8")

# for different age groups or 'successional states'
cols_agegroup <- RColorBrewer::brewer.pal(4, 'BrBG')

# * SEI -------------------------------------------------------------------

c9Names <-  c(
  'Stable CSA',
  'CSA becomes GOA',
  'CSA becomes ORA',
  'GOA becomes CSA',
  'Stable GOA',
  'GOA becomes ORA',
  'ORA becomes CSA',
  'ORA becomes GOA',
  'Stable ORA'
)

c3Names <- c("CSA", "GOA", "ORA")

c9Palette = c('#142b65', # stable core (black)
              '#b30000', #'#d7301f', # core becomes grow # reds from 9-class OrRd
              '#67001f',  # core becomes other
              '#757170', # grow becomes core
              '#99d4e7', # stable grow
              '#fc8d59', # grow becomes other
              '#000000', # other becomes core
              '#D9D9D9', # other becomes grow
              '#eee1ba') # stable other


names(c9Palette) <- c9Names

c3Palette <- c9Palette[c(1, 5, 9)]
names(c3Palette) <- c3Names
cols_c3 <- c3Palette

# * for maps --------------------------------------------------------------

# Create color ramp that represents the proportional changes in biomass...
cols_map_bio_d <- c(rev(c("#ffeda0","#fed976","#feb24c","#fc4e2a","#e31a1c",
                "#bd0026")), 'grey', 'grey',
          "#deebf7","#c6dbef","#9ecae1","#4292c6",
          "#08519c","#08306b")

index <- which(cols_map_bio_d == 'grey')[1]
cols_map_bio_d2 <- cols_map_bio_d[-index] # shorter gray area when continuous

cols_map_bio_dr <- colorRampPalette(cols_map_bio_d, interpolate = 'spline')

# Standard color ramp for current biomass:
cols_map_bio <- colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d",
                           "#238443","#005a32"), interpolate = "spline")

# colors for predictor variables of fire equation
cols_pred <- list(
  mat = c(low = '#fee0d2', high = '#67000d'),
  map = c(low = '#d0d1e6', high = '#023858'),
  psp = c(low = '#d0d1e6', high = '#023858'),
  pfg = c(low = '#f7fcb9', high = '#004529'),
  afg = c(low = '#f7fcb9', high = '#004529')
)

# for maps showing primary driver of change
cols_pred_vars <- c(MAP = "#1f78b4", MAT = "#e31a1c", PSP = "#a6cee3", 
                    Aherb = "#b2df8a", Pherb = "#33a02c")


# shapes ------------------------------------------------------------------

shapes_scen <- c("Historical" = 21, 
                 "RCP4.5\n2071-2100" = 17, "RCP4.5 2071-2100" = 17,
                 "RCP8.5\n2071-2100" = 15, "RCP8.5 2071-2100" = 15)


# linetypes ---------------------------------------------------------------

linetypes_scen <- c("Historical" = 4, 
                 "RCP4.5\n2071-2100" = 1, "RCP4.5 2071-2100" = 1,
                 "RCP8.5\n2071-2100" = 41, "RCP8.5 2071-2100" = 41)

# panel letters -----------------------------------------------------------

# so consistent letter theme is used throughout (for multiple
# panel pub quality figs)
fig_letters <- paste(letters, ")", sep = "")
names(fig_letters) <- letters

# axis labels -------------------------------------------------------------
# meaning of endings: 0= raw value, 1 = raw change, 2 = % change

lab_yrs <- "Time period"

# change
lab_change0 <- "% Change relative to light grazing"

# biomass
lab_bio0 <- expression("Biomass ("*gm^-2*")")
lab_bio1 <- expression(~Delta*" Biomass ("*gm^-2*")")
lab_bio1b <- "\u0394 Biomass" # delta biomass
lab_bio2 <- expression(~Delta*" Biomass (%)")
lab_c4on0 <- expression("Biomass ("*gm^-2*"); C4Pgrass on")
lab_c4off0 <- expression("Biomass ("*gm^-2*"); C4Pgrass off")
lab_perc_bio0 <- "% of total biomass"

# utilization
lab_util0 <- expression("Utilization ("*gm^-2*")")
lab_util0b <- "Utilization (%)"
lab_util1 <- expression(~Delta*" Utilization ("*gm^-2*")")
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

# fire area
lab_ba0 <- "Expected burned area (ha/yr)"
lab_ba1 <- "Expected burned area (%/yr)"
lab_agegroup <- "Years since fire"


# climate
lab_map <- "MAP (mm)"
lab_mat <- expression("MAT ("*degree*C*")")
lab_corrtp <- 'Seasonality (T-P correlation)'
lab_psp <- 'Prop summer ppt'
lab_psp1 <- expression(Delta~"PSP")
lab_map0 <- lab_map
lab_mat0 <- lab_mat
lab_map1 <- expression(Delta~"MAP (mm)")
lab_mat1 <- expression(Delta~"MAT ("*degree*C*")")

labs_pred1 <- list(
  'MAP' = lab_map1,
  'MAT' = lab_mat1,
  'PSP' = lab_psp1,
  'psp' = lab_psp1,
  'Pherb' = lab_bio1,
  'Aherb' = lab_bio1
)

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
lab_sei0 <- 'SEI'
lab_sei1 <- expression(Delta~"SEI")
lab_c3 <- 'SEI class'
# Q
lab_q0 <- 'Q Value'
lab_q1 <- expression(Delta~"Q")

# area
lab_areakm0 <- expression("Area ("*km^2*")")
lab_areaha0 <- "Area (million ha)"
lab_area0 <- 'Area (ha)'
lab_areaperc0 <- "Area (%)"
lab_areaperc0b <- "Area (% of region)"

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
