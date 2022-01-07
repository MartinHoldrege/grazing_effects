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


# * for maps --------------------------------------------------------------

# Create color ramp that represents the proportional changes in biomass...
cols_map_bio_d <- c(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c",
                "#bd0026","#800026")),
          "darkgrey", "#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6",
          "#2171b5","#08519c","#08306b")

# Standard color ramp for current biomass:
cols_map_bio <- colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d",
                           "#238443","#005a32"), interpolate = "spline")

# axis labels -------------------------------------------------------------

lab_yrs <- "Time period"


# biomass
lab_bio0 <- expression("Biomass ("*gm^-2*")")
lab_bio2 <- "Scaled % change in biomass" # 2 will denote the vers of label that is % change

# RCP
lab_rcp <- "Scenario"

# functions ---------------------------------------------------------------

# from 'id' column get the label to use for time period in boxplots
years2lab <- function(x) {
  # x--id (character vector), that includes the time period in it
  # returns--character vector of letter designations of the time period
  out <- case_when(
    str_detect(x, "_Current_") ~ "C",
    str_detect(x, "_2030-2060") ~ "M",
    str_detect(x, "_2070-2100") ~ "E"
  )
  if(any(is.na(out))) {
    warning("Not all time periods matched")
  }
  out
}


