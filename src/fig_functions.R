# Martin Holdrege

# Purpose: functions used in the creation of figures. Putting functions
# here that might be useful across a couple of different scripts going forward


# funs that create dfs ------------------------------------------------------


# function for summarizing df output for adding annotations to boxplot.
# this function is passed to the data argument of the geom_text() function
# in ggplot
box_anno <- function(df, var, group_by, id = "id",   mult = 0.05){
  # df--dataframe
  # var--variable of interest (string),
  # id--name of id variable (string)
  # group_by--variables to group by, PFT should come first in this vector
  # mult--how much above max y to put the text (in proportion of range)
  df %>% 
    # when calculating effect size it is sometimes infinite
    mutate(y = ifelse(is.infinite(.data[[var]]), NA_real_, .data[[var]])) %>% 
    group_by(across(all_of(group_by))) %>% 
    summarise(x = median(as.numeric(.data[[id]])),
              y = max(y, na.rm = TRUE) + 
                (max(y, na.rm = TRUE) - min(y, na.rm = TRUE))*mult,
              .groups = "drop_last") %>% 
    # for this to work group_by needs to be in the correct order (PFT first
    # followed by the 2nd variable of interest)
    mutate(y = max(y))
}



# axis functions ----------------------------------------------------------

# from 'id' column get the label to use for time period in boxplots
years2lab <- function(x) {
  # x--id (character vector), that includes the time period in it
  # returns--character vector of letter designations of the time period
  out <- case_when(
    str_detect(x, "Current") ~ "C",
    str_detect(x, "_2030-2060") ~ "M",
    str_detect(x, "_2070-2100") ~ "E"
  )
  if(any(is.na(out))) {
    warning("Not all time periods matched")
  }
  out
}


# label functions ---------------------------------------------------------

# function to create description based on which model simulation is being used
c4on_off_lab <- function(x) {
  stopifnot(length(x) == 1)
  
  out <- if (x == "c4on") {
    "Data from (normal) simulation where C4Pgrass site suitability functionality on."
  } else if (x == "c4off") {
    "Data from simulation where C4Pgrass restricted to 102 sites where it is currently present."
  } else {
    stop("Input incorrect")
  }
  out
}
