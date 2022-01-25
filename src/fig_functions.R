# Martin Holdrege

# Purpose: functions used in the creation of figures. Putting functions
# here that might be useful across a couple of different scripts going forward


# funs that misc create dfs ------------------------------------------------------


# function for summarizing df output for adding annotations to boxplot.
# this function is passed to the data argument of the geom_text() function
# in ggplot
box_anno <- function(df, var, group_by, id = "id",   mult = 0.05, y = NULL,
                     anno_same_across_panels = FALSE){
  # df--dataframe
  # var--variable of interest (string),
  # id--name of id variable (string)
  # group_by--variables to group by, PFT should come first in this vector
  # mult--how much above max y to put the text (in proportion of range)
  # y--number to specify the y coordinate of the text (usually you 
  # just let the data decide that)
  # anno_same_across_pft--logical if TRUE annotation will be at the same height 
  # for all panels
  
  out <- df %>% 
    # when calculating effect size it is sometimes infinite
    mutate(y = ifelse(is.infinite(.data[[var]]), NA_real_, .data[[var]])) %>% 
    group_by(across(all_of(group_by))) %>% 
    summarise(x = median(as.numeric(.data[[id]])),
              range = max(y, na.rm = TRUE) - min(y, na.rm = TRUE),
              y = max(y, na.rm = TRUE) +  range*mult,
              # add fixed amount in case all points are of one value
              y = ifelse(range == 0, y + 0.05, y),
              .groups = "drop_last") %>% 
    # for this to work group_by needs to be in the correct order (PFT first
    # followed by the 2nd variable of interest)
    mutate(y = max(y))
  
  # so annotation are the same across panels (ie when scales are fixed)
  if(anno_same_across_panels) {
    out <- out %>% 
      ungroup() %>% 
      mutate(y = max(y))
  }
  
  if(!is.null(y)){
    out$y = y
  }
  out
}


# funs for boxplots -------------------------------------------------------


#' compute boxplot statistics
#'
#' @param df dataframe, usually grouped
#'
#' @return dataframe with the grouping variables, and values needed to 
#' draw a boxplot. This is useful for correctly plotting a boxplot w/o
#' showing the outliers
compute_boxplot_stats <- function(df, var) {
  df %>% 
    summarise(
      # boxplot stats
      # this returns 5 values
      .s = boxplot.stats(.data[[var]])$stats,
      .groups = "keep") %>% 
    # .s is a vector of length 5 for each grouping, so extracting
    # the necessary parts
    summarise(
      ymin = .s[1],
      lower = .s[2],
      middle = .s[3],
      upper = .s[4],
      ymax = .s[5],
      .groups = "drop")
}

# function for converting compute boxplot statistics() output
# into long for use in box_anno
boxplot_stats_long <- function(df) {
  df %>% 
    pivot_longer(cols = c("ymin", "lower", "middle", "upper", "ymax"),
                 names_to = "stat",
                 values_to = "y")
}

# for when creating boxplots based on compute_boxplot_stats() output
geom_boxplot_identity <- function(...) {
  geom_boxplot(stat = "identity", 
               aes(lower=lower, upper=upper, middle=middle, 
                   ymin=ymin, ymax=ymax),
               ...)
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

# for adding a second y axis to effect size figures, that shows % change
add_sec_axis <- function(...) {
  scale_y_continuous(sec.axis = sec_axis(trans = es2pchange, 
                                         name = "% Change",
                                         ...)) 
}

# label functions ---------------------------------------------------------

# function to create description based on which model simulation is being used
c4on_off_lab <- function(x) {
  stopifnot(length(x) == 1)
  
  out <- if (x == "c4on") {
    "Data from (normal) simulation where C4Pgrass expansion on."
  } else if (x == "c4off") {
    "Data from simulation where C4Pgrass expansion off"
  } else {
    stop("Input incorrect")
  }
  out
}

# convert effect size to percent change
# assumes effect size, was calculated as ln(trmt/ctrl). For use to create
# secondary y axis.
es2pchange <- function(x, round = FALSE, digits = 0) {
  out <- (exp(x) - 1)*100
  
  if(round) {
    out <- round(out, digits = digits)
  }
  
  out
}

# testing es2pchange fun
if (FALSE){
  trmt <- 1:10
  ctrl <- trmt - 0.5
  pchange <- (trmt - ctrl)/ctrl *100 # % change
  es <- log(trmt/ctrl) # effect size
  all.equal(pchange, es2pchange(es, round = FALSE)) # should be TRUE if function working
}


