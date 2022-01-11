# Martin Holdrege


# functions that assign factor levels -------------------------------------

# convert grazing intensity to a factor
graze2factor <- function(x) {
  # x--input vector of grazing levels (either numeric or already the
  # correct strings)
  
  levels <- c("24", "41", "58", "74")
  
  labels <- c("Light", "Moderate", "Heavy", "Very Heavy")
  
  # sometimes x may be a character vector (of numbers) so this way the
  # function will work if input is the numeric or char version
  x2 <- as.character(x) 
  
  # add space back in (useful when splitting apart id variable)
  if (any(x2 == "VeryHeavy")) {
    x2 <- ifelse(x2 == "VeryHeavy", "Very Heavy", x2)
  }
  
  if(all(x2 %in% levels)) {
    out <- factor(x2, levels = levels,
                  labels = labels)
  } else if (all(x2 %in% labels)) {
    out <- factor(x2, levels = labels)
  } else {
    stop("x doesn't have the right levels")
  }

  
  out
}


# converting time period (e.g. mid, vs. end of century)
# into a factor
years2factor <- function(x) {
  # x--input vector of time periods (either in form of e.g.
  # d50yrs or 2030-2060)
  # output: factor with levels of "Current", "2030-2060", "2070-2100"
  
  x2 <- x
  
  x2[is.na(x2)] <- "Current" # in model output, NA is 1981-2010
  
  levels <- c("Current", "d50yrs", "d90yrs")
  labels = c("Current", "2030-2060", "2070-2100")
  
  # assign factor levels
  if(all(x2 %in% levels)) {
    out <- factor(x2, levels = levels,
                  labels = labels)
  } else if (all(x2 %in% labels)) {
    out <- factor(x2, levels = labels)
  } else {
    stop("x doesn't have the right levels")
  }

  out
}

# converting RCP to a factor
rcp2factor <- function(x) {
  # x--input vector of time periods
  x2 <- x
  
  x2[is.na(x2)] <- "Current" # in model output NA is the current time period
  
  levels <- c("Current", "RCP45", "RCP85")
  labels = c("Current", "RCP4.5", "RCP8.5")

  # assign factor levels
  if(all(x2 %in% levels)) {
    out <- factor(x2, levels = levels,
                  labels = labels)
  } else if (all(x2 %in% labels)) {
    out <- factor(x2, levels = labels)
  } else {
    stop("x doesn't have the right levels")
  }
  
  out
}

# Rename plant functional types into the 5 main categories used in M.E.'s 
# chapter 2

pft5_factor <- function(x) {
  # x --character vector of plant functional types
  
  levels <- c("Sagebrush", 
              "C3Pgrass",
              "C4Pgrass",
              "Cheatgrass",
              "Pforb")
  
  # using this if statement so can also have the input of this function
  # be the correct characters, but still what to convert it to an ordered
  # factor
  if(all(x %in% levels)) {
    out <- x

  } else {
    out <- case_when(
      x == "sagebrush" ~ "Sagebrush",
      x %in% c("p.cool.forb", "p.warm.forb") ~ "Pforb",
      x == "p.cool.grass" ~ "C3Pgrass",
      x == "p.warm.grass" ~ "C4Pgrass",
      x == "a.cool.grass" ~ "Cheatgrass",
      TRUE ~ NA_character_
    )
  }

  out <- factor(out,levels = levels)
  
  if (!all(levels %in% out)) {
    warning("Some PFT levels are missing")
  }

  out
}


# calculate percent change ------------------------------------------------

#' Calculate scaled percent change in biomass relative to current conditions
#' 
#' @description By default this function compares change relative to Current
#'  conditions for the same grazing treatment. However, a specific level of grazing
#'  can be chosen to take the difference from. 
#'  This function could be updated to change the reference condition
#'  to any given time period and 
#'
#' @param df --grouped dataframe
#' @param var the name of the column that taking the difference of (a string),
#' default is biomass
#' @param by --variables to group by when calculating max current biomass
#' @param ref_graze character string of the value of graze column to compare
#' to (if NULL then differences are made within grazing levels)
#' @param percent logical, if TRUE, calculates the scaled percent change,
#' otherwise calculates the actual change (raw difference)
#'
#' @return dataframe of percent change in biomass from current conditions,
#' scaled by maximum current biomass
scaled_change <- function(df, 
                          var = "biomass", 
                          by = c("PFT", "graze"),
                          ref_graze = NULL,
                          percent = TRUE) {
  stopifnot(
    is.data.frame(df),
    c("RCP", "years", "graze", by) %in% names(df)
    )
  
  if(is.null(ref_graze)) {
    current <- df %>% 
      filter(RCP == "Current")
  } else if (is.character(ref_graze) & length(ref_graze) == 1 &
             ! "graze" %in% by) {
    current <- df %>% 
      filter(RCP == "Current" & graze == ref_graze)
  } else {
    stop("ref_graze arg must be null or character string, and if not null
         'by' arg must not include 'graze'")
  }
  
  current <- current %>% 
    ungroup() %>% 
    group_by(across(all_of(by))) %>% 
    mutate(current = .data[[var]])
    
    
  current_max <- current %>% 
    summarize(max_value = max(current, na.rm = TRUE),
              .groups = "drop")
  
 
  diff_var <- paste0(var,"_diff")
  out <- left_join(df, current_max, by = by) %>%
    # joining in current so making sure are subtracting site by site
    left_join(current[, c(by, "site", "current")]) %>% 
    # percent change scaled by max biomass across sites
    # or if percent is FALSE, then just the actual difference
    ungroup() %>% 
    # note using ifelse() instead of if() here caused problems--
    mutate(!!diff_var := 
             if(percent) { 
               # % scaled change
               (.data[[var]] - .data$current)/.data$max_value*100
               } else {
                 # raw change
                 .data[[var]] - .data$current
                 }
           ) %>% 
    select(-max_value) 
  
  # removing rows for reference group
  if(is.null(ref_graze)) {
    out <- out %>% 
      filter(RCP != "Current")
  } else {
    out <- out %>% 
      filter(!(RCP == "Current" & .data$graze == ref_graze))
  }
    

  
  # checking for expected number of rows
  n_current <- if (is.null(ref_graze)) {
    sum(df$RCP == "Current") 
  } else {
    sum(df$RCP == "Current" & df$graze == ref_graze) 
    }
  
  if(nrow(df) - n_current != nrow(out)) {
    stop("Output of function has in incorrect number of rows")
  }
  
  # if this warning is thrown, consider updating function (e.g. Nan, is caused
  # by dividing by 0 so maybe replace the percent diff with 100)
  if(any(is.nan(out[[diff_var]]) | is.na(out[[diff_var]]))) {
    warning("Some diffs are NA or NaN, consider fixing output")
  }
  
  # shortening name for reverse compatibility
  if("biomass_diff" %in% names(out)) {
    out <- out %>% 
      rename(bio_diff = biomass_diff)
  }
  
  out

}

