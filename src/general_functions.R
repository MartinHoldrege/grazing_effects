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

#' Calculate scaled percent change in biomass
#'
#' @param df --grouped dataframe
#' @param by --variables to group by when calculating max current biomass
#'
#' @return dataframe of percent change in biomass from current conditions,
#' scaled by maximum current biomass
scaled_change <- function(df,  by = c("PFT", "graze")) {
  stopifnot(
    is_grouped_df(df),# should be already grouped
    c("biomass", "RCP", "years", "PFT", "graze") %in% names(df)
    )
  current_mean <- df %>% 
    ungroup() %>% 
    group_by(across(all_of(by))) %>% 
    filter(RCP == "Current") %>% 
    summarize(max_biomass = max(biomass),
              .groups = "drop")
  
  out <- left_join(df, current_mean, by = by) %>% 
    # percent change scaled by max biomass across sites
    mutate(bio_diff = (biomass - biomass[RCP == "Current"])/max_biomass*100) %>% 
    select(-max_biomass) %>% 
    filter(RCP != "Current")
  
  # checks
  n_current <- sum(df$RCP == "Current")
  if(nrow(df) - n_current != nrow(out)) {
    stop("Output of function has in incorrect number of rows")
  }
  # if this warning is thrown, consider updating function (e.g. Nan, is caused
  # by dividing by 0 so maybe replace the percent diff with 100)
  if(any(is.nan(out$bio_diff) | is.na(out$bio_diff))) {
    warning("Some bio diff are NA or NaN, consider fixing output")
  }
  
  out

}


