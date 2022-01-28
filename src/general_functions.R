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



# * PFTs ------------------------------------------------------------------

# Rename plant functional types into the 5 main categories used in M.E.'s 
# chapter 2

pft5_factor <- function(x, return_levels = FALSE) {
  # x --character vector of plant functional types
  # return_levels--logical, set to true if want function to just return
  # a vector of the 5 pft's
  
  levels <- c("Sagebrush", 
              "C3Pgrass",
              "C4Pgrass",
              "Cheatgrass",
              "Pforb")
  if(return_levels) {
    return(levels)
  }
  
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

# 3 big pft groups (excluding cheatgrass)
# especially useful when combining df based on these factor levels,
# with other df's that have other PFTs
pft3_factor <- function(x) {
  # x --character vector of plant functional types
  
  levels <- c("Shrub", 
              "Pherb",
              "Aforb")
  
  # using this if statement so can also have the input of this function
  # be the correct characters, but still what to convert it to an ordered
  # factor
  if(all(x %in% levels)) {
    out <- x
    
  } else {
    out <- case_when(
      x %in% c("shrub", "sagebrush") ~ "Shrub",
      x %in% c("p.cool.grass", "p.warm.grass",
               "p.cool.forb", "p.warm.forb") ~ "Pherb",
      x == "a.warm.forb" ~ "Aforb",
      TRUE ~ NA_character_
    )
  }
  
  out <- factor(out,levels = levels)
  
  if (!all(levels %in% out)) {
    warning("Some PFT levels are missing")
  }
  
  out
}

Pgrass_factor <- function(x) {
  # x --character vector of plant functional types
  
  levels <- "Pgrass"
  
  # using this if statement so can also have the input of this function
  # be the correct characters, but still what to convert it to an ordered
  # factor
  if(all(x %in% levels)) {
    out <- x
    
  } else {
    out <- case_when(
      x %in% c("p.cool.grass", "p.warm.grass") ~ "Pgrass",
      TRUE ~ NA_character_
    )
  }
  
  out <- factor(out,levels = levels)
  
  if (!all(levels %in% out)) {
    warning("Some PFT levels are missing")
  }
  out
}

# turning all pft levels into a factor (after they have already been
# combined based on the functions above)
pft_all_factor <- function(x) {
  levels <-  c("Sagebrush", 
              "C3Pgrass",
              "C4Pgrass",
              "Cheatgrass",
              "Pforb",
              "Shrub", 
              "Pgrass",
              "Pherb",
              "Aforb")
  
  # add a 'total' factor level, if the input includes it
  if("Total" %in% x) {
    levels <- c(levels, "Total")
  }
  
  stopifnot(x %in% levels)
  
  out <- factor(x, levels = levels)
}

# convert individual PFTs into a total category, useful for calculating
# total biomass
pft_total_factor <- function(x) {
  
  pft <- c("sagebrush", "a.cool.forb", "a.warm.forb", "p.cool.forb", 
           "p.warm.forb", "a.cool.grass", "p.cool.grass", "p.warm.grass",
           "shrub", "succulents")
  
  # x must include all these levels, otherwise won't be computing total biomass
  stopifnot(pft %in% x) 
  
  out <- ifelse(x %in% pft, "Total", NA_character_)
  out <- factor(out, levels = "Total")
  out
}

# create new variables ----------------------------------------------------

#' create id2 variable from id variable
#'
#' @param df a dataframe
#'
#' @return a dataframe with a new column id2 (a factor), which is the id variable 
#' but with grazing removed from the string
create_id2 <- function(df) {
  stopifnot(
    is.data.frame(df),
    c("RCP", "years", "id") %in% names(df)
  )
  out <- df %>% 
    arrange(RCP, years) %>% # for creating ordered factor
    # id variable grazing removed
    mutate(id2 = str_replace(id, "_[A-z]+$", ""),
           id2 = factor(id2, levels = unique(id2))) 
  out
}

#' Create reference id from from id vector
#' 
#' @description The purpose is to create an id string that is the reference
#' class that corresponds to a given id string. At the moment the reference
#' class created is current, same grazing level. For useful when subsetting vectors
#'
#' @param x character vector which is the id, including c4, RCP and years
#'
#' @return character vector, same length as x, but c4off changed
#' to c4on, and scenario changed to current. 
create_ref_id <- function(x) {
  # reference group is c4on
  str_replace(x, "^c4off_", "c4on_") %>% 
    # both changing RCP and years to current
    str_replace("_RCP\\d.\\d_\\d{4}-\\d{4}", "_Current_Current")
}


#' Remove c4 designation from id string
#'
#' @param x id string with c4on or c4off as part of the string
#'
#' @return id string w/ c4on/off removed
id_remove_c4 <- function(x){
  str_replace(x, "c4(on|off)_", "")
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
#' @param effect_size logical, whether the response variable calculated
#' is effect size instead of change (or percent change), currently calculated as 
#' ln(trmt/ctrl)
#' @param within_GCM logical, if TRUE comparisons will be made within a GCM, e.g.
#' comparing heavy to light grazing, within a GCM/RCP/years. ref_graze must be 
#' specified if set to TRUE
#'
#' @return dataframe of percent change in biomass from current conditions,
#' scaled by maximum current biomass
scaled_change <- function(df, 
                          var = "biomass", 
                          by = c("PFT", "graze"),
                          ref_graze = NULL,
                          percent = TRUE,
                          effect_size = FALSE, 
                          within_GCM = FALSE) {
  
  # checking arguments
  stopifnot(
    is.data.frame(df),
    c("RCP", "years", "graze", by) %in% names(df)
    )
  
  if(percent & effect_size) {
    stop("percent and effect_size args can't both be TRUE")
  }
  
  if(within_GCM & ! all(c("GCM", "RCP", "years") %in% by)) {
    stop("if within_GCM is TRUE the by arg needs to include GCM, RCP and years")
  }
  
  if(within_GCM & is.null(ref_graze)) {
    stop("if within_GCM is TRUE ref_graze must be provided")
  }
  if(within_GCM & "graze" %in% by) {
    stop("if within_GCM is TRUE 'by' should not included graze")
  }
  
    # creating the df that is the reference class (usually current conditions)
  if(within_GCM) {
    # comparing within GCM
    current <- df %>% 
      filter(graze == ref_graze)
  } else if (is.null(ref_graze)) {
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
  
  if (effect_size) {
    diff_var <- paste0(var,"_es") # 'es' for effect size
  } else {
    diff_var <- paste0(var,"_diff")
  }
 
  
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
             } else if (effect_size) {
               # log response ratio
               log(.data[[var]]/.data$current)
             }  else {
               # raw change
               .data[[var]] - .data$current
              }
           ) %>% 
    select(-max_value) 
  
  # removing rows for reference group
  if(within_GCM) {
    out <- out %>% 
      filter(.data$graze != ref_graze)
  } else if (is.null(ref_graze)) {
    out <- out %>% 
      filter(RCP != "Current")
  } else {
    out <- out %>% 
      filter(!(RCP == "Current" & .data$graze == ref_graze))
  }
    
  
  # checking for expected number of rows
  n_current <- if (within_GCM) {
    sum(df$graze == ref_graze)
  } else if (is.null(ref_graze)) {
    sum(df$RCP == "Current") 
  } else {
    sum(df$RCP == "Current" & df$graze == ref_graze) 
  }
  
  if(nrow(df) - n_current != nrow(out)) {
    stop("Output of function has in incorrect number of rows")
  }
  
  # if this warning is thrown, consider updating function (e.g. Nan, is caused
  # by dividing by 0 so maybe replace the percent diff with 100)
  # don't want warning for if effect_size = TRUE, b/ inevitable with effect
  # size calculations when current is 0
  if(any(is.nan(out[[diff_var]]) | is.na(out[[diff_var]])) & !effect_size) {
    warning("Some diffs are NA or NaN, consider fixing output")
  }
  
  # shortening name for reverse compatibility
  if("biomass_diff" %in% names(out)) {
    out <- out %>% 
      rename(bio_diff = biomass_diff)
  } else if ("biomass_es" %in% names(out)) {
    out <- out %>% 
      rename(bio_es = biomass_es)
  }
  
  out

}


# outliers ----------------------------------------------------------------

#' remove outliers (mostly for plotting boxplots)
#'
#' @param df dataframe, function most useful for grouped dataframes
#' @param var string, name of column to remove outliers from
#'
#' @return Data with outliers (just removed those points that would be 
#' considered outliers in a boxplot--so this isn't a rigorous outlier
#' remover, mostly for plotting), the problem is that a boxplot
#' made after outliers removed will no longer have the right quartiles
#' (especially if many outliers have been removed)
remove_outliers <- function(df, var) {
  
  stopifnot(is.data.frame(df))
  
  out <- df %>% 
    mutate(.ymax = boxplot.stats(.data[[var]])$stats[5],
           .ymin = boxplot.stats(.data[[var]])$stats[1],
           !!var := ifelse(.data[[var]] <= .ymax & .data[[var]] >= .ymin,
                           .data[[var]], NA_real_)) %>% 
    select(-.ymax, -.ymin)
  out
}




