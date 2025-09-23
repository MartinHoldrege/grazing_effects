# Martin Holdrege


# functions that assign factor levels -------------------------------------

# function factory
var2factor_factory <- function(levels) {
  function(x, return_levels = FALSE) {
    if(return_levels) {
      return(levels)
    }
    stopifnot(x %in% levels)
    factor(x, levels)
  }
}

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

# convert grazing names into naming convention that is appropriate
# for convention used in filenames
graze2filename <- function(x) {
  
  y <- graze2factor(c("24", "41", "58", "74"))
  stopifnot(all(x %in% y))
  
  factor(x, levels = levels(y),
         labels = c('grazL', 'grazM', 'grazH', 'grazVH'))

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
  labels = c("Current", "RCP45", "RCP85")

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

epoch2factor <- function(x) {
  x2 <- update_yr(x) 
  factor(x2, levels = c('2031-2060', '2071-2100'))
}

# put '.' in RCP names
update_rcp <- function(x) {
  
  if(is.factor(x)) {
    levels <- levels(x)
    labels_new <- levels %>% 
      str_replace('RCP45', 'RCP4.5') %>% 
      str_replace('RCP85', 'RCP8.5') %>% 
      str_replace('Current', 'Historical')
    
    out <- factor(as.character(x), levels = levels,
                  labels = labels_new)
 
  } else {
    out <- x %>% 
      str_replace('RCP45', 'RCP4.5') %>% 
      str_replace('RCP85', 'RCP8.5') %>% 
      str_replace('Current', 'Historical')
  }
  out

}

# doesn't create a factor, but just to update the correct
# epoch label for use in figures
update_yr <- function(x) {
  x %>% 
    str_replace('2030', '2031') %>% 
    str_replace('2070', '2071')
}

# cut depth into 3 categories
cut_depth <- function(x, two_depths = FALSE) {
  # x--numeric vector, depth of the layer
  # two_depths --logical of whether to create two depth categories
  # (default is 3--for back compatibility)
  stopifnot(x %in% lyr2depth(1:8))
  out <- if (two_depths) {
    cut(x,
        breaks = c(0, 10, 200),
        labels = c("0-10 cm", "10-150 cm"))
  } else {
    cut(x,
        breaks = c(0, 10, 40, 200),
        labels = c("0-10 cm", "10-40 cm", "40-150 cm"))
  }
  out
}

all_region_levels <- function() {
  lev_l <- list(
    'wafwa' = c("Great Plains", "Intermountain West", 
                "Southern Great Basin"),
    r1.0 = c("Western Intermountain",
             "Great Plains",
             "Southern Great Basin",
             'Eastern Intermountain'),
    r1.1 = c("NW Deserts", "NW Mountains", "Great Plains", "Southern Great Basin", 
             "Middle Rockies", "WY Basin", "Plateaus", "Southern Mountains")
  )
  lev_l
}

# determine which version of regions
which_vr <- function(x) {
  # removing the entire study area for next matching step
  x <- str_subset(x, '[Ee]ntire', negate = TRUE)
  l <- all_region_levels()
  
  lgl <- map_lgl(l, \(y) same_elem(x, y)) # which version does it match
  stopifnot(sum(lgl) == 1)
  names(l[lgl]) # name of the version
}

region_factor <- function(x = NULL, wafwa_only = FALSE, include_entire = TRUE,
                          return_levels = FALSE,
                          v = NULL) {
  
    

    lev_l <- all_region_levels()
    
    if(return_levels) {
      if(wafwa_only) {
        stopifnot(is.null(v) | v == 'r1.0')
        v = 'wafwa'
      }
      stopifnot(v %in% names(lev_l))
      levels <- lev_l[[v]]
      if(include_entire) {
        levels <- c("Entire study area", levels)
      }
      return(levels)
    }
    # using same_elem to be adjusted if it's ok for x to have levels
    # missing

  x2 <- str_subset(x, 'Entire study area', negate = TRUE) 
  levels <- if(wafwa_only) {
      lev_l$wafwa
    } else if (!is.null(v)) {
      stopifnot(v %in% names(lev_l))
      lev_l[[v]]
    } else if(same_elem(x2, lev_l$r1.0)) {
      lev_l$r1.0
    } else if(same_elem(x2, lev_l$r1.1)) {
      lev_l$r1.1
    } else {
      stop ('correct levels not present')
    }

  if(include_entire) {
    levels <- c("Entire study area", levels)
  }

  stopifnot(x %in% levels)
  factor(x, levels = levels)
}

# convert 'our' regions to the wafwa regions
# for aggregation
region2wafwa <- function(x, v = NULL) {

  x2 <- as.character(x)
  vr <- if(is.null(v)) which_vr(x2) else v
  
  eco <- load_wafwa_ecoregions(total_region = FALSE, v = vr)
  from <- eco$region
  to <- eco$wafwa_NAME
  x2 <- replace_from_to(x, from = from, to = to)

  if(is.factor(x)) {
    include_entire <- any(str_detect(x, 'Entire'))
    x2 <- region_factor(x2, wafwa_only = TRUE, include_entire = include_entire,
                        vr)
  } 
  x2
}



summary2factor <- var2factor_factory( c('low', 'median', 'high'))

driver2factor <-  function(x = NULL, include_none = FALSE,
                           return_levels = FALSE,
                           include_sagebrush = FALSE) {
  levels <- c('MAT', 'MAP', 'PSP', 'Pherb', 'Aherb')
  if(include_sagebrush) levels <- c(levels, 'Sagebrush')
  if(include_none) levels <- c(levels, 'None')
  if(return_levels) return(levels)
  x <- as.character(x)
  x <- str_replace(x, 'psp', 'PSP')
  stopifnot(x %in% levels)
  factor(x, levels)
}

relable_c3_current <- function(x) {
  if(!is.factor(x)) {
    x <- c3_factor(x)
  }
  levels <- levels(x)
  labels <- paste(levels, '(2017-2020)')
  factor(x, levels = levels, 
         labels = labels)
}

relable_graze_long <- function(x) {
  if(!is.factor(x)) {
    x <- graze2factor(x)
  }
  levels <- levels(x)
  labels <- paste(levels, 'grazing')
  factor(x, levels = levels, 
         labels = labels)
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


#' convert character character vector of all PFTs to a factor
#'
#' @param x character vector
#' @return ordered factor 
pft_all_factor <- function(x) {
  levels <-  c("Sagebrush", 
              "C3Pgrass",
              "C4Pgrass",
              "Cheatgrass",
              "Pforb",
              "Shrub", 
              "Pgrass",
              "Pherb",
              "Aforb",
              "Aherb")
  
  # add a 'total' factor level, if the input includes it
  if("Total" %in% x) {
    levels <- c(levels, "Total")
  }
  
  stopifnot(x %in% levels)
  
  out <- factor(x, levels = levels)
  out
}

# convert individual PFTs into a total category, useful for calculating
# total biomass, if return_levels = TRUE, just returns the names of the 
# PFTs that make of total biomass
pft_total_factor <- function(x, return_levels = FALSE) {
  
  pft <- c("sagebrush", "a.cool.forb", "a.warm.forb", "p.cool.forb", 
           "p.warm.forb", "a.cool.grass", "p.cool.grass", "p.warm.grass",
           "shrub", "succulents")
  
  if(return_levels) {
    return(pft)
  }
  
  # x must include all these levels, otherwise won't be computing total biomass
  stopifnot(pft %in% x) 
  
  out <- ifelse(x %in% pft, "Total", NA_character_)
  out <- factor(out, levels = "Total")
  out
}



#' convert 4 letter spp code to PFT
#'
#' @param x character vector that includes 4 letter species codes
#'
#' @return character vector of PFTs, and NAs for non-matched values
spp2pft <- function(x) {
  out <- case_when(
    x == "artr" ~ "sagebrush",
    x == "cryp" ~ "a.cool.forb",
    x == "chen" ~ "a.warm.forb",
    x == "phho" ~ "p.cool.forb", 
    x == "arfr"  ~ "p.warm.forb", 
    x == "brte" ~ "a.cool.grass", 
    x == "pssp"  ~ "p.cool.grass", 
    x == "bogr"  ~ "p.warm.grass",
    x == "chvi" ~ "shrub", 
    x ==  "oppo" ~ "succulents"
  )

  # make sure all pft levels were matched
  all_levels <- pft_total_factor(x = NULL, return_levels = TRUE)
  stopifnot(all_levels %in% out)
  
  out
}


#' Convert columns to factors
#'
#' @param df dataframe with key columns to be converted to factors
#'
#' @return dataframe with key columns converted to factors
df_factor <- function(df, v = NULL) {
  stopifnot(is.data.frame(df))
  
  nms <- names(df)
  
  if("graze" %in% nms) df$graze <- graze2factor(df$graze)
  
  if("years" %in% nms) df$years <- years2factor(df$years)
  
  if("RCP" %in% nms) df$RCP <- rcp2factor(df$RCP)
  
  if("PFT" %in% nms) df$PFT <- pft_all_factor(df$PFT)
  
  if("region" %in% nms) df$region <- region_factor(df$region, v = v)
  if("ecoregion" %in% nms) df$ecoregion <- region_factor(df$ecoregion, v = v)
  
  if('summary' %in% nms) df$summary <- summary2factor(df$summary)
  
  if('c3' %in% nms) df$c3 <- c3_factor(df$c3)
  
  df
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
#' (or if divid_by_max is false, regulare % change)
#' otherwise calculates the actual change (raw difference)
#' @param effect_size logical, whether the response variable calculated
#' is effect size instead of change (or percent change), currently calculated as 
#' ln(trmt/ctrl)
#' @param within_GCM logical, if TRUE comparisons will be made within a GCM, e.g.
#' comparing heavy to light grazing, within a GCM/RCP/years. ref_graze must be 
#' specified if set to TRUE
#' @param divide_by_max should the delta be dived by max current (i.e. to 
#' calculate scaled percent change, if not just calculating normal % change)
#'
#' @return dataframe of percent change in biomass from current conditions,
#' scaled by maximum current biomass
scaled_change <- function(df, 
                          var = "biomass", 
                          by = c("PFT", "graze"),
                          ref_graze = NULL,
                          percent = TRUE,
                          divide_by_max = TRUE,
                          effect_size = FALSE, 
                          within_GCM = FALSE) {
  
  # checking arguments
  stopifnot(
    is.data.frame(df),
    c("RCP", "years", by) %in% names(df)
    )
  
  if(within_GCM | !is.null(ref_graze)){
    stopifnot('graze' %in% names(df))
  }
  
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
             if(percent & divide_by_max) { 
               # % scaled change
               (.data[[var]] - .data$current)/.data$max_value*100
             } else if (percent & !divide_by_max) {
               # regular % change
               (.data[[var]] - .data$current)/.data$current*100
             } else if (effect_size) {
               # log response ratio
               log(.data[[var]]/.data$current)
             }  else {
               # raw change
               .data[[var]] - .data$current
              }
           ) %>% 
    select(-max_value, -current) 
  
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
  names(out) <- names(out) %>% 
    str_replace("^utilization_", "util_") %>% 
    str_replace("^biomass_", "bio_")
  
  out

}

# applying scaled_change() function two or more vars
scaled_change_2var <- function(df, 
                               vars = c("biomass", "indivs"), 
                               by = c("PFT", "graze"),
                               ref_graze = NULL,
                               percent = TRUE,
                               effect_size = FALSE, 
                               divide_by_max = TRUE,
                               within_GCM = FALSE) {
  
  stopifnot(length(vars) >= 2,
            # if not unique names then will have problems
            length(unique(vars)) == length(vars))
  if("biomass" == vars[2]) {
    # this is because biomass variables are renamed
    stop("if one of vars is biomass, it must be the first")
  }
  
  # change of first variable
  out <- scaled_change(df = df,
                       var = vars[1],
                       by = by,
                       ref_graze = ref_graze,
                       percent = percent,
                       effect_size = effect_size,
                       divide_by_max = divide_by_max,
                       within_GCM = within_GCM)
  
  # change of additional variables
  other_vars <- map(vars[2:length(vars)], function(x) {
    tmp <- scaled_change(df = df,
                  var = x,
                  by = by,
                  ref_graze = ref_graze,
                  percent = percent,
                  effect_size = effect_size,
                  divide_by_max = divide_by_max,
                  within_GCM = within_GCM)
    # last col is the diff col
    df_one_col <- tmp[ncol(tmp)]
    df_one_col
  }) %>% 
    list_cbind(name_repair = 'check_unique')
  # scaled change of 2nd variable

  # adding all the diff_var cols together
  out <- bind_cols(out, other_vars, .name_repair = 'check_unique')
  out
}

# helper: proportional change 
prop_change <- function(fut, cur) {
  stopifnot('SpatRaster' %in% class(fut))
  diff <- fut - cur
  prop <- diff/cur
  # want 0/0 to be 0 (i.e. zero change when cur is 0)
  prop <- ifel(diff == 0, 0, prop)
  prop
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

# pivoting ----------------------------------------------------------------

#' Pivot c4 column into wider format
#'
#' @param df dataframe
#'
#' @return dataframe where c4on and c4off columns are created for biomass
#' and number of individuals, the individuals columns are prefixed
#' by "indivs"
pivot_wider_c4 <- function(df) {
  
  stopifnot(c("c4", "indivs", "biomass") %in% names(df))
  
  out <- df %>% 
    pivot_wider(names_from = "c4",
                values_from = c("biomass", "indivs")) %>% 
      # to keep downstream code working biomass c4on and off columns aren't named
      # "biomass_"
      rename(c4on = biomass_c4on,
             c4off = biomass_c4off)
  out
}

# where df is a dataframe results from terra::extract
pivot_longer_extracted <- function(df, regions, values_to) {
  stopifnot(length(regions) == length(df$ID))
  df %>% 
    mutate(region = regions) %>% 
    select(-ID) %>% 
    pivot_longer(-region, names_to = 'id',
                 values_to = values_to)
  
}

# summarize functions --------------------------------------------------

#' summarize across gcms, and then pivot longer
#'
#' @param df grouped dataframe
#' @param var string, name of the variable to summarize
summarize_across_GCMs <- function(df, var) {
  df %>% 
    summarize(median = median(.data[[var]]),
              low = calc_low(.data[[var]]),
              high = calc_high(.data[[var]]), 
              .groups = 'drop') %>% 
    pivot_longer(cols = c('low', 'median', 'high'),
                 names_to = 'summary',
                 values_to = var) %>% 
    mutate(summary = summary2factor(summary))
}

#' calculate median of biomass and indivs columns
#'
#' @param df dataframe with indivs, biomass and bio_suffix, and indivs_suffix
#' columns
#' @param suffix string, the suffix (e.g. _diff, _es), of columns of interest
#' @param include_low_high also provide low/high estimates across GCMs (in addition
#' to the median)
#'
#' @return dataframe, with median values for the biomass and indivs columns
summarise_bio_indivs <- function(df, suffix = "_diff",
                                 col_names = c('biomass', "indivs", "utilization"),
                                 abbreviations = c("bio", "indivs", "util"),
                                 include_low_high = FALSE){
  
  # columns to summarise
  columns <- c(col_names, paste0(abbreviations, suffix))
  

  if(include_low_high) {
    tmp_summarize <- function(df) {
      summarise(df, across(.cols = all_of(!!columns),
                           .fns = list('median' = ~ median(.x, na.rm = TRUE),
                                       'low' = calc_low, 
                                       'high' = calc_high),
                           .names = "{.col}_{.fn}"),
                .groups = "drop") %>% 
        dplyr::rename_with(.fn = ~ str_replace(.x, "_median", ""),
                    .cols = matches("_median$"))
    }
    
    } else {
      tmp_summarize <- function(df) {
        summarise(df, across(.cols = all_of(!!columns),
                             .fns = \(x) median(x, na.rm = TRUE)),
                  .groups = "drop") 
      }
  }
  # using dtplyr because summarise is super slow with across() which is a
  # known issue: https://github.com/tidyverse/dplyr/issues/4953
  
  # note the use of !!columns below, that is a work around a bug in dtplyr
  # https://github.com/tidyverse/dtplyr/issues/164

  out <- df %>% 
    # now using lazy_dt with tmp_summarize isn't working
    # lazy_dt() %>% # i observed up to 30x speed increase by using dtplyr here
    tmp_summarize() %>% 
    as_tibble() # convert back to tibble
  
  out
}

#' add a total herbacious category
#'
#' @param df dataframe
#' @param col_names names of columns to sum across annual PFTs
#' @param group_cols names of columns to group by
#'
#' @return same dataframe as input but with rows added that provide values
#' for all annual herbacious PFTs
calc_aherb <- function(df, col_names = c('biomass', "indivs", "utilization"),
                       group_cols = c('run', 'years', 'RCP', 'graze', 'id', 
                       'site', 'GCM')) {
  
  annuals <- c('Cheatgrass', 'Aforb')
  stopifnot(annuals %in% df$PFT)
  
  if('Aherb' %in% df$PFT) {
    message('Aherb already present')
    return(df)
  }
  aherb <- df %>% 
    filter(.data$PFT %in% !!annuals) %>% 
    group_by(across(all_of(group_cols))) %>% 
    summarise(across(.cols = all_of(col_names),
                     .fns = sum),
              .groups = 'drop') %>% 
    mutate(PFT = 'Aherb')
  
  aherb$PFT <- pft_all_factor(aherb$PFT)
  
  if(!all(names(aherb) %in% names(df))) {
    stop('some col_names or group_cols likely missing')
  }
  stopifnot(levels(aherb$PFT) == levels(df$PFT))
  out <- bind_rows(aherb, df)
  out
}

# classification ----------------------------------------------------------

#' percent correctly classified
#'
#' @param x 2x2 contingency table created with table()
#'
#' @return pcc
pcc <- function(x){
  x <- addmargins(x)
  out <- (x[1, 1] + x[2,2])/x[3,3]*100
  out
}


# filtering ---------------------------------------------------------------
# misc. functions used to filter out rows

#' filter out un-needed climate scenarios, and c4off
#'
#' @param df dataframe
#' @param PFT logical, whether to also filter to only include the main
#'  PFTs interest
#' @param current logical, whether rcp and year columns should inlcude
#' current, should only be set to false in bio-diff variables, where it 
#' is the difference relative to current conditions
#' @param run string defining run parameters
#'
#' @return dataframe only including climate scenarios of interest
#' (the main analysis only includes a subset of scenarios), and focuses 
#' on a given run
filter_rcp_run <- function(df, PFT = FALSE, current = TRUE,
                          run = "fire1_eind1_c4grass1_co20_2503") {
  rcp_levs <- c("Current", "RCP45")
  year_levs <- c("Current", "2070-2100")
  run_string <- run
  
  # return the inputs (for labeling purposes)
  if(is.null(df)) {
    out <- paste(rcp_levs[2], year_levs[2],
                 run, sep = "_")
    return(out)
  }
  if (!current) {
    rcp_levs <- rcp_levs[2]
    year_levs <- year_levs[2]
  }
  
  stopifnot(rcp_levs %in% df$RCP,
            year_levs %in% df$years)
  
  
  out <- df %>% 
    filter(.data$RCP %in% rcp_levs, 
           .data$years %in% year_levs,
           .data$run == run_string) %>%
    droplevels() #drop unused factor levels
  
  # also filter by PFT
  if(PFT) {
    pft_levs <- pft5_factor(x = NULL, return_levels = TRUE)
    
    out <- out[out$PFT %in% pft_levs, ] # filter rows
    out$PFT <- pft5_factor(as.character(out$PFT)) # relevel factor
  }
  out
}

#' filter dataframe by climate scenario etc. 
#'
#' @param df dataframe
#' @param PFT character vector
#' @param run character vector
#' @param years character vector
#'
#' @return
#' dataframe (same as input if character vector inputs left NULL)
filter_scenarios <- function(df, PFT = NULL,
                             run = NULL,
                             years = NULL) {
  stopifnot(
    c("run", "years") %in% names(df)
  )
  
  out <- df %>% 
    lazy_dt() %>% # 4x faster using dtplyr
    filter(is.null(!!run) | .data$run %in% !!run,
           is.null(!!years) | .data$years %in% !!years) 
  if(!is.null(PFT)){
    out <- filter(out, .data$PFT %in% !!PFT)
  }
  
  as_tibble(out)
}

# misc functions ---------------------------------------------------------------

filter_clim_extremes <- function(df, years = '2070-2100') {
  # set which scenarios to examine
  stopifnot(years %in% df$years)
  df %>% 
    filter(RCP == "Current" | (RCP == 'RCP45' & years == !!years) |
             (RCP == 'RCP85' & years == !!years))
}


# joins -------------------------------------------------------------------

#' join stepwat data and subset cell information
#' 
#' @description This is for prepping the input data for the interpolatePoints()
#' function
#'
#' @param step_dat stepwat data (wide format)
#' @param sc_dat subset cell information (including cell numbers and 
#' site number corresponding)
#' @param subset_in_target is the matching function using subset_in_target? If false then no cell numbers are used, instead site numbers are used instead of cell numbers (otherwise the interpolatePoints function doesn't work properly)
#'
#' @return Dataframe with the two df's joined, and cell number for rownames
join_subsetcells <- function(step_dat, sc_dat,
                             subset_in_target = FALSE) {
  
  # join in cell number info
  if(subset_in_target) {
    out <- sc_dat[, c("cellnumber", "site_id")] %>% 
      rename(site = site_id) %>% 
      inner_join(step_dat, by = "site") %>% 
      dplyr::select(-site)
    

  } else {
    # in this case, the unique identifier is 
    # the site, which is used instead of the 
    # "cell number"
    out <- step_dat %>% 
      rename("cellnumber" = "site")
  }

  # rownames needed for interpolatePoints()
  out <- as.data.frame(out) # can't be a tibble b/ of rMultivariateMatching
  rownames(out) <- out$cellnumber
  stopifnot(nrow(out) == 200) # check for join issues
  
  out
}

# temporary files ---------------------------------------------------------

tmp_exists <- function(object_name, rerun, suffix = '.tif') {
  file.exists(paste0('tmp/', object_name, suffix)) & !rerun
}

read_tmp_tif <- function(object_name, dir = 'tmp') {
  rast(paste0(dir, '/', object_name, '.tif'))
}
# for saving intermediate rasters to disk to save memory
writeReadRast <- function(object, objectName, dir = 'tmp') {
  path <- file.path(dir, paste0(objectName, '.tif'))
  terra::writeRaster(object, path, overwrite = TRUE)
  rm(list = objectName, envir = .GlobalEnv)
  gc()
  terra::rast(path)
}

# misc. -------------------------------------------------------------------

# number of unique elements in a vector
lu <- function(x) length(unique(x))



#' Grazing level that causes a site's biomass to drop below threshold
#' 
#' @description Function should be applied in a summarise statement, where
#' the datafarme is grouped by scenario, PFT and site, but not by grazing
#' level. 
#'
#' @param graze factor of grazing level (length 4)
#' @param above logical, whether the given grazing level has biomas above 
#' threshold
#'
#' @return integer with values of 1 to 5, 1 means low grazing biomass
#' was below threshold, 2, means that mildes grazing level that was below
#' threshold was moderate, 5 means that no grazing level was below that threshold
cross_threshold <- function(graze, above) {
  
  stopifnot(
    length(graze) == 4,
    length(above) == 4,
    is.logical(above),
    lu(graze) == 4, # there should be 4 grazing levels
    is.factor(graze) # input needs to be an ordered factor
  )
  # grazing level (integer) (works only if graze is an ordered factor)
  graze_int <- as.integer(graze)
  
  min_graze <- if(all(above)) {
    # when no grazing level caused biomass to go below threshold 
    # (there are 4 grazing levels)
    5 
  } else {
    # the minimum (mildest), grazing level that caused biomass to go
    # below the threshold, 
    min(graze_int[!above])
  }
  
  min_graze
}



#' Percent area crossing threshold
#'
#' @param lyr name of the layer of the SpatRaster
#' @param r spatRaster object
#'
#' @return dataframe, The percent of total area belonging to each of the threshold
#' categories. This is actual area (accounting for non-equal area projections)
calc_pcent_by_thresh <- function(lyr, r) {
  
  r <- r[[lyr]]
  
  # calculating total area of all cells falling
  # into the various threshold categories
  
  area <- zonal(cellSize(r), r, fun = 'sum')
  x <- area$area/sum(area)*100
  levels <- c("Light grazing", 
              "Moderate grazing", 
              "Heavy grazing", 
              "Very Heavy Grazing",
              "Didn't cross threshold")
  
  # the raster doesn't necessarily have all 5 values
  # so just extracting the values it does have
  names(x) <- levels[area[[1]]]
  
  out <- as.data.frame(as.list(x))
  out
}


#' create dataframe of raster information
#'
#' @param rast raster (with named layers), or vector of file names
#' @param run regular expression the defines the run, the 
#' "(_\\d{4}){0,1}" part of the default regex is meant to optionally match
#' the two digit year and month code that started to be appendix to newer run ids
#' which is sometimes also folled by a version (e.g. v3)
#' @param into columns the components of the layer names will be seperated
#' into
#' @param id_noGCM logical, create a column that provides the layer id but
#' with the GCM removed (for grouping/summarizing across GCMs)
#' @return dataframe\
#' @example
#' path <- "grazing_effects/data_processed/interpolated_rasters/biomass/fire1_eind1_c4grass1_co20_Aforb_biomass_RCP45_2030-2060_Light_inmcm4.tif"
#' create_rast_info(path)
create_rast_info <- function(x,
                             run_regex = "fire\\d_eind\\d_c4grass\\d_co2\\dv{0,1}\\d{0,1}(_\\d{4}v{0,1}\\d{0,1}){0,1}",
                             into = c("PFT", "type", "RCP", "years", 
                                      "graze", "GCM"),
                             id_noGCM = FALSE) {
  
  if(isTRUE('SpatRaster' %in% class(x))) {
    y <- names(x)
  } else if (is.character(x)) {
    y <- stringr::str_replace(basename(x), ".tif$", "")
  } else {
    stop("incorrect input type (spatraster or vector needed)")
  }
  
  lookup_graze <- c("Light" = "grazL",
                    "Moderate" = "grazM",
                    "Heavy" = "grazH",
                    "VeryHeavy" = "grazVH")
  run <- stringr::str_extract(y, run_regex)  
  out <- tibble(id = y) %>% 
    mutate(run = str_extract(id, run_regex),
           to_split = str_replace(id, paste0(run_regex, "_"), "")) %>% 
    separate(col = 'to_split',
             into = into,
             sep = "_")
  if('graze' %in% names(out)) {
    out$run2 = paste(run, lookup_graze[out$graze], sep = "_")
  }
  out <- out %>%  
    df_factor() %>% 
    dplyr::select(run, all_of(into), everything())
  
  if(id_noGCM) {
    out <- out %>% 
      mutate(
        id_noGCM = str_replace(id, "_[^_]*$", ""),
        # remove the year-month from run name,
        # this is b/ different year-month runs may have been
        # done for different grazing levels, but still represent all
        # the same other run parameters
        id_noGCM = str_replace(id_noGCM, '(?<=co2[01])_\\d{4}', ''))
    
    check <- out %>% 
      group_by(id_noGCM) %>% 
      summarize(n = n())
    
    # should only be 1 (for current conditions), or 13 (for futures)
    if(!all(check$n %in% c(1, 13))) {
      stop('wrong number of id_noGCMs replicates')
    }
  }
  
  if(!all(complete.cases(out))) {
    stop('some parsings failed (NAs created)')
  }
  out
}


#' calculating matching quality
#' 
#' @description
#' for calculating matching quality for different matching criteria
#' then used in multivarmatch()
#' 
#'
#' @param match output of of multivarmatch, that contain target_cell and subset_cell columns
#' @param tc contains bioclim_vars columuns for each cell number
#' @param bioclim_vars vector of names of the bioclim_vars
#' @param crit matching criteria
#'
#' @return dataframe with cellnumber and matching quality
matchqual <- function(match, tc, bioclim_vars, crit) {
  tc2 <- tc
  
  match <- match1 %>% 
    mutate(target_cell = as.numeric(target_cell),
           subset_cell = as.numeric(subset_cell))
  
  target <- match %>% 
    left_join(tc2, by = c('target_cell' = 'cellnumber')) %>% 
    rename(cellnumber = target_cell) %>% 
    dplyr::select(cellnumber, all_of(bioclim_vars))
  
  subset <- match %>% 
    left_join(tc2, by = c('subset_cell' = 'cellnumber')) %>% 
    rename(cellnumber = subset_cell) %>% 
    dplyr::select(cellnumber, all_of(bioclim_vars))
  
  diff_norm <- target
  
  for(var in bioclim_vars) {
    # normalized squared difference
    diff_norm[[var]] <- ((target[[var]] - subset[[var]])/crit[var])^2
  }
  # euclidean distance
  diff_norm$matching_quality <- sqrt(rowSums(diff_norm[bioclim_vars]))
  
  diff_norm[, c('cellnumber', 'matching_quality')]
}

#' convert layer to depth
#'
#' @param x numeric vector (soil layer in SOILWAT2)
#'
#' @return numeric vector of depth (cm) (midpoint of layer)
#' @export
#'
#' @examples
#' lyr2depth(c(1, 1, 2:8))
lyr2depth <- function(x) {
  stopifnot(x %in% 1:8)
  # bottom of layers
  #depth_lookup <- c(10, 20, 30, 40, 60, 80, 100, 150)
  # midpoint of layer:
  depth_lookup <- c(5,  15,  25,  35,  50,  70,  90, 125)
  out <- depth_lookup[x]
  out
}



#' convert tranpsiration to a rate per cm soil (for
#' to compare unequal width soil bins)
#'
#' @param x transpiration
#' @param lyr soil layer, or soil layer category
transp_per_cm <- function(x, lyr) {
  if(is.numeric(lyr) & all(lyr %in% 1:8)) {
    lookup <- c(10, 10, 10, 10, 20, 20, 20, 50)
    width <- lookup[lyr]
  } else if(all(as.character(lyr) %in% c("0-10 cm", "10-40 cm", "40-150 cm"))){
    lookup <- c("0-10 cm" = 10 , "10-40 cm" = 30 , "40-150 cm" = 110)
    width <- lookup[as.character(lyr)]
  } else {
    stop('values in lyr not recognized')
  }
  x/width # transpiration per cm of soil depth
}

#' calculates the mode and if there are multiple, randomly returns one of them
# (by default)
#'
#' @param x vector (numeric or character)
#' @param random randomly return just one of the modes
mode_rand <- function(x, random = TRUE) {
  # adapted from: https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  out <- ux[tab == max(tab)]
  n <- length(out)
  if(random & n > 1) {
    i <- sample(1:n, size = 1)
    out <- out[i]
  }
  out
}

# multiplies the number of pixels
# by the area per pixel--at the moment all are 1km2 pixels (100ha)
# (using equal area projection)
pixel2area <- function(x) {
  x*100
}




weight2area <- function(x) {
  x*100 # each pixel has an area of 100 ha
}

# extracts the first letter each word, to create abbreviations
# keeps existing abbrevionts (WY Basin becomes 'WYB')
words2abbrev <- function(x) {
  str_split(x, "\\s+") %>%
    map_chr(~ {
      abbrev <- map_chr(.x, function(word) {
        if (str_detect(word, "^[A-Z]{2,}$")) word else str_to_upper(str_sub(word, 1, 1))
      })
      str_c(abbrev, collapse = "")
    })
}

same_elem <- function(x, y) {
  all(x %in% y) & all(y %in% x)
}

replace_from_to <- function(x, from, to) {
  out <- to[match(x, from)]
  out[is.na(out)] <- x[is.na(out)]
  out
}

names_replace <- function(x, pattern, replacement) {
  nms <- names(x)
  names(x) <- stringr::str_replace(nms, pattern, replacement)
  x
}
