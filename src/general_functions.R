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
              "Aforb")
  
  # add a 'total' factor level, if the input includes it
  if("Total" %in% x) {
    levels <- c(levels, "Total")
  }
  
  stopifnot(x %in% levels)
  
  out <- factor(x, levels = levels)
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
df_factor <- function(df) {
  stopifnot(is.data.frame(df))
  
  nms <- names(df)
  
  if("graze" %in% nms) df$graze <- graze2factor(df$graze)
  
  if("years" %in% nms) df$years <- years2factor(df$years)
  
  if("RCP" %in% nms) df$RCP <- rcp2factor(df$RCP)
  
  if("PFT" %in% nms) df$PFT <- pft_all_factor(df$PFT)
  
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
  if("biomass_diff" %in% names(out)) {
    out <- out %>% 
      rename(bio_diff = biomass_diff)
  } else if ("biomass_es" %in% names(out)) {
    out <- out %>% 
      rename(bio_es = biomass_es)
  }
  
  out

}

scaled_change_2var <- function(df, 
                               vars = c("biomass", "indivs"), 
                               by = c("PFT", "graze"),
                               ref_graze = NULL,
                               percent = TRUE,
                               effect_size = FALSE, 
                               within_GCM = FALSE) {
  
  stopifnot(length(vars) ==2)
  if("biomass" == vars[2]) {
    # this is because biomass variables are renamed
    stop("if one of vars is biomass, it must be the first")
  }
  
  # scaled change of first variable
  out <- scaled_change(df = df,
                       var = vars[1],
                       by = by,
                       ref_graze = ref_graze,
                       percent = percent,
                       effect_size = effect_size,
                       within_GCM = within_GCM)
  
  # scaled change of 2nd variable
  df_var2 <- scaled_change(df = df,
                           var = vars[2],
                           by = by,
                           ref_graze = ref_graze,
                           percent = percent,
                           effect_size = effect_size,
                           within_GCM = within_GCM)
  names <- names(df_var2)
  diff_var2 <- names[length(names)] # name of the last column
  
  # make sure will not be overwriting an existing column
  stopifnot(!(diff_var2) %in% names(out))
  
  # adding 2nd diff var column, to first dataframe
  out[[diff_var2]] <- df_var2[[diff_var2]]
  
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



# summarize functions -----------------------------------------------------

#' calculate median of biomass and indivs columns
#'
#' @param df dataframe with indivs, biomass and bio_suffix, and indivs_suffix
#' columns
#' @param suffix string, the suffix (e.g. _diff, _es), of columns of interest
#'
#' @return dataframe, with median values for the biomass and indivs columns
summarise_bio_indivs <- function(df, suffix = "_diff"){
  
  # columns to summarise
  columns <- c("biomass", "indivs", paste0("bio", suffix), 
            paste0("indivs", suffix))
  
  # using dtplyr because summarise is super slow with across() which is a
  # known issue: https://github.com/tidyverse/dplyr/issues/4953
  
  # note the use of !!columns below, that is a work around a bug in dtplyr
  # https://github.com/tidyverse/dtplyr/issues/164
  out <- df %>% 
    lazy_dt() %>% # i observed up to 30x speed increase by using dtplyr here
    summarise(across(.cols = all_of(!!columns),
                     .fns = median, na.rm = TRUE),
              .groups = "drop") %>% 
    as_tibble() # convert back to tibble
  
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
#' @param c4 string, c4on or c4off
#'
#' @return dataframe only including climate scenarios of interest
#' (the main analysis only includes a subset of scenarios), and focuses 
#' on c4on simulations (by default)
filter_rcp_c4 <- function(df, PFT = FALSE, current = TRUE,
                          c4string = "c4on") {
  rcp_levs <- c("Current", "RCP85")
  year_levs <- c("Current", "2030-2060")
  
  if (current) {
    rcp_levs <- c("Current", "RCP85")
    year_levs <- c("Current", "2030-2060")
  } else {
    rcp_levs <- "RCP85"
    year_levs <- "2030-2060"
  }
  
  stopifnot(rcp_levs %in% df$RCP,
            year_levs %in% df$years)
  out <- df %>% 
    filter(.data$RCP %in% rcp_levs, 
           .data$years %in% year_levs,
           .data$c4 == c4string) %>%
    droplevels() #drop unused factor levels
  
  # also filter by PFT
  if(PFT) {
    pft_levs <- pft5_factor(x = NULL, return_levels = TRUE)
    
    out <- out[out$PFT %in% pft_levs, ] # filter rows
    out$PFT <- pft5_factor(as.character(out$PFT)) # relevel factor
  }
  out
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
#'
#' @return Dataframe with the two df's joined, and cell number for rownames
join_subsetcells <- function(step_dat, sc_dat) {
  
  # join in cell number info
  out <- sc_dat[, c("cellnumber", "site_id")] %>% 
    rename(site = site_id) %>% 
    inner_join(step_dat, by = "site") %>% 
    dplyr::select(-site)
  
  # rownames needed for interpolatePoints()
  rownames(out) <- out$cellnumber
  
  stopifnot(nrow(out) == 200) # check for join issues
  
  out
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
#' @param run regular expression the defines the run
#' @param into columns the components of the layer names will be seperated
#' into
#'
#' @return dataframe\
#' @example
#' path <- "grazing_effects/data_processed/interpolated_rasters/biomass/fire1_eind1_c4grass1_co20_Aforb_biomass_RCP45_2030-2060_Light_inmcm4.tif"
#' create_rast_info(path)
create_rast_info <- function(x,
                             run_regex = "fire\\d_eind\\d_c4grass\\d_co2\\d",
                             into = c("PFT", "type", "RCP", "years", 
                                      "graze", "GCM")) {
  
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
             sep = "_") %>% 
    mutate(run2 = paste(run, lookup_graze[graze], sep = "_")) %>% 
    df_factor() %>% 
    select(run, all_of(into), everything())
  
  if(!all(complete.cases(out))) {
    stop('some parsings failed (NAs created)')
  }
  out
}
