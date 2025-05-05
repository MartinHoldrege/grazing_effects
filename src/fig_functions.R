# Martin Holdrege

# Purpose: functions used in the creation of figures. Putting functions
# here that might be useful across a couple of different scripts going forward

# a number of the functions here require the 'src/fig_params.R' script
# to also be sourced. 


# params used by some functions below ----------------------------------------

# these params are defined here instead of in fig_params.R b/ they are
# only used in the functions below. 

# location of legends in boxplots
legend_pos_box1 <- "top"

line_loc <- c(5.5, 10.5, 15.5) # locations to draw vertical lines on boxplot
outlier.size <-  0.5 # outlier size used in boxplot functions below

# ggplot themes -----------------------------------------------------------

theme_custom1 <- function() {
  theme_bw() %+replace%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.background = element_blank())
}

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

# base of absolute biomass boxplots 
box1 <- function(var = "biomass", y = lab_bio0, add_facet_wrap = TRUE,
                 group_by = c("PFT", "graze")) {
  # first plotting text, so it doesn't overplot data
  out <- list(
    geom_text(data = ~box_anno(., var = var, group_by = group_by,
                               mult = 0.05),
              aes(x, y, label = graze, fill = NULL), 
              size = 2.5),
    geom_boxplot(outlier.size = outlier.size), # not showing outliers as points
    scale_fill_manual(values = cols_rcp, name = "Scenario"),
    scale_x_discrete(labels = id2year),
    geom_vline(xintercept = line_loc, linetype = 2),
    theme(legend.position = legend_pos_box1,
          axis.text = element_text(size = 7)),
    labs(x = lab_yrs,
         y = y)
  )
  
  if(add_facet_wrap) {
    # ncol_box currently defined in the script doing the sourcing
    out[["wrap"]] <- facet_rep_wrap(~ PFT, scales = "free", ncol = ncol_box)
  }
  out
}

# boxplot of change in biomass (scaled percent and effect size), 
# for each of the 5 main
# PFTs, by, RCP, grazing treatment and time period, relative to
# current conditions for that grazing level
box2 <- function(axis_data, # axis data can be a different data frame than
                 # is plotted, and that way axes can be constant across 
                 # multiple figures
                 var = "bio_diff", # y variable
                 mult = 0.05,
                 box_identity = FALSE, # should stat = identity be used to make 
                 # the figure
                 subtitle = "Comparing to current conditions within a grazing level",
                 xintercept = 2.5,
                 repeat.tick.labels = FALSE,
                 # whether box_anno annotation should be at the same height for
                 # all panels (useful if scales are fixed)
                 anno_same_across_panels = FALSE, 
                 scales = "free_y",
                 ...){
  
  
  out <- list(
    geom_hline(yintercept = 0, alpha = 0.3, linetype = 1),
    facet_rep_wrap(~PFT, scales = scales, ncol = ncol_box,
                   repeat.tick.labels = repeat.tick.labels),
    scale_fill_graze(),
    # so just display the RCP
    scale_x_discrete(labels = id2year),
    theme(legend.position = legend_pos_box1),
    geom_vline(xintercept = xintercept, linetype = 2),
    labs(x = lab_yrs, subtitle = subtitle),
    # text and empty points based on a different dataframe, so that
    # axis limits among multiple figures can be the same
    geom_text(data = box_anno(axis_data, var = var, 
                              group_by = c("PFT", "RCP"),
                              id = "id2", mult = mult,
                              anno_same_across_panels = anno_same_across_panels),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5),
    geom_point(data = axis_data, aes_string(y = var), color = NA)
  )
  
  # create the boxplot using stat 'identity instead
  if(box_identity) {
    out$box <- geom_boxplot_identity(position = "dodge")
    # identity = boxplot (i.e. regular data input)
  } else {
    out$box <- geom_boxplot(position = "dodge",
                            outlier.size = outlier.size,
                            ...)
  }
  out
}

# boxplots showing change in biomass (and effect size) relative to current time 
# period and a given grazing intensity ('reference class')
box3 <- function(axis_data, var = "bio_diff") {
  list(
    geom_hline(yintercept = 0, alpha = 0.3, linetype = 1),
    geom_boxplot(position = position_dodge(preserve = "single"),
                 outlier.size = outlier.size),
    facet_rep_wrap(~PFT, scales = "free", ncol = ncol_box,
                   repeat.tick.labels = "y"),
    scale_fill_graze(),
    # so just display the RCP
    scale_x_discrete(labels = id2year),
    theme(legend.position = legend_pos_box1),
    geom_vline(xintercept = line_loc2, linetype = 2),
    labs(x = lab_yrs),
    geom_text(data = ~box_anno(axis_data, var = var, 
                               group_by = c("PFT", "RCP"),
                               id = "id2"),
              aes(x, y, label = RCP, fill = NULL),
              size = 2.5),
    geom_point(data = axis_data, color = NA)
  )
}

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

weighted_box1 <- function(df, y_string, ylab = NULL, subtitle = NULL) {
  ggplot(df, aes(graze, .data[[y_string]], fill = summary)) +
    geom_boxplot(aes(weight = weight), position = position_dodge2(preserve = 'single'),
                 outlier.size = 0.25, outlier.alpha = 0.5) + 
    facet_grid(region~rcp_year) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x.top = element_text(size = rel(0.7))) +
    labs(x = lab_graze,
         y = ylab,
         subtitle = subtitle,
         fill = 'Summary across GCMs')
}

weighted_violin1 <- function(df, y_string, ylab = NULL, subtitle = NULL) {
  ggplot(df, aes(graze, .data[[y_string]], fill = summary)) +
    geom_violin(aes(weight = weight), position = position_dodge2(preserve = 'single'),
                 outlier.size = 0.25, outlier.alpha = 0.5) + 
    facet_grid(region~rcp_year) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x.top = element_text(size = rel(0.7))) +
    labs(x = lab_graze,
         y = ylab,
         subtitle = subtitle,
         fill = 'Summary across GCMs')
}


# modify/combine groups of plots ------------------------------------------



# combine side by side groups of facet_grid plots
combine_grid_panels1 <- function(l, remove_y_axis_labels = TRUE) {
  n <- length(l)
  stopifnot(length(l) > 1)
  # remove duplicated y axis labels (b/ putting figures side by side)
  for (i in 2:n) {
    l[[i]] <- l[[i]] + theme(axis.title.y = element_blank())
  }
  
  # only the right most figure will have strip text (to avoid repeat's)
  for (i in 1:(n-1)) {
    l[[i]] <- l[[i]] + theme(strip.text.y = element_blank())
  }
  
  if(remove_y_axis_labels) {
    # this is useful if y axis limits are identical across panels
    for (i in 2:n) {
      l[[i]] <- l[[i]] + theme(axis.text.y = element_blank())
    }
  }
  g <- wrap_plots(l, nrow = 1) +
    plot_layout(guides = 'collect', axis_titles = 'collect')
  g2 <- g&theme(legend.position = 'bottom')
  g2
  
}

#' remove axis labels from list of ggplots
#'
#' @param plot_list list of ggplot object
#' @param index_keep_y indices (of the list) of ggplot objects to keep the y axis title for
#' @param index_keep_y_sec indices to keep the secondary y axis for
remove_y_titles <- function(plot_list, 
                            index_keep_y = NULL, 
                            index_keep_y_sec = NULL){
  ind <- 1:length(plot_list)
  stopifnot(is.list(plot_list))
  if(!is.null(index_keep_y)) {
    stopifnot(index_keep_y <= max(ind) & index_keep_y >= 1)
    ind_y <- ind[ind != index_keep_y]
  } else {
    ind_y <- ind
  }
  plot_list[ind_y] <- purrr::map(plot_list[ind_y], function(x) {
    x + labs(y = NULL)
  })
  
  if(!is.null(index_keep_y_sec)) {
    stopifnot(index_keep_y_sec <= max(ind) & index_keep_y_sec >= 1)
    ind_sec <- ind[ind != index_keep_y_sec]
  } else {
    ind_sec <- ind
  }
  plot_list[ind_sec] <- purrr::map(plot_list[ind_sec], function(x) {
    x + theme(axis.title.y.right = element_blank())
  })
  
  plot_list
}

# funs for scatterplots ---------------------------------------------------

# given a ggplot (g) create scatter plots y ~ MAP and y ~ MAT
climate_scatter <- function(g, include_psp = FALSE) {
  out <- list()
  # vs MAP
  out[["MAP"]] <- g +
    geom_point(aes(x = PPT)) +
    geom_smooth(aes(x = PPT), method = "loess", se = FALSE) +
    labs(x = lab_map)
  
  # vs MAT
  out[["MAT"]] <- g +
    geom_point(aes(x = Temp)) +
    geom_smooth(aes(x = Temp), method = "loess", se = FALSE) +
    labs(x = lab_mat)
  
  # average correlation between monthly precip and temp (seasonality metric)
  out[["CorrTP2"]] <- g +
    geom_point(aes(x = CorrTP2)) +
    geom_smooth(aes(x = CorrTP2), method = "loess", se = FALSE) +
    labs(x = lab_corrtp)
  
  if(include_psp) {
    out[["psp"]] <- g +
      geom_point(aes(x = psp)) +
      geom_smooth(aes(x = psp), method = "loess", se = FALSE) +
      labs(x = lab_psp)
  }

  out
}

# base of scatter plots for change in light grazing relative to light grazing
# vs climate
scatter_light <- function(pft, # for subtitle
                          run, # c4on/off for subtitle
                          axis_data, # for setting axis limits
                          subtitle_response = "biomass") {
  list(
    geom_hline(yintercept = 0, linetype = 2),
    facet_rep_wrap(~RCP + years),
    scale_color_graze(),
    labs(caption = paste(run,
                         "\nReference class is light grazing under current",
                         "conditions"),
         subtitle = paste("Change in", pft, subtitle_response)),
    theme(legend.position = c(0.85, 0.15),
          axis.text = element_text(size = 7)),
    # for setting axis limits based on bigger data frame, so multiple
    # figs can have same y axis lims
    geom_point(data = axis_data,
               x = NA)
  )
}
# axis functions ----------------------------------------------------------

# from 'id' column get the label to use for time period in boxplots

id2year <- function(x) {
  # x--id (character vector), that includes the time period in it
  # returns--character vector of letter designations of the time period
  out <- case_when(
    str_detect(x, "Current") ~ "H",
    str_detect(x, "_2030-2060") ~ "M",
    str_detect(x, "_2070-2100") ~ "E"
  )
  if(any(is.na(out))) {
    warning("Not all time periods matched")
  }
  out
}

# return the grazing level from an id string
id2graze <- function(x, abbrev = FALSE) {
  # x--id (character vector), that includes the grazing level in it
  # returns--character vector of the grazing level
  out <- case_when(
    str_detect(x, "Light") ~ "Light",
    str_detect(x, "Moderate") ~ "Moderate",
    str_detect(x, "Very[ ]{0,1}Heavy") ~ "Very Heavy",
    str_detect(x, "Heavy") ~ "Heavy"
    
  )
  if(any(is.na(out))) {
    warning("Not all time years matched")
  }
  
  # return just first letters
  if(abbrev) {
    out <- purrr::map_chr(str_extract_all(out, "\\b\\w"), \(x) paste(x, collapse = ""))
  }
  out
}

# for adding a second y axis to effect size figures, that shows % change
add_sec_axis <- function(name = "% Change", breaks = waiver(),...) {
  scale_y_continuous(sec.axis = sec_axis(trans = es2pchange, 
                                         name = name,
                                         ...)) 
}

# label functions ---------------------------------------------------------

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

#' paste rcp and years together
#'
#' @param rcp character vector of RCPs
#' @param years character vector of year ranges
#' @param add_letters add letters (for figure labeling)
#' @param include_parenth put year in parantheses
#'
#' @return character vector
rcp_label <- function(rcp, years, add_letters = FALSE,
                      include_parenth = TRUE, add_newline = FALSE) {
  years <- as.character(epoch2factor(years)) # convert 2070 to 2071 etc. so labels are correct
  rcp <- update_rcp(rcp)
  if(add_newline) {
    space <- '\n'
  } else {
    space <- " "
  }
  if (include_parenth) {
    x1 <- ifelse(rcp == "Current" | rcp == "Historical", 
                 "(Historical)", paste0("(",rcp,",", space, years, ")"))
  } else {
    x1 <- ifelse(rcp == "Current" | rcp == "Historical", 
                 "Historical", paste0(rcp, space, years))
  }
  
  
  
  if(add_letters) {
    levs <- sort(unique(x1))
    labels <- paste(fig_letters[1:length(levs)], levs)
    x2 <- factor(x1, levels = levs, labels = labels)
  } else {
    x2 <- x1
  }
  x2
}



# color functions ---------------------------------------------------------

# grazing colors (cols_graze defined in the fig_params.R script)
scale_color_graze <- function() {
  scale_color_manual(values = cols_graze, name = "Grazing")
}

scale_fill_graze <- function(exclude = NULL) {
  col <- cols_graze
  if(!is.null(exclude)) {
    stopifnot(exclude %in% names(col))
    col <- col[exclude != names(col)] 
  }

  scale_fill_manual(values = col, name = "Grazing")
}

#' create values vector for ggplot, with midpoint, not 
#' at mathematical midpoint, useful for not centered diverging color pallets
#'
#' @param colors vector of colors
#' @param midpoint 
#' @param limits vector contain the min and max
#'
#' @returns
#' vector 
#' 
#' @examples
#' colors <- c("#bd0026", "#e31a1c", "#fc4e2a", "#feb24c", "#fed976", "#ffeda0", 
#' "grey", "grey", "#deebf7", "#c6dbef", "#9ecae1", "#4292c6", "#08519c", 
#' "#08306b")
#' values_about_midpoint(colors, midpoint = 0, limits = c(-1, 4))
values_about_midpoint <- function(colors, midpoint, limits) {
  stopifnot(midpoint > limits[1] & midpoint < limits[2],
            length(limits) == 2,
            length(colors) >= 3
  ) 
  
  n <- length(colors)
  
  # odd number of colors, put middle color at the midpoint
  if(n %% 2 == 1) {
    seq1 <- seq(from = limits[1], to = midpoint, length.out = (n+1)/2)
    seq2 <- seq(from = midpoint, to = limits[2], length.out = (n+1)/2)
    values <- scales::rescale(c(seq1, seq2[-1]))
  } else {
    # even number, putting the middle two colors centered around the midpoint
    seq1 <- seq(from = limits[1], to = midpoint, length.out = n/2)
    seq2 <- seq(from = midpoint, to = limits[2], length.out = n/2)
    
    delta <- min(c(diff(seq1), diff(seq2)))/2 # half a the min distance between numbers
    seq1 <- seq1 - delta # moving numbers down 'half' a unit
    seq2 <- seq2 + delta # moving numbers up 'half' a unit
    values <- scales::rescale(c(seq1, seq2))
  }
  
  values
}
  


# write functions ---------------------------------------------------------

png20panel <- function(filename) {
  png(filename, res = 600, height = 5*1.8 + 1, width = 4*1.8 - 0.3, units = 'in')
}
