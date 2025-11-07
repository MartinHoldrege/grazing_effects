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


# dependencies ------------------------------------------------------------

source("src/fig_params.R")

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

weighted_box1 <- function(df, y_string, ylab = NULL, subtitle = NULL,
                          strip.text.y = element_text(size = rel(0.7))) {
  ggplot(df, aes(graze, .data[[y_string]], fill = summary)) +
    geom_boxplot(aes(weight = weight), position = position_dodge2(preserve = 'single'),
                 outlier.size = 0.25, outlier.alpha = 0.5) + 
    facet_grid(region~rcp_year) +
    theme(panel.spacing.x = unit(0, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x.top = element_text(size = rel(0.7)),
          strip.text.y = strip.text.y) +
    scale_fill_smry() +
    labs(x = lab_graze,
         y = ylab,
         subtitle = subtitle)
}

# boxplot (same as weighted box, but w/ fivenumber summary)
box_fivenum1 <- function(df, ylab = NULL, subtitle = NULL) {
  ggplot(df, aes(graze, ymin = min, lower = lower, middle = middle, upper = upper,
                 ymax = max, fill = summary)) +
    geom_boxplot(position = position_dodge2(preserve = 'single'),
                 stat = 'identity',
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





#' side-by side boxplots of absolute and change values by PFT, by grazing intensity
#'
#' @param df_abs dataframe of historical conditions
#' @param y_abs name of y variable in df_abs
#' @param y_diff y variable in df_diff dataframse
#' @param df_diff_cref difference data frame (climate reference)
#' @param df_diff_gref difference data frame (grazing reference)
#' @param ylab_diff_cref 
#' @param ylab_diff_gref 
#' @param ylab_abs label for abs panels
box_abs_diff <- function(df_abs, 
                         df_diff_cref,
                         df_diff_gref,
                         y_abs = 'biomass',
                         y_diff = 'bio_diff',
                         ylab_abs = y_abs,
                         ylab_diff_cref = y_diff,
                         ylab_diff_gref = y_diff,
                         scales_cref = 'free_y',
                         scales_gref = 'fixed'
) {
  
  stopifnot(unique(sort(df_abs$PFT)) == unique(sort(df_diff_gref$PFT)))

  base <- function(add_hline = FALSE, scales = 'free_y') {
    out <- list(
      geom_boxplot(aes(weight = weight), coef = 10,
                   position = position_dodge2(preserve = 'single')),
      scale_fill_smry(),
      facet_grid(PFT~rcp_year, scales = scales),
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.spacing.x = unit(0, "lines"),
            plot.subtitle = element_text(hjust = 0.5))

    )
    if(add_hline) {
      out <- c(out, list(geom_hline(yintercept = 0, alpha = 0.5, linetype = 2)))
    }
    out
  }
  
  blank_strip <- function() {
    theme(strip.background.y = element_blank(), strip.text.y = element_blank())
  }

  g1 <- ggplot(df_abs, aes(graze, .data[[y_abs]], fill = summary)) +
    base() +
    blank_strip() +
    labs(x = lab_graze,
         y = ylab_abs) +
    # add space on the left (for tags)
    # (and adding space to the right, to be symetric)) 
    #scale_x_discrete(expand = expansion(add = c(1, 1))) +
    guides(fill = 'none') +
    expand_limits(y = 0)
  

  g2  <- ggplot(df_diff_gref, aes(graze, .data[[y_diff]], fill = summary))+
    base(add_hline = TRUE, scales = scales_gref) +
    blank_strip() +
    labs(x = lab_graze,
         y = ylab_diff_gref,
         subtitle = 'Grazing effect') 

  g3  <- ggplot(df_diff_cref, aes(graze, .data[[y_diff]], fill = summary))+
    base(add_hline = TRUE, scales = scales_cref) +
    labs(x = lab_graze,
         y = ylab_diff_cref,
         subtitle = 'Climate effect')  
  
  # letters for each panel (running rowwise left to right)
  ncol1 <- lu(df_abs$rcp_year)
  ncol2 <- lu(df_diff_gref$rcp_year)
  ncol3 <- lu(df_diff_cref$rcp_year)
  nrow <- lu(df_abs$PFT)
  ncol <- ncol1 + ncol2 + ncol3

  m <- matrix(data = fig_letters[1:(ncol*nrow)], ncol = ncol, byrow = TRUE)
  let1 <- sort(as.vector(m[, 1:ncol1]))
  let2 <- sort(as.vector(m[, (ncol1 + 1):(ncol1 + ncol2)]))
  let3 <- sort(as.vector(m[, (ncol1 + ncol2 + 1):ncol(m)]))
  
  # adding letters to panels

  g1b <- tag(g1, let1)
  g2b <- tag(g2, let2)
  g3b <- tag(g3, let3)

  comb <- g1b + g2b + g3b + plot_layout(widths = c(ncol1, ncol2, ncol3),
                          guides = 'collect')
  comb&theme(legend.position = "bottom")
}

# *stacked bar ------------------------------------------------------------

# Stacked bars (% left axis, area right axis), dodged by grazing scenario.
# One row per ecoregion; within the row, two panels: Current vs chosen RCP.
# Each panel gets its own right-side axis transform using its area_region.


#' @title Build a stacked bar panel of SEI classes
#' @description Creates a stacked bar chart for one region × climate facet (e.g. Current or RCP scenario),
#' with left axis showing % of region and right axis showing area in ha.
#' @param df_panel Data frame filtered to one region and one climate facet. Must contain:
#'   \code{graze}, \code{c3_percent} (0–100), \code{c3} (SEI class), and \code{area_region}.
#' @param area_total Numeric, total area of the region in ha.
#' @param legend.position.inside Numeric vector of length 2 giving legend position inside plot.
stacked_panel <- function(df_panel,  area_total = NULL,
                          legend.position.inside = c(0.08, 0.2)) {
  # df_panel: filtered to one region & one climate facet (Current or future)
  # expects columns: graze, c3_percent, c3, area_region
  df_panel$c3 <- forcats::fct_rev(df_panel$c3) # so CSA at bottom of stack
  
  scale_y <- if(!is.null(area_total)) {
      function() {
        scale_y_continuous(
          expand = expansion(mult = c(0, 0.02)),
          sec.axis = sec_axis(~ . * area_total / 100,
                              name = NULL,
                              labels = \(x) x/1e6 # convert to millions ha
          )  
        )
      }
    } else {
     function() {
        scale_y_continuous(
        )
    }
    }
  x_text <- (lu(df_panel$graze) + 1)/2 


  ggplot(df_panel, aes(x = graze, y = c3_percent)) +
    geom_col(aes(fill = c3), position = 'stack') +
    adjust_graze_x_axis() +
    scale_fill_manual(values = cols_c3, name = lab_c3) +
    scale_y() +
    geom_text(# to avoid overplotting the word many times
              data = distinct(df_panel[, 'rcp_year']), 
              mapping = aes(x = x_text, y = 85, label = rcp_year),
              family = "sans",
              size = 8, size.unit = 'pt', fontface = 1) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      # shrink subtitle
      plot.subtitle = element_text(size = 9, margin = margin(b = 1)),
      # move facet strip text closer to plot
      # strip.text = element_text(size = 8, margin = margin(b = 1, t = 1),
      #                           vjust = -1),
      strip.text = element_blank(),
      # reduce panel spacing between facets
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      # shrink plot margins
      plot.margin = margin(2, 3, 2, 2),
      theme(text = element_text(face = "plain"))
    )+
    facet_wrap(~rcp_year) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(y = lab_areaperc0)
}

# main builder for one RCP across all regions 
#' @title Build stacked SEI class panels for all regions
#' @description Builds a patchwork of stacked_panel plots for each region, showing Current vs selected RCP scenario.
#' @param rcp Character, the RCP scenario to plot (must be present in \code{sei_pcent1}).
#' @return A patchwork object combining all regions  with a shared legend.
make_stacked_c3_panels <- function(df, rcp, vr, include_sec_y = FALSE) {
  
  # filter to median summary and to Current vs selected rcp
  df0 <- df %>%
    filter(summary == "median",
           RCP %in% c("Current", rcp)) %>%
    arrange(region, graze, RCP, years)
  
  # collect region order & labels
  regions <- levels(df0$region)
  region_label <- region_label_factory(v = vr)
  # build per-region patchwork rows
  plots <- purrr::map(regions, function(rr) {
    
    # split to current and selected future for this region
    d_reg <- filter(df0, region == rr)
    if(include_sec_y) {
      area_total <- unique(d_reg$area_region)
      stopifnot(length(area_total) == 1)
    } else {
      area_total <- NULL
    }

    # two panels for this region
    g1 <- stacked_panel(df_panel = d_reg, area_total = area_total) 
    
    title <- region_label(rr)
    g1 + 
      labs(subtitle = title)
  })
  
  if(include_sec_y){
    if(nr == 9) {
      remove_x = 1:6
      remove_y = c(2, 3, 5, 6, 8, 9)
    } else if (nr == 5) {
      remove_x = c(2, 3)
      remove_y = c(2, 3, 5)
    } else {
      remove_x = NULL
      remove_y = NULL
    }
    plots <- map(plots, \(x) x + labs(x = NULL, y = NULL))
    plots <- remove_axis_labels(plots, remove_x = remove_x, 
                                remove_y = remove_y)
    ylab_sec <- lab_areaha0
    xlab <-  lab_graze
    ylab <-  lab_areaperc0b
    axes <- NULL
  } else {
    ylab_sec <- NULL
    xlab <-  NULL
    ylab <-  NULL
    axes <- 'collect'
  }

  comb <- combine_panels_labs(plots, xlab = xlab, ylab = ylab, 
                      ylab_sec = ylab_sec,
                      axes = axes) 
  
  comb
  
}

# dot plots ---------------------------------------------------------------


#' @title Build dotplot of burned area
#' @description Creates dotplot  region × climate facet (e.g. Current or RCP scenario),
#'
#' @param df_panel Data frame filtered to one region and one climate facet. 
#' @param ylim y limits
fire_dot_panel <- function(df_panel, ylim) {
  
  
  x_text <- (lu(df_panel$graze) + 1)/2 
  y_text <- 0.95*max(ylim)
  df_text <- df_panel %>% 
    group_by(rcp_year) %>% 
    summarize(y_max = max(c(area_high_perc, area_median_perc),
                          na.rm = TRUE)) %>% 
    mutate(y_text = ifelse(y_max > !!y_text, 0.1*max(ylim), y_text),
           x_text = x_text)
  
  ggplot(df_panel, aes(x = graze)) +
    geom_errorbar(aes(ymin = area_low_perc, ymax = area_high_perc),
                  width = 0) +
    geom_point(aes(y = area_median_perc), 
               position = position_dodge(width = 0.5)) +
    adjust_graze_x_axis() +
    geom_text(# to avoid overplotting the word many times
      data = df_text, 
      mapping = aes(x = x_text, y = y_text, label = rcp_year),
      family = "sans",
      size = 8, size.unit = 'pt', fontface = 1) +
    theme(
      # shrink subtitle
      plot.subtitle = element_text(size = 10, margin = margin(b = 1)),
      strip.text = element_blank(),
      # reduce panel spacing between facets
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      # shrink plot margins
      plot.margin = margin(2, 3, 2, 2)
    )+
    facet_wrap(~rcp_year) +
    labs(y = lab_ba1) +
    coord_cartesian(ylim = ylim)
}

#' @title Build  fire dotplots, 1 panel per region
#'
#' @param df dataframe
#' @param vr region verion
#' @param rcp Character, the future RCP scenario to plot .
make_firedot_panels <- function(df, vr, rcp) {
  
  # filter to median summary and to Current vs selected rcp
  df0 <- df %>%
    filter(RCP %in% c("Current", rcp)) %>%
    arrange(ecoregion, graze, RCP, years)
  
  max <- max(with(df0, c(area_high_perc, area_median_perc)),
             na.rm = TRUE)
  ylim <- c(0, max)
  
  # collect region order & labels
  regions <- levels(df0$ecoregion)
  region_label <- region_label_factory(v = vr)
  # build per-region patchwork rows
  plots <- purrr::map(regions, function(rr) {
    
    # split to current and selected future for this region
    d_reg <- filter(df0, ecoregion == rr)
    
    # two panels for this region
    g1 <- fire_dot_panel(df_panel = d_reg,ylim = ylim) 
    
    title <- region_label(rr)
    g1 + 
      labs(subtitle = title)
  })
  
  
  comb <- combine_panels_labs(plots, xlab = NULL, ylab = NULL, 
                              ylab_sec = NULL,
                              axes = 'collect') 
  
  comb
  
}

# bar charts --------------------------------------------------------------

base_c9_area <- function(pattern_var = 'graze', legend_title = 'Grazing') {
  out <- list(
    bar_pattern = ggpattern::geom_bar_pattern(aes(pattern = .data[[pattern_var]],
                                       pattern_density = .data[[pattern_var]],
                                       pattern_angle = .data[[pattern_var]],
                                       fill = c9),
                                   stat = 'identity',
                                   position = position_dodge(),
                                   pattern_fill = '#636363',
                                   pattern_color = '#636363',
                                   pattern_spacing = 0.02,
                                   color = 'white',
                                   pattern_key_scale_factor = 0.5# relative density in the legend
    ),
    scale_fill_manual(values = c9Palette, guide = 'none'),
    ggpattern::scale_pattern_manual(values = c("stripe", "none", "stripe", "stripe"),
                                    name = legend_title),
    ggpattern::scale_pattern_density_manual(values = rep(0.01, 4),
                                            name = legend_title),
    ggpattern::scale_pattern_angle_manual(values = c(45, 0, 0, -45),
                                          name = legend_title),
    theme(axis.text.x = element_text(angle = 75, hjust = 1),
          legend.position = 'right')
  )

  out
}




# functions for tradeoff figure -------------------------------------------

base_tradeoff <- function(group = 'rcp_year', linetypes_scen = c(1, 1, 1),
                          xlab = '% Core Sagebrush Area') {
  list(
    geom_path(aes(group = .data[[group]], linetype = rcp_year), color = "blue", 
              linewidth = 0.5, alpha = 0.5),
    #geom_path(aes(group = graze, color = graze), linewidth = 0.5, alpha = 0.5)
    geom_point(aes(color = graze, shape = rcp_year)),
    scale_color_manual(values = cols_graze, name = lab_graze),
    scale_shape_manual(values = shapes_scen, name = 'Scenario'),
    scale_linetype_manual(name = 'Scenario', values = linetypes_scen),
    labs(x = xlab,
         y = 'Expected area burned (%/year)'),
    expand_limits(x = 0)
  )
}
# simple barchart (for inset) of percent area by SEI class
area_bar <- function(df_smry_region, area_perc_region, limits) {
  
  stopifnot(lu(df_smry_region$region) == 1,
            lu(df_smry_region$area) == 3)
  
  
  ggplot(area_perc_region, aes(c3, area_perc, fill = c3)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c3Palette) +
    labs(x = NULL,
         y = lab_areaperc0) +
    guides(fill = 'none') +
    theme(axis.ticks = element_blank(),
          axis.title.y = element_text(size = 7),
          axis.text = element_blank(),
          plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(ylim = limits)
}

# mapping simple c3 area bar chart over regions
area_bar_map <- function(df_smry) {
  
  # percent of total area in each SEI class
  area_perc <- df_smry %>% 
    group_by(region, c3) %>% 
    summarize(area = unique(area),
              .groups = 'drop_last') %>% 
    mutate(area_perc = area/sum(area)*100)
  limits <- c(0, max(area_perc$area_perc))
  
  regions <- levels(df_smry$region) %>% 
    setNames(nm = .)
  
  map(regions, function(region) {
    df_smry_region <- df_smry %>% 
      filter(region == !!region)
    area_perc_region <- area_perc %>% 
      filter(region == !!region)
    area_bar(df_smry_region, area_perc_region, limits)
  })
}

# make the tradeoff lines/point plot 
tradeoff_lines <- function(df_smry_region, df_gcm_region, xlim, ylim,
                           v,
                           gcm_path = TRUE) {
  
  cols2 <- setNames(c3Palette, relable_c3_current(names(c3Palette)))
  region <- unique(df_smry_region$region)
  stopifnot(lu(region) == 1)
  
  g <- ggplot(df_smry_region, aes(SEI_mean, expected_ba_perc, color = c3_cur)) 
  
  if(gcm_path) {
    g <- g + geom_path(data = df_gcm_region, aes(group = rcp_year_c3_gcm, linetype = rcp_year),
                       linewidth = 0.2, alpha = 0.5) 
  }
  region_labeller <- region_labeller_factory(v = vr)
  subtitle <- region_labeller(region)[[1]]
  g +  
    geom_path(aes(group = rcp_year_c3, linetype = rcp_year),
              linewidth = 1, alpha = 1) +
    scale_color_manual(values = c(cols2), name = 'Historical SEI class') +
    ggnewscale::new_scale_color() +
    geom_point(aes(color = graze, shape = rcp_year)) +
    scale_color_manual(values = cols_graze, name = 'Grazing scenario') +
    scale_linetype_manual(name = 'Climate scenario', values = unname(linetypes_scen)[1:2]) +
    scale_shape_manual(values = unname(shapes_scen)[1:2], name = 'Climate scenario') +
        labs(x = 'Mean SEI',
         y = 'Expected burned area (%/year)',
         subtitle = subtitle)  +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    make_legend_small()
}

# map over region
tradeoff_lines_map <- function(df_smry, df_gcm, v, gcm_path = TRUE) {
  x <- df_smry$SEI_mean
  y <- df_smry$expected_ba_perc
  if(gcm_path) {
    x <- c(x, df_gcm$SEI_mean)
    y <- c(y, df_gcm$expected_ba_perc)
  }
  xlim <- range(x)
  ylim <- range(y)
  regions <- levels(df_smry$region) %>% 
    setNames(nm = .)
  
  map(regions, function(region) {
    df_smry_region <- df_smry %>% 
      filter(region == !!region)
    df_gcm_region <- df_gcm %>% 
      filter(region == !!region)
    tradeoff_lines(df_smry_region, df_gcm_region, xlim = xlim, ylim = ylim,
                   v = v, gcm_path = gcm_path)
  })
}


tradeoff_add_inset <- function(plot, inset) {
  plot + inset_element(inset, 
                       left = 0.73,
                       bottom = 0.73,
                       right = 0.98,
                       top = 0.98,
                       align_to = 'panel')
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


combine_5_panels <- function(plots, axes = 'collect') {
  stopifnot(length(plots) == 5)
  l1 <- c(plots[1:3], list(patchwork::guide_area()),
          plots[4:5])
  wrap_plots(l1, nrow = 2, guides = 'collect') +
    # thise axis collect only works if plots are not
    # already patchwork combinations (e.g. plot with inset)
    plot_layout(axes = axes,
                axis_titles = 'collect')
}

combine_more_panels <- function(plots, axes = 'collect',
                                guides = 'collect') {
  stopifnot(length(plots) > 5)
  l1 <- plots
  comb <- wrap_plots(l1, nrow = 3, guides = guides) +
    # thise axis collect only works if plots are not
    # already patchwork combinations (e.g. plot with inset)
    plot_layout(axes = axes,
                axis_titles = 'collect')
  
  comb&theme(legend.position = 'right')
}

# combine panels, and add seperate plots that are
# the y and x axis labels
# this requires the axis labels of the individual
# plots to already have been removed as needed
combine_panels_labs <- function(plots, xlab = NULL, ylab = NULL, ylab_sec = NULL, ...) {

  
  f <- if(length(plots) == 5) combine_5_panels else combine_more_panels
  panels <- f(plots, ...)
  
  if(is.null(xlab) & is.null(ylab) & is.null(ylab_sec)) {
    return(panels)
  }
  y <- ylab_plot(ylab)
  x <- xlab_plot(xlab)
  
  # Combine y label + panels
  if(!is.null(ylab_sec)) {
    y_sec <- ylab_plot(ylab_sec, secondary = TRUE)
    if(length(plots) == 5) {
      main <- wrap_plots(y, panels, y_sec, ncol = 3, widths = c(0.02, 1, 0.02))
    } else {
      panels <- f(plots, guides = NULL, ...)
      main <- patchwork::wrap_plots(
        y, panels, y_sec, patchwork::guide_area(),
        ncol = 4,
        widths = c(0.02, 1, 0.02, 0.14),
        guides = "collect"
      ) & theme(legend.position = "right")
    }

  } else {
    main <- wrap_plots(y, panels, ncol = 2, widths = c(0.02, 1))
  }
  # Combine with x label
  final <- wrap_plots(main, x, ncol = 1, heights = c(1, 0.02))
  final
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

# remove axis labels for elements of the list of plots
# based on their index,
# defaults work for 5 panel region layout
remove_axis_labels <- function(plots, 
                               remove_x = c(2, 3),
                               remove_y = c(2, 3, 5)) {
  plots[remove_x] <- map(plots[remove_x], function(g) {
    g + theme(axis.text.x = element_blank())
  })
  
  plots[remove_y] <- map(plots[remove_y], function(g) {
    g + theme(axis.text.y.left = element_blank())
  })
  # remove all axis titles
  plots <- map(plots, function(g) {
    g + labs(x = NULL, y = NULL)
  })
  
}

xlab_plot <- function(xlab = '', ...) {
  ggplot() + 
    theme_void() + 
    labs(x = xlab) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      axis.title.x = element_text(vjust = 0.5),
      ...
    )
}

ylab_plot <- function(ylab = '', secondary = FALSE, ...) {
  ggplot() +
    theme_void() +
    labs(y = ylab) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      axis.title.y = element_text(
        angle = if (secondary) 270 else 90
      ),
      ...
    )
}


tag <- function(g, tag_pool) {
  egg::tag_facet(g, tag_pool = tag_pool,
                 open = "",
                 close = "",
                 hjust = -0.1,
                 size = 3) + 
    theme(strip.text = element_text())
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


# crossplots --------------------------------------------------------------

# make a single panel, plantting IQR of two variables
# (for plotting climate variables agains each other)
crossplot_1panel <- function(df_wide, var1, var2, linewidth =0.3,
                             colors = cols_GCM1,
                             shapes = shapes_GCM1) {
    g <- ggplot(data = df_wide) 
  # adding vertical and horizontal lines, as long as it doesn't
  # extend the axes
  if(min(df_wide[[paste0(var1, '_p25')]]) < 0) {
    g <- g + geom_vline(
      mapping = aes(xintercept = 0),
      alpha = 0.3, linetype = 2, linewidth = linewidth
    ) 
  }
  if(min(df_wide[[paste0(var2, '_p25')]]) < 0) {
    g <- g + geom_hline(
      mapping = aes(yintercept = 0),
      alpha = 0.3, linetype = 2, linewidth = linewidth
    ) 
  }
  labeller <- driver_labeller(delta = TRUE)
  g <- g + 
    geom_segment(aes(x = .data[[paste0(var1, '_median')]], 
                     y = .data[[paste0(var2, '_p25')]], 
                     yend = .data[[paste0(var2, '_p75')]],
                     color = GCM), linewidth = linewidth) +
    geom_segment(aes(x = .data[[paste0(var1, '_p25')]], 
                     xend = .data[[paste0(var1, '_p75')]], 
                     y = .data[[paste0(var2, '_median')]],
                     color = GCM), 
                 linewidth = linewidth) +
    geom_point(aes(x = .data[[paste0(var1, '_median')]], 
                   y = .data[[paste0(var2, '_median')]], 
                   color = GCM, shape = GCM)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(x = labeller(var1),
         y = labeller(var2))
  g

}

# facet functions ---------------------------------------------------------

facet_manual_region <- function(v = 'r1.0', 
                                legend.position = NULL,
                                legend.position.inside = NULL,
                                color_strips = FALSE,
                                # letters to use for facet prefixes
                                region_letters = NULL,
                                scales = 'fixed') {
  
  if(v == 'r1.0') {
    design <- "
               ABC
               #DE
               "
    if(is.null(legend.position)) legend.position <-  'inside'
    if(is.null(legend.position.inside)) legend.position.inside = c(0.05, 0)

  } else if(v == 'r1.1') {
    design <- "
               ABC
               DEF
               GHI
               "
    legend.position = 'right'
  } else {
    stop("function doesn't support this version of regions, must update")
  }

  
  if(color_strips) {
    # ordering, and only keeping the appropriate 
    cols <- get_cols_ecoregion(v = v)
    strip_theme_region <- ggh4x::strip_themed(
      text_x = ggh4x::elem_list_text(color = c('black', cols))
    )
  } else {
    strip_theme_region <- ggh4x::strip_vanilla()
  }
  region_labeller <- region_labeller_factory(region_letters, v = v)
  list(
    ggh4x::facet_manual(~region, design = design, labeller = region_labeller,
                        strip = strip_theme_region,
                        scales = scales),
    theme(legend.position = legend.position, 
          legend.position.inside = legend.position.inside,
          strip.text = element_text(hjust = 0))
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

adjust_graze_x_axis <- function(xlab = lab_graze) {
  list(
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    labs(x = xlab))
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

region_label_factory <- function(region_letters = NULL,
                                 v = NULL) {
  if(is.null(region_letters)) {
    region_letters <- fig_letters # from the global environment
  }
  
  function(x) {
    x <- region_factor(x, v = v)
    levels <- levels(x)
    labels <- paste(region_letters[1:length(levels)], levels)
    out <- factor(as.character(x), levels = levels,
                  labels = labels)
    as.character(out)
  }

}

region_labeller_factory <- function(region_letters = NULL,
                                    v = NULL) {
  region_label <- region_label_factory(region_letters = region_letters,
                                       v = v)
  ggplot2::as_labeller(x = region_label)
}

region_labeller <- region_labeller_factory(fig_letters)
  
driver_labeller <- function(delta = FALSE) {
  lookup_md <- c(
    "MAT" = "MAT (\u00b0C)",
    "MAP" = "MAP (mm)",
    "PSP" = "PSP",
    "psp" = "PSP",
    "Aherb" = "Aherb (g/m<sup>2</sup>)",
    "Pherb" = "Pherb (g/m<sup>2</sup>)",
    "Sagebrush" = "Sagebrush (g/m<sup>2</sup>)"
  )
  if(delta) {
    nms <- names(lookup_md)
    lookup_md <- paste("\u0394", lookup_md)
    names(lookup_md) <- nms
  }
  
  as_labeller(lookup_md)
}

#' create legend labels
#'
#' @param x numeric vector of break points
#'
#' @return character vector, where last category is just
#' > x[n-1] instead of showing a range

#' @examples
#' label_creator(1:5)
label_creator <- function(x, convert2percent = FALSE) {
  if(convert2percent) x <- x*100
  
  n <- length(x)
  labels <- vector(mode = 'character', length = n-1)
  
  for (i in 1:(n-1)) {
    if(i < n -1) {
      labels[i] <- paste(x[i], "to", x[i+1])
    } else {
      labels[i] <- paste(">", x[i])
    }
  }
  labels[1] <- paste0('< ', x[2])
  
  labels
}


# * plots -----------------------------------------------------------------

# 'ggplot' label plot for left side
label_left_plot <- function(x) {
  ggplot() +
    theme_void() +
    annotate("text", x = 1, y = 1, label = x, angle = 90, size = 3) +
    theme(plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
}

# 'ggplot' label plot for top
label_top_plot <- function(x) {
  ggplot() +
    theme_void() +
    annotate("text", x = 1, y = 1, label = x, size = 3) +
    theme(plot.margin = unit(c(0, 0, 0, 0), units = 'in'))
}


# legends -----------------------------------------------------------------

make_legend_small <- function(pointSize = 2, textSize = 7, spaceLegend = 0.2,
                              lineLength = 20) {
  list(guides(shape = guide_legend(override.aes = list(size = pointSize)),
              color = guide_legend(override.aes = list(size = pointSize))),
       theme(legend.title = element_text(size = textSize), 
             legend.text  = element_text(size = textSize),
             legend.key.size = unit(spaceLegend, "lines"),
             legend.key.width = unit(lineLength, 'pt')))
  
}

smaller_legend <- function() {
  theme(
    legend.title = element_text(size = rel(0.7)),  # Reduce legend title size
    legend.text = element_text(size = rel(0.7)),   # Reduce legend text size
    legend.key.size = unit(0.7, "lines"),          # Reduce legend key size (relative to lines)
    legend.spacing = unit(0.7, "lines")           # Reduce spacing between legend items
  )
}

# color functions ---------------------------------------------------------

# grazing colors (cols_graze defined in the fig_params.R script)
scale_color_graze <- function() {
  scale_color_manual(values = cols_graze, name = "Grazing")
}

scale_fill_smry <- function() {
  scale_fill_manual(values = cols_smry,
                    name = 'Summary across GCMs')
}

scale_fill_graze <- function(exclude = NULL) {
  col <- cols_graze
  if(!is.null(exclude)) {
    stopifnot(exclude %in% names(col))
    col <- col[exclude != names(col)] 
  }

  scale_fill_manual(values = col, name = "Grazing")
}

scale_fill_c9 <- function(...) {
  scale_fill_manual(values = unname(c9Palette),
                    labels = names(c9Palette),
                    na.value = 'transparent',
                    na.translate = FALSE,
                    name = 'SEI class change',
                    drop = FALSE,
                    ...)
}

scale_fill_c3 <- function(...) {
  scale_fill_manual(values = unname(c3Palette),
                    na.value = 'transparent',
                    labels = names(c3Palette),
                    na.translate = FALSE,
                    name = 'SEI class',
                    ...)
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


ggsave_tradeoff <- function(g, prefix, width = 7, height = 4.5,
                            xvar = "csa") {
  ggsave(paste0("figures/sei/tradeoff/", xvar, "-scd-adj-vs-ba_perc_dotplot_", 
                prefix, '_', suffix, ".png"), 
         plot = g, dpi = 600,
         width = width, height = height)
}

ggsave_delta_prob <- function(plot, rcp, smry) {
  filename <- paste0('figures/fire/delta-prob/delta-prob_vs_delta-pred_', 
                     smry, "_", vr_name, rcp, yr_lab, "_", runv, '.png')
  ggsave(
    filename = filename,
    plot = plot,
    width = 13, height = 10*mr,
    dpi = 600
  )
}

# create 9 color matrix ---------------------------------------------------
# Creating a 3x3 colored matrix of current and future SEI classes,


color_matrix <- function(xlab = 'Future',
                         ylab = 'Current') {
  
  
  c3_levels <- c('CSA', 'GOA', 'ORA')
  
  # c9 levels, with there associated current and future SEI categories
  df_c9 <- tibble(
    c9Name  = c9Names, 
    c9Value = 1:9,
    current = rep(c3_levels, each = 3), # Current SEI (3 levels)
    future = rep(c3_levels, 3) # future SEI (3 levels)
  ) %>% 
    mutate(
      current = factor(current, levels = rev(c3_levels)),
      future = factor(future, levels = c3_levels)
    )
  
  # Adding label column (category of change)
  df_c9$label <- NA
  df_c9$label[c(1, 5, 9)] <- "Stable"
  df_c9$label[c(2, 3, 6)] <- "Decline"
  df_c9$label[c(4, 7, 8)] <- "Increase"
  
  # color of text in color matrix
  text_color <- rep('black', 9)
  text_color[c(1, 2, 3, 4, 7)] <- 'white' # background is dark
  names(text_color) <- c9Names
  stopifnot(names(c9Palette) == df_c9$c9Name)
  c9Palette
  text_color
  c9Names
  ggplot(df_c9, aes(future, current)) +
    geom_tile(aes(fill = c9Name)) +
    geom_text(aes(label = label, color = c9Name), size = 2) +
    theme_minimal() +
    scale_x_discrete(position = 'top') +
    scale_fill_manual(values = c9Palette) +
    labs(x = xlab,
         y = ylab) +
    scale_color_manual(values = text_color)+
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          text = element_text(size = 8),
          # alternatively make the background white
          plot.background = element_rect(fill = 'transparent', 
                                         color = 'transparent'),
          plot.margin = unit(c(0, 0, 0, 0), "in"),
          aspect.ratio = 1)
}


