# Functions in this script written by Daniel Schlaepfer

#----- Load geographic data ------

#' @export
load_sagebiome_polygon <- function(path, crs) {
  newRR3::import_sagebiome_polygon(dir_dataraw = path) |>
    sf::st_read() |>
    sf::st_transform(crs = crs) |>
    sf::st_make_valid()
}


#' @export
load_NEadmin1_polygon <- function(path, crs) {
  newRR3::import_NEadmin1_polygon(dir_dataraw = path) |>
    readRDS() |>
    sf::st_transform(crs = crs) |>
    sf::st_make_valid()
}


#' @export
load_ecoregions_polygon <- function(path, crs) {

  res <- newRR3::import_EPAecoregionL3_polygon(dir_dataraw = dir_dataraw_tu) |>
    readRDS() |>
    sf::st_transform(crs = crs_map) |>
    # union polygons by L3 ecoregions
    dplyr::group_by(!!rlang::sym("NA_L3NAME")) |>
    dplyr::summarize(do_union = TRUE) |>
    sf::st_make_valid()

  fname_aggecoregions <- file.path(
    path,
    "EPAecoregionL3",
    "EcoregionsL3_Simplified.csv"
  )

  if (file.exists(fname_aggecoregions)) {
    agged_ecoregions <- data.frame(
      utils::read.csv(fname_aggecoregions),
      stringsAsFactors = FALSE
    )

    merge(
      res,
      agged_ecoregions,
      by = "NA_L3NAME",
      all.x = TRUE
    )

  } else {
    res
  }
}


#------ Colors ------

#' Color used for `no data` in Chambers et al. 2023
colors_NA_Chambers2023 <- function() {
  "#ffffbf"
}

#' Colors used for `R&R` levels in Chambers et al. 2023
colors_RR_Chambers2023 <- function(rr_levels = c("L", "ML", "M", "H+MH")) {
  stopifnot(requireNamespace("RColorBrewer"))
  res <- RColorBrewer::brewer.pal(length(rr_levels), "RdYlBu")
  names(res) <- rr_levels
  res
}

#' Okabe-Ito colors
#' @export
colors_okabeito <- function() {
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

#' Colors for a 4 x 4 convolution matrix
colors_4x4 <- function() {
  stopifnot(requireNamespace("biscale"))

  # purple-green scale
  col_4biscale <- biscale::bi_pal(
    pal = "PurpleGrn",
    dim = 4,
    flip_axes = FALSE,
    rotate_pal = TRUE,
    preview = FALSE
  )

  # set all ref == scen (no change) to gray
  col_4biscale[c("1-1", "2-2", "3-3", "4-4")] <- col_4biscale[["4-4"]]

  # rotate colors so that darkest ones are used for largest change
  tmpbs <- col_4biscale
  from_to <- rbind(
    # upper-left triangle
    c(from = "1-2", to = "1-4"),
    c("2-3", "1-2"),
    c("2-4", "2-3"),
    c("1-4", "2-4"),
    # lower-right triangle
    c("2-1", "4-1"),
    c("3-2", "2-1"),
    c("4-2", "3-2"),
    c("4-1", "4-2")
  )
  for (kft in seq_len(nrow(from_to))) {
    col_4biscale[from_to[kft, "to"]] <-
      tmpbs[from_to[kft, "from"]]
  }

  if (FALSE) {
    biscale::bi_pal(col_4biscale, dim = 4, preview = TRUE)
  }

  col_4biscale
}


#' Color scale to display `R&R` indicators
#' @export
colorscale_rr <- function(
  name,
  levels,
  aesthetics = "fill",
  ...
) {
  stopifnot(requireNamespace("ggplot2"))

  ggplot2::scale_fill_manual(
    name = if (missing(name)) NULL else name,
    aesthetics = aesthetics,
    values = colors_RR_Chambers2023(levels),
    breaks = levels,
    na.value = colors_NA_Chambers2023(),
    ...
  )
}

#' Add a `no data` (mask) color to legend
#' @references \url{https://stackoverflow.com/a/42380983}
#' @export
add_NA_colorlegend <- function() {
  #stopifnot(requireNamespace("ggplot2"))

  # trick ggplot into drawing a separate color scale by using aes() to
  # map a column of empty strings to color
  ggplot2::geom_point(
    data = data.frame(X = 0:2, Y = 0:2),
    ggplot2::aes(x = "X", y = "Y", color = ""),
    show.legend = TRUE
  ) +
    # don't draw colors
    ggplot2::scale_color_manual(values = c("No data" = NA)) +
    # ensure the new color legend is drawn with the correct color
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1
      ),
      color = ggplot2::guide_legend(
        order = 2,
        override.aes = list(fill = colors_NA_Chambers2023())
      )
    )
}


#' Continous approach to coloring RR indicators
#' @export
colorscale_rr2 <- function(
  name,
  levels,
  aesthetics = "fill",
  ...
) {
  stopifnot(requireNamespace("ggplot2"))

  ggplot2::scale_fill_gradientn(
    name = if (missing(name)) NULL else name,
    limits = c(1, length(levels)),
    labels = factor(levels, levels = levels),
    colors = colors_RR_Chambers2023(levels),
    na.value = colors_NA_Chambers2023(),
    ...
  )
}


colors_robustsimpledelta <- function() {
  tmp <- colors_4x4()

  c(
    robust_decrease = tmp[["3-1"]],  #tmp[["4-1"]],
    nonrobust_decrease = tmp[["3-2"]],
    robust_stable = "gray50", # "gray60",
    nonrobust_stable = "gray80", # "gray90", #tmp[["1-1"]],
    robust_increase = tmp[["1-3"]], #tmp[["1-4"]],
    nonrobust_increase = tmp[["2-3"]]
  )
}

#' Color scale to display robust simple deltas with `ggplot2`
#'
#' @export
colorscale_robustsimpledelta <- function(
  name,
  levels,
  aesthetics = "fill",
  ...
) {
  stopifnot(requireNamespace("ggplot2"))

  level_colors <- colors_robustsimpledelta()

  stopifnot(names(level_colors) %in% levels)

  ggplot2::discrete_scale(
    aesthetics = aesthetics,
    scale_name = "manual",
    name = if (missing(name)) NULL else name,
    limits = factor(levels, levels = levels),
    labels =
      levels |>
      gsub(pattern = "_", replacement = " ", x = _, fixed = TRUE) |>
      gsub(pattern = "nonrobust", replacement = "other", x = _, fixed = TRUE),
    palette = function(n) level_colors,
    na.value = colors_NA_Chambers2023(),
    drop = FALSE,
    ...
  )
}


#' Color scale to display robust deltas grouped by original level with `ggplot2`
#'
#' @export
colorscale_robustdelta_by_ref <- function(
  name,
  levels,
  aesthetics = "fill",
  ...
) {
  stopifnot(requireNamespace("ggplot2"))

  tmp <- colors_robustsimpledelta()

  tmp_robust_levels <-
    strsplit(levels, split = "_", fixed = TRUE) |>
    vapply(function(x) paste0(x[[1L]], "_", x[[2L]]), FUN.VALUE = NA_character_)

  ids <- match(tmp_robust_levels, names(tmp))
  level_colors <- tmp[ids]
  names(level_colors) <- levels

  ggplot2::scale_fill_manual(
    name = if (missing(name)) NULL else name,
    aesthetics = aesthetics,
    values = level_colors,
    breaks = factor(levels, levels = levels),
    labels = gsub("_", " ", levels, fixed = TRUE),
    na.value = colors_NA_Chambers2023(),
    ...
  )
}


#' Color scale to display a combination of fixed/flat color and continous scale
#'
#' @param position_flat
#'   1: values between limits[1] and threshold: fixed `color_flat`
#'      values between threshold and limits[2]: viridis "plasma"
#'   2: values between limits[1] and threshold: viridis "plasma"
#'      values between threshold and limits[2]: fixed `color_flat`
#' @export
colorscale_flatcontinuous <- function(
  name,
  limits,
  threshold,
  aesthetics = "fill",
  color_flat = "antiquewhite",
  position_flat = 1,
  viridis_pal_option = "C", # C = plasma; D = viridis,
  reverse = FALSE,
  begin_hue = 0,
  end_hue = 1,
  ...
) {
  position_flat <- as.integer(position_flat)

  stopifnot(
    threshold >= limits[1],
    threshold <= limits[2],
    position_flat %in% 1:2,
    requireNamespace("viridis")
  )

  vfix <- if (position_flat == 1L) {
    unique(c(limits[1], threshold))
  } else if (position_flat == 2L) {
    unique(c(threshold, limits[2]))
  }
  nfix <- length(vfix)

  colors <- c(
    if (position_flat == 1L) rep(color_flat, nfix),
    viridis::viridis_pal(
      alpha = 1,
      begin = begin_hue,
      end = end_hue,
      direction = if (xor(position_flat == 1L, reverse)) 1 else -1,
      option = viridis_pal_option
    )(6),
    if (position_flat == 2L) rep(color_flat, nfix)
  )

  color_values <- if (position_flat == 1L) {
    c(
      vfix,
      threshold + sqrt(.Machine[["double.eps"]]),
      seq(threshold, limits[2], length.out = 6)[-1]
    )
  } else {
    c(
      seq(limits[1], threshold, length.out = 6)[-6],
      threshold - sqrt(.Machine[["double.neg.eps"]]),
      vfix
    )
  }

  ggplot2::scale_fill_gradientn(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    limits = limits,
    na.value = "gray",
    oob = scales::oob_squish,
    colors = colors,
    # this vector gives the position (between 0 and 1) for each color in the colors vector
    values = scales::rescale(color_values, to = c(0, 1)),
    ...
  )
}


#' Binned color scale where part of the scale is flat (one consistent color)
#'
#' @param breaks Integer vector with breaks including lower and upper limit.
#' @export
colorscale_flatbinned <- function(
  name,
  breaks,
  limits = range(breaks),
  threshold,
  aesthetics = "fill",
  color_flat = "antiquewhite",
  position_flat = 1,
  viridis_pal_option = "C", # C = plasma; D = viridis,
  reverse = FALSE,
  begin_hue = 0,
  end_hue = 1,
  ...
) {
  position_flat <- as.integer(position_flat)

  stopifnot(
    threshold >= limits[1],
    threshold <= limits[2],
    position_flat %in% 1:2,
    requireNamespace("viridis")
  )

  colors <- vector("character", length = length(breaks))

  is_flat <- if (position_flat == 1L) {
    breaks <= threshold
  } else if (position_flat == 2L) {
    breaks >= threshold
  }

  colors[is_flat] <- color_flat

  colors[!is_flat] <- viridis::viridis_pal(
    alpha = 1,
    begin = begin_hue,
    end = end_hue,
    direction = if (xor(position_flat == 1L, reverse)) 1 else -1,
    option = viridis_pal_option
  )(sum(!is_flat))

  ggplot2::scale_fill_stepsn(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    limits = limits,
    na.value = "gray",
    oob = scales::oob_squish,
    colors = colors,
    breaks = breaks,
    ...
  )
}



#' @export
colorscale_pfirst <- function(
  name,
  aesthetics = "fill",
  limits = NULL,
  ...
) {
  colorscale_flatcontinuous(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    limits = c(0, 1),
    threshold = 0.25,
    position_flat = 1,
    ...
  )
}


#' @export
colorscale_novelty_flatbinned <- function(
  name,
  aesthetics = "fill",
  limits = c(0, 3),
  ...
) {
  colorscale_flatcontinuous(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    limits = limits,
    threshold = 1,
    position_flat = if (FALSE) 2L else 1L,
    ...
  )
}

#' @export
colorscale_novelty_diverging <- function(
  name,
  aesthetics = "fill",
  limits = c(0, 3),
  mid = 1,
  ...
) {
  colorspace::scale_fill_continuous_diverging(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    palette = "Vik",
    na.value = "gray",
    limits = limits,
    mid = mid,
    ...
  )
}

#' @export
colorscale_values_c <- function(
  name,
  aesthetics = "fill",
  limits = c(NA, NA),
  ...
) {
  ggplot2::scale_fill_viridis_c(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    limits = limits,
    na.value = "gray",
    ...
  )
}


#' @export
colorscale_deltas_c <- function(
  name,
  aesthetics = "fill",
  limits = c(NA, NA),
  mid = 0,
  # TRUE: green:decrease <> purple:increase
  # FALSE: purple:decrease <> green:increase
  rev = TRUE,
  ...
) {
  colorspace::scale_fill_continuous_diverging(
    aesthetics = aesthetics,
    name = if (missing(name)) NULL else name,
    palette = "Purple-Green",
    rev = rev,
    na.value = "gray",
    mid = mid,
    limits = limits,
    ...
  )
}


#' @export
colorscale_predictors19 <- function(
  name,
  aesthetics = "fill",
  n = 19L,
  ...
) {
  stopifnot(requireNamespace("ggplot2"))

  cols19 <- grDevices::palette.colors(
    n = n,
    palette = "alphabet"
  )

  args <- list(...)
  args[["levels"]] <- NULL # remove a "levels" argument

  args[["aesthetics"]] <- aesthetics
  args[["values"]] <- unname(cols19)
  args[["na.value"]] <- "gray"

  if (missing(name) || is.null(name)) {
    # doesn't work, legend title becomes "NULL": args[["name"]] <- list(NULL)
    args <- c(args, list(name = NULL))
  } else {
    args[["name"]] <- name
  }

  do.call(ggplot2::scale_fill_manual, args)
}


#----- Beautifications ------
#' Create nice labels from across-model tags
#'
#' @examples
#' make_acrmod_labels(c("RCP45_2029-2064", "RCP45_2064-2099"))
#' make_acrmod_labels("historical_1950-2005")
#'
#' @export
make_acrmod_labels <- function(x) {
  strsplit(x, split = "_", fixed = TRUE) |>
    lapply(
      function(x) {
        x1 <- x[[1L]]
        substr(x1, 1L, 1L) <- toupper(substr(x1, 1L, 1L))
        paste0(x1, " (", x[[2L]], ")")
      }
    ) |>
    unlist()
}


create_negmask <- function(poly_not_erase, bbox) {
  xbox <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(bbox[c(1, 2, 1, 4, 3, 4, 3, 2, 1, 2)], ncol = 2, byrow = TRUE)
      )
    ),
    crs = sf::st_crs(poly_not_erase)
  )

  if (!is.null(poly_not_erase)) {
    sf::st_difference(xbox, sf::st_union(sf::st_combine(poly_not_erase)))
  }
}


#' @export
add_tag_as_label <- function(tag) {
  stopifnot(requireNamespace("ggplot2"))

  ggplot2::geom_label(
    data = data.frame(
      xpos = -Inf,
      ypos = Inf,
      vjustvar = 1,
      hjustvar = 0,
      text = tag[[1L]]
    ),
    ggplot2::aes(
      x = xpos,
      y = ypos,
      hjust = hjustvar,
      vjust = vjustvar,
      label = text
    ),
    label.r = grid::unit(0, "lines"),
    fontface = "bold",
    label.size = NA # remove border
  )
}

#' @export
add_tag_as_richlabel <- function(tag) {
  stopifnot(requireNamespace("ggtext"))

  ggtext::geom_richtext(
    data = data.frame(
      xpos = -Inf,
      ypos = Inf,
      vjustvar = 1,
      hjustvar = 0,
      text = tag
    ),
    ggplot2::aes(
      x = xpos,
      y = ypos,
      hjust = hjustvar,
      vjust = vjustvar,
      label = text
    ),
    label.r = grid::unit(0, "lines"),
    fontface = "bold",
    label.size = NA, # remove border
    fill = "#ffffff00" # transparent label background
  )
}

#' @export
add_tag_as_richlabel_rightadjusted <- function(tag) {
  stopifnot(requireNamespace("ggtext"))

  ggtext::geom_richtext(
    data = data.frame(
      xpos = Inf,
      ypos = Inf,
      vjustvar = 1,
      hjustvar = 1,
      text = tag
    ),
    ggplot2::aes(
      x = xpos,
      y = ypos,
      hjust = hjustvar,
      vjust = vjustvar,
      label = text
    ),
    label.r = grid::unit(0, "lines"),
    fontface = "bold",
    label.size = NA # remove border
  )
}



#' Remove inner labels from plot panels
#'
#' @param ks_panels_do_remove An integer vector of panel numbers for which
#' inner labels actually should be removed. `NULL` (default) removes
#' inner labels from all panels.
#' @param remove_xlabs Remove x-axis label and x-axis tick labels.
#' @param remove_xticks Remove x-axis ticks.
#'
#' @export
panels_remove_inner_labels <- function(
  list_panels,
  ks_panels_do_remove = NULL,
  remove_xlabs = NULL,
  remove_xticks = NULL,
  remove_ylabs = NULL,
  remove_yticks = NULL,
  remove_facetlabs = NULL,
  remove_inner_plotmargins = NULL,
  n_panels = NULL,
  byrow = TRUE
) {
  stopifnot(
    requireNamespace("ggplot2"),
    as.logical(byrow)
  )

  if (!is.null(n_panels)) {
    ntmp <- prod(n_panels)

    ids <- if (
      isTRUE(remove_xlabs) ||
        isTRUE(remove_xticks) ||
        isTRUE(remove_inner_plotmargins)
    ) {
      # exclude panels in bottom-most row
      tmp <- seq(from = ntmp - n_panels[[2L]] + 1L, to = ntmp, by = 1L)
      seq_len(ntmp)[-tmp]
    }
    remove_xlabs <- if (isTRUE(remove_xlabs)) ids
    remove_xticks <- if (isTRUE(remove_xticks)) ids
    remove_plotmargins_bottom <- if (isTRUE(remove_inner_plotmargins)) ids

    ids <- if (
      isTRUE(remove_ylabs) ||
        isTRUE(remove_yticks) ||
        isTRUE(remove_inner_plotmargins)
    ) {
      # exclude panels in left-most column
      tmp <- 1L + (seq_len(n_panels[[1L]]) - 1L) * n_panels[[2]]
      seq_len(ntmp)[-tmp]
    }
    remove_ylabs <- if (isTRUE(remove_ylabs)) ids
    remove_yticks <- if (isTRUE(remove_yticks)) ids
    remove_plotmargins_left <- if (isTRUE(remove_inner_plotmargins)) ids

    if (isTRUE(remove_inner_plotmargins)) {
      # exclude panels in top-most row
      tmp <- seq_len(n_panels[[2L]])
      remove_plotmargins_top <- seq_len(ntmp)[-tmp]

      # exclude panels in right-most column
      tmp <- seq_len(n_panels[[1L]]) * n_panels[[2]]
      remove_plotmargins_right <- seq_len(ntmp)[-tmp]
    } else {
      remove_plotmargins_top <- remove_plotmargins_right <- NULL
    }


    if (isTRUE(remove_facetlabs)) {
      stop("Logical argument 'remove_facetlabs' is not implemented.")
    }
  }


  #--- Loop over panels and remove axis elements
  ks_all <- seq_along(list_panels)

  ks_used <- if (is.null(ks_panels_do_remove)) {
    ks_all
  } else {
    ks_panels_do_remove[ks_panels_do_remove %in% ks_all]
  }

  for (k0 in ks_used) {

    if (isTRUE(remove_inner_plotmargins)) {
      tmp <- list_panels[[k0]]$theme$plot.margin
      dflt_plot_margins <- as.numeric(tmp)
      dflt_margin_units <- unique(grid::unitType(tmp))
      stopifnot(length(dflt_margin_units) == 1L)
    } else {
      dflt_plot_margins <- NULL
      dflt_margin_units <- NULL
    }

    # Remove x-axis label and x-axis tick labels
    if (k0 %in% remove_xlabs) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank()
        )
    }

    # Remove x-axis tick marks
    if (k0 %in% remove_xticks) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank()
        )
    }

    # Remove facet labels
    if (k0 %in% remove_facetlabs) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          strip.text.x.top = ggplot2::element_blank()
        )
    }

    # Remove y-axis label and y-axis tick labels
    if (k0 %in% remove_ylabs) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank()
        )
    }

    # Remove y-axis tick marks
    if (k0 %in% remove_yticks) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          axis.ticks.y = ggplot2::element_blank()
        )
    }

    # Remove inner plot margins
    if (isTRUE(remove_inner_plotmargins)) {
      list_panels[[k0]] <- list_panels[[k0]] +
        ggplot2::theme(
          plot.margin = ggplot2::margin(
            t = if (k0 %in% remove_plotmargins_top) 0 else dflt_plot_margins[[1L]],
            r = if (k0 %in% remove_plotmargins_right) 0 else dflt_plot_margins[[2L]],
            b = if (k0 %in% remove_plotmargins_bottom) 0 else dflt_plot_margins[[3L]],
            l = if (k0 %in% remove_plotmargins_left) 0 else dflt_plot_margins[[4L]],
            unit = dflt_margin_units
          )
        )
    }
  }

  list_panels
}


#' Create a polygon from raster that can be plotted as cross-hatching layer
#' @export
create_crosshatching <- function(
  xstars,
  x,
  rat_values,
  dx = units::set_units(4L, "km"),
  fun_crosshatch = function(x) ifelse(x >= 0.9, NA, 1L)
) {
  #--- create stars object of nonrobust area
  tmp <- prepare_stars(
    xstars = xstars,
    x = x,
    rat_values = rat_values
  )
  #--- create spatially lower resolution
  tmpa <- stars::st_warp(
    src = tmp,
    dest = stars::st_as_stars(sf::st_bbox(tmp), dx = dx),
    method = "mode",
    use_gdal = TRUE,
    no_data_value = -.Machine[["double.xmax"]]
  )
  #--- convert to polygon
  tmpa |>
    stars::st_apply(MARGIN = 1:2, FUN = fun_crosshatch) |>
    sf::st_as_sf(as_points = FALSE, merge = TRUE) |>
    sf::st_make_valid()
}



#------ ggplot2 visuals ------
#' Clean theme for `ggplot2`
#' Based on `egg::theme_article()`
#' @export
ggplot2_clean_theme <- function() {
  stopifnot(requireNamespace("ggplot2"))

  line_color <- "gray25"

  ggplot2::theme_bw() +
    ggplot2::theme(
      line = ggplot2::element_line(color = line_color),
      rect = ggplot2::element_rect(fill = NA, color = NA),
      text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_blank(), # modified by MH
      axis.text.y = ggplot2::element_blank(),# modified by MH
      axis.text.x = ggplot2::element_blank(),# modified by MH
      legend.key = ggplot2::element_rect(color = NA, fill = NA),
      #legend.key.height = ggplot2::unit(3, "mm"),
      panel.border = ggplot2::element_rect(color = line_color, fill = NA),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
}


#' @export
ggplot2_map_theme <- function(legend_xoffset = 0, legend_yoffset = 0) {
  stopifnot(requireNamespace("ggplot2"))

  ggplot2::theme(
    # remove axis title label
    axis.title = ggplot2::element_blank(),
    # rotate y-axis tick mark labels to be vertical
    axis.ticks = ggplot2::element_blank(), # modified by MH
    axis.text.y = ggplot2::element_blank(),# modified by MH
    axis.text.x = ggplot2::element_blank(),# modified by MH
    # remove area around a plot
    plot.margin = ggplot2::margin(0, 0, 0, 0),
    # set "tag" to bold text in top-left corner
    plot.tag = ggplot2::element_text(face = "bold", hjust = 0, vjust = 0),
    plot.tag.position = c(0.1, 0.95),
    # beautify legends: larger font, white background, lower-right corner
    legend.title = ggplot2::element_text(size = ggplot2::rel(1)),
    legend.text = ggplot2::element_text(size = ggplot2::rel(1)),
    legend.position = c(0.8 + legend_xoffset, 0.4 + legend_yoffset), # lower-right corner of map
    legend.justification = c("left", "top"),
    legend.background = ggplot2::element_rect(fill = "white", linewidth = 0),
    legend.margin = ggplot2::margin(5, 5, 5, 5)
  )
}

#------ Plots ------


#' Create stacked bar plot for spatial tabulations
#' @export
plot_spatialtabulation_stackedbars <- function(data, var_y, colorscale) {
  stopifnot(
    requireNamespace("ggplot2"),
    requireNamespace("rlang"),
    c("values", var_y, "Response", "Var") %in% colnames(data)
  )

  ggplot2::ggplot(data = data) +
    ggplot2::aes(
      x = values,
      y = !!rlang::sym(var_y),
      fill = Response
    ) +
    ggplot2::facet_wrap(ggplot2::vars(Var)) +
    ggplot2::geom_col(position = ggplot2::position_fill(reverse = TRUE)) +
    colorscale(levels = levels(data[["Response"]])) +
    ggplot2::xlab("Proportion of Area") +
    ggplot2::ylab(NULL) +
    newRR3::ggplot2_clean_theme() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      #axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 1),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text.x.top = ggplot2::element_text(size = 11), # base_size = 11
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line(color = "gray")
    )
}

#' Create stacked bar plot for spatial tabulations
#' @export
plot_spatialtabulation_stackedbars2 <- function(data, var_y, colorscale) {
  stopifnot(
    requireNamespace("ggplot2"),
    requireNamespace("rlang"),
    c("values", var_y, "Response", "Var") %in% colnames(data)
  )

  ggplot2::ggplot(data = data) +
    ggplot2::aes(
      x = values,
      y = !!rlang::sym(var_y),
      fill = Response
    ) +
    ggplot2::facet_wrap(ggplot2::vars(Var)) +
    ggplot2::geom_bar(
      position = ggplot2::position_stack(reverse = TRUE),
      stat = "identity"
    ) +
    colorscale(levels = levels(data[["Response"]])) +
    ggplot2::xlab("Proportion of Area") +
    ggplot2::ylab(NULL) +
    newRR3::ggplot2_clean_theme() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      #axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 1),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text.x.top = ggplot2::element_text(size = 11), # base_size = 11
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line(color = "gray")
    )
}


#' Create stacked bar plot for spatial tabulations, grouped for low, med, high
#'
#' @references for variable-width stacked bars
#' @references \url{https://stackoverflow.com/questions/5948604/variable-width-bars-in-ggplot2-barplot-in-r}
#' @references \url{https://stackoverflow.com/questions/9556359/different-width-for-each-bar-using-geom-bar-with-position-fill-in-ggplot2}
#'
#' @export
plot_spatialtabulation_stackedbars3 <- function(data, var_y, colorscale) {
  stopifnot(
    requireNamespace("ggplot2"),
    requireNamespace("rlang"),
    c("values", var_y, "Response", "Var") %in% colnames(data)
  )
  ggplot2::ggplot(data = data) +
    ggplot2::aes(
      x = values,
      y = factor(acrmod_level),
      fill = Response,
      width = width / max(width, na.rm = TRUE)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(!!rlang::sym(var_y)),
      cols = ggplot2::vars(Var)
    ) +
    ggplot2::geom_bar(
      position = ggplot2::position_fill(reverse = TRUE),
      stat = "identity"
    ) +
    colorscale(levels = levels(data[["Response"]])) +
    ggplot2::xlab("Proportion of Area") +
    ggplot2::ylab(NULL) +
    newRR3::ggplot2_clean_theme() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      #axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = ggplot2::element_text(hjust = 1),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text.x.top = ggplot2::element_text(size = 11), # base_size = 11
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line(color = "gray")
    )
}

#' @export
add_coords <- function(bbox, crs, expand_bbox = TRUE) {
  ggplot2::coord_sf(
    xlim = bbox[c("xmin", "xmax")],
    ylim = bbox[c("ymin", "ymax")],
    crs = crs,
    expand = expand_bbox
  )
}

#' Create one map
#' @export
plot_map <- function(
  xstars,
  add_coords = TRUE,
  crs = sf::st_crs(xstars),
  bbox = sf::st_bbox(xstars),
  expand_bbox = TRUE,
  maintitle = NULL,
  subtitle = NULL,
  panel_tag = NULL,
  show_legend = TRUE,
  legend_xoffset = 0.,
  legend_yoffset = 0.,
  st_poly_hatched = NULL,
  st_geom_sb = NULL,
  st_geom_state = NULL,
  st_geom_eco = NULL,
  ... # passed to geom_stars
) {
  stopifnot(
    requireNamespace("ggplot2"),
    requireNamespace("stars")
  )


  #--- * Plot raster & title ------
  tmp <- ggplot2::ggplot() +
    stars::geom_stars(
      data = xstars,
      show.legend = show_legend,
      ...
    )

  if (!is.null(maintitle)) {
    tmp <- tmp +
      ggplot2::labs(
        title = paste0(
          maintitle,
          if (!is.null(subtitle)) {
            paste0("\n(", subtitle, ")")
          }
        )
      )
  }

  if (!is.null(panel_tag)) {
    tmp <- tmp +
      ggplot2::labs(tag = panel_tag)
  }


  #--- * Beautify plot ------
  if (!is.null(st_geom_sb)) {
    tmp <- tmp +
      # cover up raster areas outside of sagebrush biome polygon
      ggplot2::geom_sf(
        data = create_negmask(poly_not_erase = st_geom_sb, bbox = bbox),
        fill = "white",
        size = 0,
        color = "white"
      )
  }

  if (!is.null(st_geom_state)) {
    tmp <- tmp +
      ggplot2::geom_sf(
        data = st_geom_state,
        size = 0.25,
        linetype = "dashed",
        fill = NA
      )
  }

  if (!is.null(st_geom_state)) {
    tmp <- tmp +
      ggplot2::geom_sf(data = st_geom_eco, size = 0.5, fill = NA)
  }

  if (!is.null(st_geom_sb)) {
    tmp <- tmp +
      # add outline of sagebrush biome polygon
      ggplot2::geom_sf(data = st_geom_sb, size = 1, fill = NA)
  }

  if (!is.null(st_poly_hatched)) {
    stopifnot(requireNamespace("ggpattern"))
    tmp <- tmp +
      # add hatched, transparent polygon (e.g., non-robust signal)
      ggpattern::geom_sf_pattern(
        data = st_poly_hatched,
        pattern = "stripe",
        color = NA,
        fill = NA,
        pattern_angle = 45,
        pattern_density = 0.01,
        pattern_spacing = 0.015,
        pattern_size = 0.3,
        pattern_res = 600,
        pattern_colour = "black",
        pattern_fill = "black"
      )
  }

  if (add_coords) {
    tmp <- tmp +
      add_coords(bbox = bbox, crs = crs, expand_bbox = expand_bbox)
  }

  tmp +
    ggplot2_clean_theme()
}


#' Create minimal plot of weighted (count * density)
#' @export
inset_densitycountplot <- function(
  x,
  limits = c(NA, NA),
  x_binned = FALSE,
  weight = 1,
  consistency = NULL,
  consistency_lvls = NULL,
  add_vertical0 = TRUE,
  show.legend = FALSE
) {

  if (is.null(consistency)) {
    dflt_consistency <- if (is.null(consistency_lvls)) {
      "Overall"
    } else {
      consistency_lvls[[1L]]
    }
    consistency_lvls <- dflt_consistency
  }

  tmp_inset <- ggplot2::remove_missing(
    df = data.frame(
      x = x,
      w = if (is.null(weight)) 1L else weight,
      grp = factor(
        if (is.null(consistency)) dflt_consistency else consistency,
        levels = consistency_lvls
      )
    ),
    finite = TRUE,
    na.rm = TRUE # remove missing/non-finite values silently
  )


  # remove group level with less than two elements
  ttmp <- table(tmp_inset[["grp"]])
  if (any(ttmp <= 2)) {
    id <- names(ttmp)[which(ttmp <= 2)]
    if (length(id) == 1L) {
      tmp_inset <- tmp_inset[!(tmp_inset[["grp"]] %in% id), , drop = FALSE]
    }
  }


  # create plot if enough data
  if (nrow(tmp_inset) > 2) {
    coi <- colors_okabeito()

    res <- ggplot2::ggplot(tmp_inset) +
      ggplot2::aes(
        x = x,
        # `after_stat(count)`:
        #   * y for `geom_density`: "density * number of points"
        #   * y for `geom_bar` (default): "number of cases"
        y = ggplot2::after_stat(count),
        weight = w,
        fill = grp
      ) +
      ggplot2::scale_fill_discrete(
        type = if (length(consistency_lvls) > 1) {
          c(coi[[2L]], coi[-2L]) # other: blue - robust:orange
        } else {
          coi[[3L]] # green
        },
        limits = factor(consistency_lvls)
      )

    res <- if (isTRUE(x_binned)) {
      res +
        ggplot2::geom_bar(
          alpha = 0.5,
          color = NA,
          position = "identity", # don't stack bars
          show.legend = show.legend,
          na.rm = TRUE # remove missing values silently
        ) +
        ggplot2::scale_x_binned()
    } else {
      res +
        ggplot2::geom_density(
          alpha = 0.5,
          color = NA,
          n = 1024, # default 512 produced visible edges for some novelty metrics
          show.legend = show.legend,
          na.rm = TRUE # remove missing values silently
        )
    }

    if (isTRUE(add_vertical0)) {
      res <- res +
        # add vertical dashed line at zero
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
    }

    res +
      # zoom (and not clip) into `limits`
      ggplot2::coord_cartesian(xlim = limits) +
      # beautify plot theme
      ggplot2_clean_theme() +
      ggplot2::theme(
        # remove most area around a plot (but leave small space for text at x-axis ticks)
        plot.margin = ggplot2::margin(0, 5.5, 0, 5.5),
        # make plot background white
        plot.background = ggplot2::element_rect(fill = "white", color = "gray25"),
        panel.border = ggplot2::element_rect(color = NA, fill = NA),
        # remove x-axis title
        axis.title.x = ggplot2::element_blank(),
        # slightly reduce size of x-axis tick text
        axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.75)),
        # remove y-axis
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        # legend: top-left corner
        legend.position = c(0.05, 0.95),
        legend.justification = c("left", "top"),
        legend.background = ggplot2::element_rect(fill = "white", linewidth = 0)
      )
  }
}
