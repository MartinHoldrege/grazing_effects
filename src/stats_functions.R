# stats functions

# slope -> angle in degrees
slope_to_deg <- function(m) {
  out <- atan(m) * 180 / pi
  # keep NA as NA; handle +/-Inf explicitly (atan handles it too, but this is clearer)
  out[is.infinite(m)] <- sign(m[is.infinite(m)]) * 90
  out[is.na(m)] <- NA_real_
  out
}

calc_slope <- function(x, y) {
  mod <- lm(y ~ x)
  coefs <- coef(mod)
  stopifnot(length(coefs) ==2)
  slope <- coefs[[2]]
  slope
}

#' fit lm and convert slope to degrees
#'
#' @param x 
#' @param y 
#'
#' @returns
# slope of the model y~x in terms of degrees (-90 to 90)
calc_slope_deg <- function(x, y) {
  slope <- calc_slope(x, y)
  slope_to_deg(slope)
}
