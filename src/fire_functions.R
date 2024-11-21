# Contains functions for predicting fire probability


# fire probability, predicted by equation in Holdrege et al. (2024)
predict_fire <- function(mat, map, psp, afg, pfg, run_checks = TRUE) {
  # these checks don't work if input is a raster,
  # so set to false in that case
  if(run_checks) {
    if(any(mat > 70 | mat < - 70)) {
      stop('mat should be in C')
    }
    if (any(psp<0 | psp>1)) {
      stop('psp has a range of 0-1')
    }
  }

  mat <-  mat + 273.15
  y = -92.86 + 14.01*log10(afg + 1) - 0.8117*(log10(afg + 1))^2 +
    0.5054*pfg^0.5 - 0.03730*pfg +0.02672*mat + 49.25*log10(map + 1) - 
    8.236*log10(map + 1)^2 -7.505*log10(psp + 0.001) - 3.118*log10(psp + 0.001)^2 -
    4.047*log10(afg + 1)*log10(map + 1)
  
  p <- 1/(1 + exp(-y))
  p
}

# testing
if (FALSE) {
  afg <- seq(0, 120, by = 1)
  p <- predict_fire(mat = 283, map = 250, psp = 0.2, afg = afg,
                    pfg = 30)
  plot(afg, p*100) 
}


# (approximate) derivative of the fire equation, with respect to one
# predictor
# note b/ the logistic equation part of the fire equation, the derivative
# does depend on the value of the other variables
# deriv_var is a string (name of variable), the other inputs are vectors of 
# the predictor variables
calc_fire_deriv <- function(deriv_var, mat, map, psp, afg, pfg, ...)  {
  l1 <- list(mat = mat, map = map, psp = psp, afg = afg, pfg = pfg)
  l2 <- l1
  stopifnot(deriv_var %in% names(l1))
  
  # determine width of the change
  tmp <- l1[[deriv_var]]
  range <- if(is.numeric(tmp)) {
    range(tmp)
  } else if ("SpatRaster" %in% class(tmp)){
    range(values(tmp), na.rm = TRUE)
  } else {
    stop('incorrect data type')
  }

  # makes +/- change be 1% of the range of the variable
  change <- (range[2] - range[1])*0.005 
  l1[[deriv_var]] <- l1[[deriv_var]] - change
  l2[[deriv_var]] <- l2[[deriv_var]] + change
  first <- predict_fire(mat = l1$mat, map = l1$map, psp = l1$psp, afg = l1$afg,
                        pfg = l1$pfg,
                        ...)
  second <- predict_fire(mat = l2$mat, map = l2$map, psp = l2$psp, afg = l2$afg, 
                         pfg = l2$pfg,
                         ...)
  
  # delta per 1% change in pfg
  out <- second - first
  names(out) <- deriv_var
  out
}


#' wrapper to use for raster inputs
#'
#' @param f function that takes the 5 inputs used by the fire equation
#' @param r raster with properly named layers
fire_rast_wrapper <- function(f, r, ...) {
  f(mat = r[['mat']], map = r[['map']], psp = r[['psp']], afg = r[['afg']], 
    pfg = r[['pfg']], ...)
}
