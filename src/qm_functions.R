# Description: functions for quantile matching
# focusing on functions that do quantile matching (and linear interpolation)
# based on knowing a limited number of quantiles. The purpose
# is to keep it simple enough that an analogous function can be implemented
# in STEPWAT2 in C. 


#' conduct quantile matching
#'
#' @param x vector of values to create matched values for
#' @param from vector of values (in units of x), that provide values for a given set of quantiles
#' those quantiles should include 0 and 1 (i.e. the min and max)
#' @param to vector of values (in units of the desired output) that provides values for the same
#' set of quantiles
#'
#' @return vector with same length as x 
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 11
#' quant <- seq(0, 1,length.out = n)
#' from <- sort(abs(rnorm(n = n, mean = 4, sd = 2))) # values corresponding 
#' to <- from^4/200
#' qm_from_quantiles(c(1, 2, 5), from = from, to = to)
#' # testing
#' all(qm_from_quantiles(from, from, to) == to)
#' all((qm_from_quantiles(from + 0.001, from, to) - to) < 0.01)
qm_from_quantiles <- function(x, from, to, ascending = TRUE) {
  
  stopifnot(length(from) == length(to),
            # make sure from/to are in ascending order
            all(diff(from) >= 0))
  
  if (ascending) {
    stopifnot(all(diff(to) >=0))
  }
  
  
  low <- which(x <= min(from))
  high <- which(x>= max(from))
  
  
  x_subset <- if (length(c(low, high) > 0)) {
    to_subset <- (1:length(x))[-c(low, high)]
    x[to_subset] 
  } else {
    to_subset <- 1:length(x)
    x
  }
  
  out <- vector(mode = 'numeric', length = length(x))
  out[] <- NA
  if(length(x_subset) > 0) {
    start <- purrr::map_dbl(x_subset, function(z) max(which(z >= from)))
    
    # linear interpolation
    frac <- (x_subset - from[start])/(from[start + 1] -  from[start])
    
    y_delta <- to[start + 1] - to[start]
    
    y_new <- to[start] + y_delta*frac # quantile of the x value
    out[to_subset] <- y_new
  }
  
  # replace values lower and higher than any of the reference
  # breaks then replace with either the min or max
  if(length(low) > 0) {
    out[low] <- to[1]
  }
  
  if(length(high) > 0) {
    out[high] <- to[length(to)]
  }
  out
}

# same arguments as qm_from_quantiles, but ~20x faster
# this was written by chat GPT and tested and checked by me
qm_from_quantiles_fast <- function(x, from, to, ascending = TRUE) {
  stopifnot(length(from) == length(to),
            all(diff(from) >= 0))
  if (ascending) stopifnot(all(diff(to) >= 0))
  
  # Index of the lower bin
  idx <- findInterval(x, from, rightmost.closed = TRUE)
  
  # Clamp out-of-bounds
  idx[idx < 1] <- 1
  idx[idx >= length(from)] <- length(from) - 1
  
  x0 <- from[idx]
  x1 <- from[idx + 1]
  y0 <- to[idx]
  y1 <- to[idx + 1]
  
  # Linear interpolation
  frac <- (x - x0) / (x1 - x0)
  frac[!is.finite(frac)] <- 0  # handle divide-by-zero
  y <- y0 + frac * (y1 - y0)
  
  # Cap extremes
  y[x <= min(from)] <- to[1]
  y[x >= max(from)] <- to[length(to)]
  
  y
}
if(FALSE) {
  set.seed(123)
  x <- sort(runif(100, min = -1, max = 11))
  from <- c(0, 2, 4, 6, 8, 10)
  to <- c(0, 1, 2, 4, 8, 10)
  
  old_y <- qm_from_quantiles(x, from, to)
  new_y <- qm_from_quantiles_fast(x, from, to)
  
  all.equal(old_y, new_y, tolerance = 1e-8)
  
  library(microbenchmark)
  
  set.seed(123)
  x <- sort(runif(1e5, min = 0, max = 10))
  from <- c(0, 2, 4, 6, 8, 10)
  to <- c(0, 1, 2, 4, 8, 10)
  
  microbenchmark::microbenchmark(
    original = qm_from_quantiles(x, from, to),
    fast     = qm_from_quantiles_fast(x, from, to),
    times = 10
  )
}







qm_quant_factory <- function(from, to) {
  function(x) {
    qm_from_quantiles(x, from = from, to = to)
  }
}

# these are the numbers from the
# "data_processed/qm/quantiles_for_qm_0.5match_fire1_eind1_c4grass1_co20_2311_grazLMH.csv"
# file, same values used for quantile mapping in stepwat2
# that file was created in scripts/bio_matching/02_quantile_matching_simple.R
qm_Aherb <- qm_quant_factory(
  from = c(6.49, 14.48, 18.38, 23.77, 37.96, 45.45, 53.98, 63.03, 69.16, 
           76.48, 83.03, 87.28),
  to = c(0, 0.64, 2.9, 4.31, 6.11, 8.56, 12.21, 18.04, 29.35, 58.62, 
         82.68, 172.44))

qm_Pherb <- qm_quant_factory(
  from = c(13.71, 21.8, 41.19, 46.63, 51.77, 57.89, 67.89, 81.59, 98.22, 
           114.58, 123.15, 140.67),
  to = c(0, 6.52, 18.93, 25.66, 33.57, 42.89, 54.42, 69.71, 88.57, 116.52, 
         136.07, 309.82))

# values from the "data_processed/qm/quantiles_for_qm_sagebrush_0.5match_fire1_eind1_c4grass1_co20_2503_LMH.csv"
# calculated in the scripts/bio_matching/qm_sagebrush.R script
# quantile mapping from stepwat biomass to rcmap cover
qm_Sagebrush_bio2cov <- qm_quant_factory(
  from = c(94.9, 99.8, 104, 130, 165, 171, 171, 183, 186, 196, 200, 202, 
           244, 268, 286, 297, 306, 314, 320, 326, 330, 333, 337, 338, 344, 
           349, 352, 356, 359, 364, 369, 371, 373, 375, 378, 382, 384, 386, 
           392, 396, 398, 401, 406, 409, 413, 414, 421, 423, 433, 441, 449, 
           463, 492, 515, 585, 641),
  to = c(0, 0, 0, 0, 1, 1, 1.67, 2, 2, 2.67, 3, 3, 3.33, 3.67, 4, 4, 
         4.33, 4.67, 5, 5, 5, 5.33, 5.67, 6, 6, 6.33, 6.33, 6.67, 7, 7, 
         7.33, 7.67, 8, 8.33, 8.67, 9, 9.33, 9.67, 10, 10.3, 10.7, 11, 
         11.7, 12, 12.3, 13, 13.3, 14, 14.7, 15.3, 16, 17.3, 20, 21.3, 
         25.7, 33.7))














