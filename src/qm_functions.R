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
qm_from_quantiles <- function(x, from, to) {
  
  stopifnot(length(from) == length(to),
            # make sure from/to are in ascending order
            all(diff(from) >= 0),
            all(diff(to) >= 0))
  
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
    out[low] <- min(to)
  }
  
  if(length(high) > 0) {
    out[high] <- max(to)
  }

  out
}




