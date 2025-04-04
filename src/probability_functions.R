# functions related to probability calculations 


#' probability have being between a and b years post fire
#' 
#' @description
#' This function is meant to answer the question, given a particular
#' fire probability, what is the probability that a observation (location in space or time)
#' has not had fire for between a and b years
#' 
#'
#' @param p probability of fire in a given year, a vector, or a spatraster
#' @param a start of the period of interest (inclusive)
#' @param b end of the period of interest
prob_geometric_rule <- function(p, a = 0, b) {
  stopifnot(a <= b,
            a >= 0,
            any(class(p) %in% c('numeric', 'SpatRaster')))
  
  # making sure p could be a probability
  if(class(p) == 'numeric') {
    stopifnot(p>= 0 & p<=1)
  } else {
    mm <- as.numeric(minmax(p))
    stopifnot(mm[1] >= 0 & mm[2] <=1)
  }

  Pnf <- 1 - p # probability of no fire in a given year
  
  # probability of having a sequence of between a and b years of no fire
  # formula based on the 2nd proof provided here:
  # https://math.stackexchange.com/questions/971761/calculating-sum-of-consecutive-powers-of-a-number
  # note that instead of summing from 0 to n (as in that example),
  # I'm want to be able to sum from some number (a) to b, so the 
  # 'previous' (a - 1) sums are subtracted
  prob_a_b <- (1 - Pnf^(b + 1))/(1 - Pnf) - (1 - Pnf^a)/(1 - Pnf)
  
  prob_a_b*p # the sequence must occur after a fire
}


# helper funs -------------------------------------------------------------

create_age_groups <- function() {
  # for area by age class group
  # these intervals of years
  # to answer the question of (e.g) what is the expected area of the region
  # that will have burned between 0 and 10 years ago
  age_groups <- list(
    c(a = 0, b = 10),
    c(a = 11, b = 40),
    c(a = 41, b = 100),
    c(a = 101, b = Inf))
  
  names(age_groups) <- map_chr(age_groups, function(x) {
    out <- if(is.finite(x['b'])) {
      paste0(x['a'], '-',  x['b'])
    } else {
      paste0(x['a'], ' +')
    }
    out
  })
  
  age_groups
  
}





# testing -----------------------------------------------------------------


if (FALSE) {
  n <- 1e6
  
  p1 <- 0.01
  ps <- rbeta(n = n, shape1 = 0.01, shape2 = 0.99)
  mean(ps)
  
  # showing that even when each probability is different (but has same mean)
  # the lengths of no fire periods are very close to
  # those when the probability is fixed at the mean)
  x <- runif(n = n)
  hasNoFire1 <- x > p1
  res1 <- rle(hasNoFire1)
  # see below, this doesn't deal with 0 length sequences, so is approximate
  mean(res1[["lengths"]][res1[["values"]]] > 100) 
  hist(res1[["lengths"]][res1[["values"]]])
  mean(res1[["lengths"]][res1[["values"]]])
  
  mean(hasNoFire1)
  
  hasNoFire2 <- x > ps
  #hasNoFire2 <- rbinom(n = n, size = 1, prob = 1 - ps)
  hist(ps); summary(ps)
  res2 <- rle(hasNoFire2)
  hist(res2[["lengths"]][res2[["values"]]])
  mean(res2[["lengths"]][res2[["values"]]])
  # see below, this doesn't deal with 0 length sequences, so is approximate
  mean(res2[["lengths"]][res2[["values"]]] > 100) # fraction of time/space 
  mean(hasNoFire2)
  
  # showing that the function recreates the simulated values
  
  k <- 100
  lengths100 <- vector(length = k)
  lengths26_50 <- vector(length = k)
  for (i in 1:k) {
    x <- runif(n = n)
    hasNoFire1 <- x > p1
    res1 <- rle(hasNoFire1)
    lts <- res1[["lengths"]]
    lts[(lts == 1 & !res1[["values"]])] <- NA # single year that burned
    
    # this is year that burned, when the previous year burned
    # so this could be considered a unburned sequence of length 0.
    lts[(lts > 1 & !res1[["values"]])] <- 0
    lts <- lts[!is.na(lts)]
    
    lengths100[i] <- mean(lts > 100)
    # probability of a sequence having had no fire for 26-50 years
    lengths26_50[i] <- mean(lts >= 26 & lts <= 50)
  }
  
  # probability of being in a place/time where fire has not occurred for
  # over 100 years
  p100 <- 1 - prob_geometric_rule(p1, a = 0, b = 100)
  hist(lengths100)
  abline(v = p100, col = 'blue') # theoretical mean
  abline(v = mean(lengths100), col = 'red', lty = 2) # simulated mean
  
  # probability of being in a place/time where fire has has not occurred
  # in 26-50 years
  p26_50 <- prob_geometric_rule(p1, a = 26, b = 50)
  hist(lengths26_50); 
  abline(v = p26_50, col = 'blue') # theoretical mean
  abline(v = mean(lengths26_50), col = 'red', lty = 2) # simulated mean
}


