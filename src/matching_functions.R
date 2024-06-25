library(rMultivariateMatching)
library(dplyr)

# Next steps
# use stats::kmeans to find k centroids
# find the nearest gridmet cells to the centroids
# find the coverage of 


  
#' Bin values into n equal sized bins
#' 
#' @description breaks x into n bins, where the bins by default
#' are 10 equal spaced categories between the 1st and 99th percentile of
#' x (values below 1st percentile and above 99th percentile get put in the
#' first and last bin, respectivel)
#'
#' @param x numeric vector
#' @param xlab string--label to be pasted into names of the bins
#' @param n_categories number of bins to create
#'
#' @return vector
#' @examples
#' x <- rnorm(20)
#' cut_bins(x)
cut_bins <- function(x, xlab = 'x', n_categories = 10) {
  lower <- quantile(x, 0.01)
  upper <- quantile(x, 0.99)
  
  break_points <- seq(lower, upper, length.out = n_categories + 1)
  
  # making the first ad last cut points infinity so that 
  # all data at the tails is included
  break_points[1] <- -Inf
  break_points[length(break_points)] <- Inf
  
  labels <- paste0(xlab, "-", 1:n_categories)
  
  out <- cut(x, breaks = break_points, labels = labels)
  out <- as.character(out)
  out
}
  
#' Bin columns of a dataframe
#'
#' @param df dataframe
#' @param cols columns to bin
#' @param n_categories number of groups to break each column into
#'
#' @return dataframe with new column bin_all column which is the
#' concatenation of all the bins
#' @examples
#' n = 1000
#' df <- tibble(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = runif(n),
#'   x4 = 3*x3 + runif(n, 0, 0.1),
#'   x5 = runif(n, 100, 200)
#' )
#' bin_df(df, names(df))
bin_df <- function(df, cols, n_categories = 10) {
  new_cols <- paste0('bin_', cols)
  
  out <- df
  for (i in seq_along(cols)) {
    col <- cols[i]
    new <- new_cols[i]
    out[[new]] <- cut_bins(out[[col]],
                           xlab = col,
                           n_categories = n_categories)
  }
  out <- tidyr::unite(out, col = 'bin_all', all_of(new_cols), 
               sep = "_")
  out
}  

#' finding evenly spaced centers for starting kmeans
#'
#' @param df dataframe
#' @param bin_vars variables to use for binning the dataset
#' @param n number of centers
#' @param n_categories number of equal interval groups to break in bin_vars into
starting_centers <- function(df, bin_vars, n, n_categories = 1) {
  
  binned <- bin_df(df, cols = bin_vars, n_categories = n_categories)
  n_possible <- length(unique(binned$bin_all))
  
  while(n > n_possible) {
    n_categories <-  n_categories + 1
    binned <- bin_df(df, cols = bin_vars, n_categories = n_categories)
    
    n_possible <- length(unique(binned$bin_all))
    print(n_possible)
  }
  

  centers <- group_by(binned, bin_all) %>% 
    slice_sample(n = 1) %>% 
    ungroup() %>% 
    slice_sample(n = n) %>% 
    select(-bin_all)
    
  centers

}

#' fitting kmeans in situations where convergence is difficult
#'
#' @param df data frame, which only columns to be used for kmeans
#' @param k number of clusters
#' @param cols_for_centers column names of variables for selecting initial
#' centers (samples from an even grid)
#' @param prop_seq vector of proportion of dataset to fit kmeans to, first fit 
#' to a small fraction of the dataset, then use those centers to fit to a larger
#' fraction of the data set and so on
#' @param n_categories number of categories to break each cols_for_centers 
#' into for initial centers
#' @param iter.max max iterations for kmeans
#' @param verbose 
#'
#' @return
#' kmeans ojbection
#' @export
#' 
#' @examples
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Restrict data to matching variables of interest
#' df <- allvars[,c("bioclim_01","bioclim_04",
#'    "bioclim_09","bioclim_12","bioclim_15","bioclim_18")]
#'    cols_for_centers = c("bioclim_01", "bioclim_12", "bioclim_18")
#'  criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'  for (i in 1:ncol(df)) {
#'    df[, i] <- df[, i]/criteria[i]
#'  }
fit_kmeans <- function(df, k, 
                       cols_for_centers = names(df),
                       prop_seq = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.75,1),
                       random_centers = FALSE,
                       digits = NULL,
                       n_categories = 5, 
                       iter.max = 100,
                       verbose = TRUE
) {

  if(!is.null(digits)) {
    df <- purrr::map_dfc(df, \(x) signif(x, digits = digits))
  }
  # set.seed(1234) # for testing

  nr <- nrow(df)
  rows <- 1:nr
  rows_sample <- sample(rows, size = round(prop_seq[1]*nr))
  
  # selecting evenly spaced starting points, to improve convergence
  
  if(random_centers) {
    centers <- k
  } else {
    centers <- starting_centers(df[rows_sample, ],
                                bin_vars = cols_for_centers, 
                                n = k, 
                                n_categories = n_categories)
  }
  algorithm <- "Hartigan-Wong" #
  kmns <- stats::kmeans(df[rows_sample, ], centers = centers,
                        algorithm = algorithm,
                        iter.max = iter.max,
                        nstart = 1)
  # options(warn = 1) # for debugging
  for(i in 2:length(prop_seq)) {
    if(prop_seq[i] < 1) {
      size_new <- round((prop_seq[i] - prop_seq[i - 1]) *nr)
      # adding new rows to the previous sample
      # this way avoid the problem of potential centers with zero points nearby
      rows_sample <- c(rows_sample, sample(rows[-rows_sample], size = size_new))
    } else {
      rows_sample <- rows
    }
    
    if(verbose) {
      print(paste('proportion of data used:', prop_seq[i]))
    }
    kmns <- stats::kmeans(df[rows_sample, ], 
                          algorithm = algorithm,
                          centers = kmns$centers,
                          iter.max = iter.max,
                          nstart = 1) 
    #if(kmns$ifault==4) warning('convergence issues')
  }
  
  kmns
}
  
#' select observations closest to kmeans centers
#'
#' @param df dataframe
#' @param kmns kmeans object
#' @param criteria_vars names of variables (columns for calculating distances)
#'
#' @return
#' subset df, these are the rows that are nearest neighbors to the centers in 
#' kmns object, also column 'distance' added (i.e. the select 'sites')
select_nearest_neighbors <- function(df, kmns, criteria_vars) {
  centers <- kmns$centers
  stopifnot(criteria_vars %in% colnames(centers),
            criteria_vars %in% names(df))
  
  dat <- as.matrix(df[criteria_vars])
  query <- as.matrix(centers[, criteria_vars])
  nn <- FNN::get.knnx(data = dat, 
                      query = query, k = 1)
  index <- nn$nn.index
  
  # removing duplicates (1 observation might be closest to two centers)
  not_dup <- !duplicated(index)

  out <- df[index[not_dup], ]
  # distance of observation to kmean center
  out$distance <- as.numeric(nn$nn.dist)[not_dup] 
  out
 
}


#' Title
#'
#' @param df dataframe of 'target cells'
#' @param sites selected sites (output from select_nearest_neighbors)
#' @param criteria_vars variables (column names) used for calculating distance
#'
#' @return
#' proportion of rows in dataframe that have a euclidean distance < cuttoff (default)
#' from the sites
#' @examples\
#' # Load targetcells data for Target Cells
#' data(targetcells)
#'
#' # Create data frame of potential matching variables for Target Cells
#' allvars <- makeInputdata(targetcells)
#'
#' # Restrict data to matching variables of interest
#' criteria_vars <- c("bioclim_01","bioclim_04",
#'    "bioclim_09","bioclim_12","bioclim_15","bioclim_18")
#' df <- allvars[,criteria_vars]
#' criteria <- c(0.7,42,3.3,66,5.4,18.4)
#'  for (i in 1:ncol(df)) {
#'    df[, i] <- df[, i]/criteria[i]
#'  }
#'  kmns <- fit_kmeans(df, k = 400,
#'     cols_for_centers = c("bioclim_01", "bioclim_12", "bioclim_18"), digits = 4)
#'  # note--for select_nearest_neighbors df argument could be from a different dataset
#'  # e.g. gridmet, where kmeans was fit to daymet. 
#'  sites <- select_nearest_neighbors(df = df, kmns = kmns,
#'      criteria_vars = criteria_vars)
#'  prop_coverage(df, sites = sites, criteria_vars = criteria_vars)
prop_coverage <- function(df, sites, criteria_vars,
                          cuttoff = 1) {
  dat <- as.matrix(sites[criteria_vars])
  query <- as.matrix(df[criteria_vars])
  
  nn <- FNN::get.knnx(data = dat, 
                  query = query, k = 1)
  
  out <- mean(as.numeric(nn$nn.dist) < cuttoff)
  out
}

kpoints2 <- function(target, potential_sites, k,
                     criteria_vars,
                     cols_for_centers = criteria_vars,
                     prop_seq = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.75,1),
                     random_centers = FALSE,
                     digits = 4,
                     n_categories = 5, 
                     iter.max = 100,
                     verbose = FALSE,
                     cuttoff = 1) {
  # find 'center's via kmeans 
  kmns <- fit_kmeans(df = target, k = k,
                      cols_for_centers = cols_for_centers, 
                      prop_seq = prop_seq,
                      random_centers = random_centers,
                      digits = digits,
                      n_categories =  n_categories, 
                      iter.max = iter.max,
                      verbose = verbose)
   
   sites <- select_nearest_neighbors(df = potential_sites, kmns = kmns,
       criteria_vars = criteria_vars)
   prop <- prop_coverage(df = target, sites = sites, criteria_vars = criteria_vars)
   out <- list(coverage = coverage,
               sites = sites,
               ifault = kmns$ifault)
   out
}


 