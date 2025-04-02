# Functions used for interpolation of STEPWAT2 results

#' which columns to upscale (based on what has already be done)
#'
#' @param df dataframe with columns to upscale
#' @param path path to folder where already upscale .tifs live
#' @param pattern for filtering which files to look at in path
#' @param min_size minimum file size to consider that file to have already
#' been succesfully upscaled
#' @param rerun should files be upscaled again even if the files were already
#' created?
#'
#' @return vector of column names
which_todo <- function(df, path, pattern, min_size, rerun = FALSE) {
  # what have already been upscaled?
  already_done1 <- list.files(path, pattern = pattern, full.names = TRUE)  
  
  already_done2 <- already_done1[file.size(already_done1) > min_size] %>% 
    basename()
  already_done2 <- str_replace(already_done2, ".tif", "")
  
  todo1 <- names(df)
  
  # recreate existing files?
  if(rerun) {
    todo2 <- todo1
  } else {
    todo2 <- todo1[!todo1 %in% already_done2]
  }
  
  todo2
}


#' create list of column numbers to loop through with dopar (when have
#' very wide dataframe)
#'
#' @param df dataframe of columns where want upscale all but the 
#' first column
#' @param by how many columns should get passed to 
#' rMultivariateMatching::interpolatePoints() at once
#'
#' @return list where each element consists of a vector of column numbers
col_nums_parallel <- function(df, by = 50) {
  # making list of columns to use each time through dopar loop
  stopifnot(is.data.frame(df))
  n <- ncol(df)
  
  vecs <- seq(from = 2, to = n, by = by)
  vecs_l <- map(vecs, function(x) {
    to <- x + by - 1
    
    # for the last set of columns
    if(to > n) {
      to <- n
    }
    x:to # column numbers
  })
  
  # check that all columns accounted for
  stopifnot(unlist(vecs_l) == 2:n)
  vecs_l
}


#' run the interpolatePoints function
#'
#' @param df dataframe passed to interpolatePoints
#' @param match dataframe that is output by multivarmatch
#' @param template raster template 
#' @param path path of where .tif files will be put
#' @param min_size minimum file size to consider an existing .tif not corrupt
#' (applies when rerun is FALSe)
#' @param rerun recreate tif files that already exist?
#' @param exclude_poor_matches limit interpolation to well matched areas?
#' # Run the following before using this function:
#' num.cores <- parallel::detectCores() 
#' registerDoParallel(num.cores)
#' @returns
#' doesn't return anything saves rasters to tif files
run_interpolatePoints <- function(df, 
                                  match, 
                                  template,
                                  path, 
                                  min_size, 
                                  rerun,
                                  exclude_poor_matches = TRUE) {
  if(!dir.exists(path)) {
    dir.create(path)
  }
  
  todo2 <- which_todo(df = df,
                      path = path,
                      pattern = ".tif$",
                      min_size = min_size,
                      rerun = rerun)
  
  # if all have been upscaled then todo just includes 'cellnumbers'
  if(length(todo2) > 1) {
    
    if(test_run) {
      todo2 <- todo2[1:3]
    }
  
    df2 <- df[, todo2]
    
    # making list of columns to use each time through dopar loop
    vecs_l <- col_nums_parallel(df2, by = 20)
    
    message(paste(path, 'start'))
    message(Sys.time())
    
    # separately running interpolation on different sets of columns
    foreach (x = vecs_l) %dopar% {
      rMultivariateMatching::interpolatePoints(
        matches = match,
        output_results = df2[, c(1, x)], 
        exclude_poor_matches = exclude_poor_matches,
        subset_cell_names = "subset_cell",
        quality_name = "matching_quality",
        matching_distance = 1.5,
        raster_template = template,
        plotraster = FALSE,
        saveraster = TRUE,
        filepath = path,
        overwrite = TRUE
      )
    }
    message(paste(path, 'done'))
    message(Sys.time())
  }

}
