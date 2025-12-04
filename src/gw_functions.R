# Functions for calculating GCM-wise (as opposed to pixelwise) 
# summaries of GCMs (e.g. results for the GCM with median of a variable,
# as opposed to using median (across gcms) at each pixel fora particular variable)



#' Return site (pixel) level data for lo, median and high GCMs
#'
#' @param df dataframe
#' @param summarize_vars vector of variables interested in
#' @param group_vars grouping variables
#'
#' @returns original data, except just for the low, median and high GCMs
gw_select_smry_site <- function(df, summarize_vars, 
                                group_vars = c('years', 'RCP', 'graze', 
                                               'PFT', 'region')) {
  vars <- summarize_vars
  gcm_smry <- df %>% 
    group_by(across(all_of(group_vars)), .data$GCM) %>% 
    summarise(across(all_of(vars), 
                     .fns = \(x) weighted_median(x, weight = weight)),
              .groups = 'drop')
  dfs_smry <- map(vars, function(var) {
    
    df2 <- df[c(group_vars, 'site', 'weight', 'GCM', var)]
    long <- gcm_smry %>% 
      group_by(across(all_of(group_vars))) %>% 
      summarise(low = calc_low(.data[[var]]),
                high = calc_high(.data[[var]]),
                median = median(.data[[var]]),
                .groups = 'drop') %>% 
      pivot_longer(cols = c('low', 'high', 'median'),
                   names_to = 'summary')
    
    #stopifnot(!is.na(long$value) | (is.na(long$value) & long$RCP == 'Current'))
    long2 <- long %>% 
      filter(!is.na(value))
    
    long3 <- long2 %>% 
      left_join(gcm_smry[c(group_vars, 'GCM', var)],
                by = c(group_vars, 'value' = unname(var))) %>% 
      mutate(summary = summary2factor(summary))
    
    test <- long3 %>% 
      group_by(across(c(all_of(group_vars), 'summary'))) %>% 
      summarise(n = n(),
                .groups = 'drop') %>% 
      pull(n) %>% 
      unique()
    
    # if fails would need to update function to select one of 
    # multiple 'equivalent' GCMs
    stopifnot(lu(test) == 1)
    
    long3 %>% 
      select(-value) %>% 
      left_join(df2, by = c(group_vars, 'GCM'))
    
  })
  
  if(length(dfs_smry) == 1) {
    out <- dfs_smry[[1]]
  } else {
    out <- reduce(dfs_smry, left_join)
  }
  out
  
}
