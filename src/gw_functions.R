# Functions for calculating GCM-wise (as opposed to pixelwise) 
# summaries of GCMs (e.g. results for the GCM with median of a variable,
# as opposed to using median (across gcms) at each pixel fora particular variable)



#' Return site (pixel) level data for lo, median and high GCMs
#'
#' @param df dataframe
#' @param summarize_vars vector of variables interested in
#' @param group_vars grouping variables
#' @param smry_fun function that takes two arguments: x and weight vectors
#'
#' @returns original data, except just for the low, median and high GCMs
gw_select_smry_site <- function(df, summarize_vars, 
                                group_vars = c('years', 'RCP', 'graze', 
                                               'PFT', 'region'),
                                smry_fun = weighted_median) {
  vars <- summarize_vars
  gcm_smry <- df %>% 
    group_by(across(all_of(group_vars)), .data$GCM) %>% 
    summarise(across(all_of(vars), 
                     .fns = \(x) smry_fun(x, weight)),
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

c3_to_long <- function(sei_byGCM) {
  c3_gcm2 <- sei_byGCM %>% 
    mutate(percent_ora = 100 - percent_csa - percent_goa) %>% 
    select(-matches('SEI_'), -group, -type, -id) %>% 
    pivot_longer(matches('percent_'), 
                 values_to = 'c3_percent',
                 names_to = 'c3') %>% 
    mutate(c3 = str_replace(c3, 'percent_', ''),
           c3 = str_to_upper(c3),
           c3_area = area*c3_percent/100) %>% 
    select(-area) %>% 
    df_factor()
  c3_gcm2
}


# this function creates the file that
# is the same as 'sei_pcent1' object previously
# used but with GCM-wise calculations
summarize_c3_gw <- function(sei_byGCM) {
  c3_gcm2 <- c3_to_long(sei_byGCM)
  
  join_vars <- c('region', 'RCP', 'years', 'graze')
  
  # which GCM has low, median or high values of 'mean' SEI
  # using mean SEI to sort gcms instead of median, because
  # because in some regions and scenarios median SEI is at or near
  # 0 so doesn't give seperation between GCMs
  sei_med_smry <- sei_byGCM %>% 
    group_by(across(all_of(join_vars))) %>% 
    summarize_across_GCMs(var = 'SEI_mean') %>% 
    left_join(sei_byGCM[c(join_vars, 'GCM', 'SEI_mean')],
              by = c(join_vars, 'SEI_mean')) %>% 
    filter(!is.na(GCM))
  
  test <- sei_med_smry %>% 
    group_by(across(all_of(join_vars)), summary) %>% 
    summarise(n = n(), .groups = 'drop')
  
  stopifnot(all(test$n == 1))
  
  # gcm-wise c3 class area
  sei_pcent1 <- c3_gcm2 %>% 
    right_join(sei_med_smry, by = c(join_vars, 'GCM')) %>% 
    mutate(rcp_year = rcp_label(RCP, years, include_parenth = FALSE,
                                add_newline = TRUE)) %>% 
    df_factor()
  
  # check if that got a 'unique' gcm for low, median and high
  test <- sei_pcent1 %>% 
    group_by(across(all_of(join_vars)), summary, c3) %>% 
    mutate(n = n(), .groups = 'drop') %>% 
    pull(n)
  
  stopifnot(all(test == 1))
  
  sei_pcent1
}


# summarize burned area
summarize_ba_gw <- function(df_ba_gcm, df_area) {
  df_area <- rename(df_area, area_total = area)
  df_ba_gcm2 <- left_join(df_ba_gcm, df_area,
                          by = 'ecoregion')
  
  df_ba_gcm3 <- df_ba_gcm2 %>% 
    group_by(ecoregion, RCP, years, graze, area_total) %>% 
    summarize_across_GCMs(var = 'area')
  
  stopifnot(with(df_ba_gcm3, any(!is.na(area) & RCP != 'Current')))
  
  df_ba_gcm3 %>% 
    pivot_wider(values_from = c('area'),
                names_from = 'summary',
                names_prefix = 'area_') %>% 
    mutate(across(c('area_low', 'area_median', 'area_high'),
                  .fns = list('perc' = \(x) x/area_total*100)),
           rcp_year = rcp_label(RCP, years, include_parenth = FALSE)) %>% 
    df_factor()
}
