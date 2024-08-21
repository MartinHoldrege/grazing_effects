# code from Daniel S
#--- Paths ------
# database
source('src/paths.R')
dir_results <- "data_processed/site_means"

#--- Connect to database ------
rSOILWAT2::dbW_setConnection(
  dbFilePath = file.path(
    path_sw,
    "dbWeatherData_STEPWAT2_200sites.sqlite3"
  )
)


#--- Obtain table with sites -----
xsites <- rSOILWAT2::dbW_getSiteTable() |> as.data.frame()
Nsites <- nrow(xsites)

monthly <- tidyr::expand_grid(site = xsites$Site_id,
                              month = 1:12) 
monthly$PPT <- NA
monthly$Tmean <- NA
#--- Loop over sites ------
xsites[, c("MAT_C", "MAP_mm", "CorrTP", "CorrTP2")] <- NA

pb <- utils::txtProgressBar(max = Nsites, style = 3L)

for (k in seq_len(Nsites)) {
  #--- * Extract weather from database ------
  wdl <- rSOILWAT2::dbW_getWeatherData(
    Label = xsites[k, "Label"],
    Scenario = 'Current'
  )

  #--- * Calculate annual temperature and precipitation amount ------
  clim <- rSOILWAT2::calc_SiteClimate(wdl)
  xsites[k, "MAT_C"] <- clim[["MAT_C"]]
  xsites[k, "MAP_mm"] <- 10 * clim[["MAP_cm"]]
  monthly$PPT[monthly$site == k] <- 10* clim[["meanMonthlyPPTcm"]]
  monthly$Tmean[monthly$site == k] <- clim[["meanMonthlyTempC"]]

  #--- * Calculate monthly correlation between precipitation and temperature ---
  wdt <- rSOILWAT2::dbW_weatherData_to_monthly(
    wdl,
    funs = rSOILWAT2::weather_dataAggFun()[c("Tmax_C", "Tmin_C", "PPT_cm")]
  ) |>
    as.data.frame()

  wdt[["Tmean_C"]] <- rowMeans(wdt[, c("Tmax_C", "Tmin_C")])

  xsites[k, "CorrTP2"] <- by(
    data = wdt[c("Tmean_C", "PPT_cm")],
    INDICES = wdt[["Year"]],
    FUN = function(x) stats::cor(x[, 1], x[, 2])
  ) |>
    as.vector() |>
    mean()

wdt
  # correlation of monthly means
  # 'type 1' corrTP (not what we want to use)
  # xsites[k, 'CorrTP'] <- wdt |>
  #   dplyr::group_by(Month) |>
  #   dplyr::summarize(Tmean_C = mean(Tmean_C),
  #             PPT_cm = mean(PPT_cm)) |>
  #   dplyr::summarize(CorrTP = stats::cor(Tmean_C, PPT_cm)) |> 
  #   dplyr::pull(CorrTP)
  # utils::setTxtProgressBar(pb, value = k)
}

close(pb)



#--- Write result to disk ---
utils::write.csv(
  xsites,
  file = file.path(dir_results, "dbWeather_200sites.csv"),
  row.names = FALSE
)

utils::write.csv(
  monthly,
  file = file.path(dir_results, "dbWeather_200sites_monthly.csv"),
  row.names = FALSE
)

# hist(xsites[["CorrTP"]])


#--- Clean up ------
rSOILWAT2::dbW_disconnectConnection()
