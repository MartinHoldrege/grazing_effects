# functions for the 06_sei_attribution.R script
# (and others as needed)


prepare_info_dom <- function(info) {
  out <- info %>% 
    mutate(id_noPFT = str_replace(id, '_[[:alpha:]]+_biomass', '_dom-driver'),
           id_noGCMPFT = str_replace(id_noGCM, '_[[:alpha:]]+_biomass', '_dom-driver')) %>% 
    group_by(run, graze, RCP, years, GCM) %>% 
    group_split() %>% 
    map(\(x) arrange(x, PFT))
  stopifnot(map_dbl(out, nrow) == 3)
  out
}

#' create raster where grid cells provide the dominant driver of change
#'
#' @param info_l list of dataframs from prepare_info_dom
#' @param q_prop_rel raster from relativize_q_prop()
dom_driver <- function(info_l, q_prop_rel) {
  
  out <- map(info_l, function(df) {
    stopifnot(nrow(df) == 3)
    r <- q_prop_rel[[df$id]]
    # warnings here are normal
    out <- terra::which.max(r)
    names(out) <- unique(df$id_noPFT)
    out
  })
  terra::rast(out)
}

prepare_info_dom_mode <- function(info_l) {
  info_dom_l <- bind_rows(info_l) %>% 
    select(-PFT, -id, -id_noGCM, -matches('id_ref')) %>% 
    distinct() %>% 
    group_by(id_noGCMPFT) %>% 
    group_split()
}

# calculate dominant drivers across GCMs and then take mode
# across GCMs
dom_driver_mode <- function(info_l, q_prop_rel) {
  q_dom <- dom_driver(info_l, q_prop_rel)
  
  info_dom_l <- prepare_info_dom_mode(info_l)
  
  out <- map(info_dom_l, \(df) {
    stopifnot(nrow(df) == 13 | (nrow(df) == 1 & df$GCM == 'Current'))
    r <- terra::modal(q_dom[[df$id_noPFT]])
    names(r) <- unique(df$id_noGCMPFT)
    r
  }) 
  terra::rast(out)
}

#' summarize a spatraster across GCMs
#'
#' @param info dataframe containing id and id_noGCM columns
#' @param r raster
#' @param fun summary function name to be passed to app() 
#'
#' @returns
#' a spatraster w/ fewer bands
raster_gcm_smry <- function(info, r, fun = 'median') {
  info_l <- info%>% 
    group_by(id_noGCM) %>% 
    group_split()
  r_l <- map(info_l, \(df) {
    stopifnot(nrow(df) == 13 | (nrow(df) == 1 & df$GCM == 'Current'),
              lu(df$id_noGCM) == 1)
    r <- r[[df$id]]
    out <- terra::app(r, fun)
    names(out) <- unique(df$id_noGCM)
    out
  })
  terra::rast(r_l)
}

#' @description
#' the dominant driver layers are masked, using that mask and applying
#' it to the q proportional change layers. Importantly both inputs need to
#' be of the same type of change (i.e. comparing climate scenarios, or
#' comparing grazing scenarios within a climate scenario)
#' 
#'
#' @param r raster of proportional q changes
#' @param dom raster of dominant drivers (where value of 4 means no dominant driver,
#' b/ below threshold change )
mask_using_dom_layer <- function(r, dom) {
  
  into <- c("variable", "type", "RCP", "years",  "graze")
  if(any(str_detect(names(dom), 'median'))) {
    info1 <- create_rast_info(dom, into = c(into, 'summary')) %>% 
      filter(summary == 'median')
  } else {
    info1 <- create_rast_info(dom, into = into) 
  }

  info2 <- create_rast_info(r, into = c("variable", "type", "RCP", "years", 
                                        "graze"))
  
  info_comb <- info2 %>% 
    left_join(info1, by = c('run', 'run2', 'RCP', "years", 'graze'),
              suffix = c('', '_mask'))
  
  stopifnot(all(!is.na(info_comb$id_mask)),# check succesful join
            nrow(info_comb) == nrow(info2)
            ) 
  
  mask <- terra::ifel(dom == 4, NA, 1)
  
  terra::mask(r[[info_comb$id]], mask[[info_comb$id_mask]])
}
