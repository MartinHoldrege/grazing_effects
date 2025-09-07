# functions for the 06_sei_attribution.R script
# (and others as needed)


prepare_info_dom <- function(info) {
  info %>% 
    mutate(id_noPFT = str_replace(id, '_[[:alpha:]]+_biomass', '_dom-driver'),
           id_noGCMPFT = str_replace(id_noGCM, '_[[:alpha:]]+_biomass', '_dom-driver')) %>% 
    group_by(run, graze, RCP, years, GCM) %>% 
    group_split() %>% 
    map(\(x) arrange(x, PFT))
}
