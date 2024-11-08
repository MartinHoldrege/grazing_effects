# Contains functions for predicting fire probability


# fire probability, predicted by equation in Holdrege et al. (2024)
predict_fire <- function(mat, map, psp, afg, pfg) {
  if(any(mat > 70 | mat < - 70)) {
    stop('mat should be in C')
  }
  if (any(psp<0 | psp>1)) {
    stop('psp has a range of 0-1')
  }
  mat <-  mat + 273.15
  y = -92.86 + 14.01*log10(afg + 1) - 0.8117*(log10(afg + 1))^2 +
    0.5054*pfg^0.5 - 0.03730*pfg +0.02672*mat + 49.25*log10(map + 1) - 
    8.236*log10(map + 1)^2 -7.505*log10(psp + 0.001) - 3.118*log10(psp + 0.001)^2 -
    4.047*log10(afg + 1)*log10(map + 1)
  
  p <- 1/(1 + exp(-y))
  p
}

# testing
if (FALSE) {
  afg <- seq(0, 120, by = 1)
  p <- predict_fire(mat = 283, map = 250, psp = 0.2, afg = afg,
                    pfg = 30)
  plot(afg, p*100) 
}


