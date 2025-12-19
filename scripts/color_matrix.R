# create color matrix for legends

library(tidyverse)
source('src/fig_params.R')
source('src/fig_functions.R')
source('src/SEI_functions.R')

matrix_c12 <- color_matrix_c12()

width = 2
height = 2
filename <- 'figures/sei/c12_color_matrix.png'
ggsave(
  plot = matrix_c12,
  filename = filename, 
  height = height,
  width = width,
  dpi = 600
)
ggsave(
  plot = matrix_c12,
  filename = str_replace(filename, '.png', '.pdf'), 
  height = height,
  width = width,
  device = cairo_pdf,
  family = 'sans'
)
