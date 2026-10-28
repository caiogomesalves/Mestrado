#----Pacotes----

library(tidyverse)
library(ETAS)
library(spatstat)
library(readxl)
library(parallel)

#----Ajuste----

# Base de dados:
jpoly <- list(long = c(134.0, 137.9, 143.1, 144.9, 147.8,
                       137.8, 137.4, 135.1, 130.6),
              lat = c(31.9, 33.0, 33.2,
                      35.2, 41.3, 44.2, 40.2, 38.0, 35.4))

japan.cat <- catalog(japan.quakes, study.start = "1953-05-26",
                     study.end = "1990-01-08", region.poly = jpoly,
                     mag.threshold = 4.5)

# Chute inicial parâmetros:
theta_japan <- c("A" = 0.5, "alpha" = 1.5, "c" = 0.25,
                 "p" = 1.5, "D" = 0.1, "q" = 2, "gamma" = 1)

x_grid_japan <- seq(range(japan.cat$region.poly$long)[1], range(japan.cat$region.poly$long)[2], length.out = 128)
y_grid_japan <- seq(range(japan.cat$region.poly$lat)[1], range(japan.cat$region.poly$lat)[2], length.out = 128)

aj_japan <- ajustar_etas_proprio(japan.cat, theta_japan, 4.5,
                                 x_grid_japan, y_grid_japan, num_cores = 8)
