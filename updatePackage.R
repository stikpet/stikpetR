library(devtools)
library(roxygen2)

setwd("H:/PeterStatistics/Packages/R")
devtools::create("stikpetR")

#update documentation
devtools::document("stikpetR")

library(stikpetR)


install_github("stikpet/Rpackage")