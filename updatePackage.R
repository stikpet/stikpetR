library(devtools)
library(roxygen2)
library(rcmdcheck)


setwd("H:/PeterStatistics/Packages/R")

#devtools::create("stikpetR")

#add GPL-3 licence
#setwd("H:/PeterStatistics/Packages/R/stikpetR")
#use_gpl_license(version = 3, include_future = TRUE)



library(stikpetR)

#update documentation
setwd("H:/PeterStatistics/Packages/R")
devtools::document("stikpetR")

#build manual
#thanks a lot to https://tex.stackexchange.com/questions/125274/error-font-ts1-zi4r-at-540-not-found
setwd("H:/PeterStatistics/Packages/R/stikpetR")
build_manual(pkg = ".")

setwd("H:/PeterStatistics/Packages/R/stikpetR")
check_man()


setwd("H:/PeterStatistics/Packages/R")
rcmdcheck("stikpetR")

