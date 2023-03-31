library(devtools)
library(roxygen2)

setwd("H:/PeterStatistics/Packages/R")
#devtools::create("stikpetR")

#update documentation
devtools::document("stikpetR")

library(stikpetR)

#usethis::use_readme_rmd()

#install_github("stikpet/stikpetR")

build_manual(pkg = "stikpetR")

#tinytex::tlmgr_install("makeindex")
#tinytex::install_tinytex(TRUE)
#install.packages("texlive-fonts-extra")

library(rcmdcheck)
rcmdcheck("stikpetR")



library(TinyTex)
