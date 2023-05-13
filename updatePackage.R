library(devtools)
library(roxygen2)
library(rcmdcheck)

library(stikpetR)

#Note to self 1
#If a global library needs to be installed:
#add it WITHOUT comma to the help_import_global file
#then add it to the depends WITH comma in the dependencies
#update the documentation

#update documentation
setwd("H:/PeterStatistics/Packages/R")
devtools::document("stikpetR")

#build manual
#thanks a lot to https://tex.stackexchange.com/questions/125274/error-font-ts1-zi4r-at-540-not-found
setwd("H:/PeterStatistics/Packages/R/stikpetR")
build_manual(pkg = ".")


#update GitHub
#see https://rfortherestofus.com/2021/02/how-to-use-git-github-with-r/


#Create a vignette
#usethis::use_vignette("single_binary")

#first clean the environment, restart R and then use:
install_github("stikpet/stikpetR", build_vignettes = TRUE)

help(package = "stikpetR")

#Troubleshooting manual
setwd("H:/PeterStatistics/Packages/R/stikpetR")
check_man()

#in case manual doesn't work.
setwd("H:/PeterStatistics/Packages/R")
rcmdcheck("stikpetR")


#add GPL-3 licence
#setwd("H:/PeterStatistics/Packages/R/stikpetR")
#use_gpl_license(version = 3, include_future = TRUE)


#setwd("H:/PeterStatistics/Packages/R")
#devtools::create("stikpetR")