library(devtools)
library(roxygen2)
library(rcmdcheck)

library(stikpetR)
#devtools::install_github("stikpet/stikpetR")

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

#in case manual doesn't work and test all examples.
setwd("H:/PeterStatistics/Packages/R")
rcmdcheck("stikpetR")

#Note that \lt and \gt don't work in equations, just use < and >

# DOCUMENTATION
#roxygen2 documentation
#' @param data vector with the scores to determine the mode from
#' @param allEq optional indicator on what to do if maximum frequency is equal for more than one category. Either `"none"` (default), or `"all"`
#' 
#' @seealso 
#' \code{\link{me_mode_bin}}, to determine the mode with binned data
#' 
#' @section Alternatives:
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)

#add GPL-3 licence
#setwd("H:/PeterStatistics/Packages/R/stikpetR")
#use_gpl_license(version = 3, include_future = TRUE)


#setwd("H:/PeterStatistics/Packages/R")
#devtools::create("stikpetR")

# DOCUMENTATION ASCII FIXES
# see documentation errors
setwd("H:/PeterStatistics/Packages/R/stikpetR")
results <- lapply(rd_files, tools::checkRd)
names(results) <- basename(rd_files)
for (file in names(results)) {
  if (length(results[[file]]) > 0) {
    cat("Issues in:", file, "\n")
    print(results[[file]])
  }
}



# preview files
setwd("H:/PeterStatistics/Packages/R/stikpetR")
preview_roxygen_issues <- function(file) {
  lines <- readLines(file, encoding = "UTF-8")
  roxy_lines <- grep("^#'", lines, value = TRUE)
  
  if (any(grepl("[^\x01-\x7F‘’“”éÉàèê–—•]", roxy_lines))) {
    cat("Needs cleaning:", file, "\n")
    bad_lines <- grep("[^\x01-\x7F‘’“”éÉàèê–—•]", roxy_lines, value = TRUE)
    print(bad_lines)
  }
}
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, preview_roxygen_issues))

# replace non-ASCII text in all files
fix_roxygen_quotes <- function(file) {
  lines <- readLines(file, encoding = "UTF-8")
  for (i in seq_along(lines)) {
    if (grepl("^#'", lines[i])) {
      lines[i] <- gsub("[‘’]", "'", lines[i])  # Single quotes
      lines[i] <- gsub("[“”]", '"', lines[i])  # Double quotes
      lines[i] <- gsub("á", "a", lines[i])
      lines[i] <- gsub("ä", "a", lines[i])
      lines[i] <- gsub("à", "a", lines[i])
      lines[i] <- gsub("æ", "a", lines[i])
      lines[i] <- gsub("ā", "a", lines[i])
      
      lines[i] <- gsub("ç", "c", lines[i])
      
      lines[i] <- gsub("é", "e", lines[i])
      lines[i] <- gsub("é", "e", lines[i])
      lines[i] <- gsub("è", "e", lines[i])
      lines[i] <- gsub("É", "E", lines[i])
      
      lines[i] <- gsub("ğ", "g", lines[i])
      
      lines[i] <- gsub("í", "i", lines[i])
      lines[i] <- gsub("ı", "i", lines[i])
      lines[i] <- gsub("İ", "I", lines[i])
      
      lines[i] <- gsub("ñ", "n", lines[i])
      lines[i] <- gsub("ö", "o", lines[i])
      lines[i] <- gsub("Ö", "O", lines[i])
      lines[i] <- gsub("ø", "o", lines[i])
      lines[i] <- gsub("ó", "o", lines[i])
      lines[i] <- gsub("ó", "o", lines[i])
      lines[i] <- gsub("ő", "o", lines[i])
      
      lines[i] <- gsub("Š", "S", lines[i])
      lines[i] <- gsub("š", "s", lines[i])
      lines[i] <- gsub("ş", "s", lines[i])
      
      lines[i] <- gsub("ü", "u", lines[i])
      lines[i] <- gsub("ú", "u", lines[i])
      lines[i] <- gsub("×", "x", lines[i])
      
      
      lines[i] <- gsub("–|—", "-", lines[i])   # dashes
      lines[i] <- gsub("•", "*", lines[i])     # bullet
    }
  }
  writeLines(c(lines, ""), file, useBytes = TRUE)
}

# Apply to all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, fix_roxygen_quotes))



# check for ASCII anywhere including the code
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  non_ascii_lines <- grep("[^\x01-\x7F]", lines)
  if (length(non_ascii_lines) > 0) {
    cat("\nIn file:", file, "\n")
    for (i in non_ascii_lines) {
      cat(sprintf("Line %d: %s\n", i, lines[i]))
    }
  }
}



