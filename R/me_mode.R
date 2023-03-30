#' Mode
#' 
#' @param data the scores to determine the mode from
#' @param allEq indicator on what to do if maximum frequency is equal for more than one category (see details)
#' @returns 
#' A list with:
#' \item{modes}{the mode(s)}
#' \item{modeFreq}{frequency of the mode}
#' 
#' @description 
#' The mode is a measure of central tendency and defined as “the abscissa corresponding to the ordinate of maximum frequency” (Pearson, 1895, p. 345). 
#' A more modern definition would be “the most common value obtained in a set of observations” (Weisstein, 2002). 
#' The word mode might even come from the French word 'mode' which means fashion. Fashion is what most people wear, 
#' so the mode is the option most people chose.
#' 
#' If one category has the highest frequency this category will be the modal category and if two or more categories have 
#' the same highest frequency each of them will be the mode. If there is only one mode the set is sometimes called unimodal, 
#' if there are two it is called bimodal, with three trimodal, etc. For two or more, thse term multimodal can also be used.
#' 
#' An advantage of the mode over many other measures of central tendency (like the median and mean), is that it can be determined for already 
#' nominal data types. 
#' 
#' @details 
#' One small controversy exists if all categories have the same frequency. 
#' In this case none of them has a higher occurence than the others, so none of them would be the mode 
#' (see for example Spiegel & Stephens, 2008, p. 64, Larson & Farber, 2014, p. 69). 
#' This is used when *allEq="none"* and the default.
#' 
#' On a rare occasion someone might argue that if all categories have the same frequency, 
#' then all categories are part of the mode since they all have the highest frequency. 
#' This is used when *allEq="all"*.
#' 
#' @seealso 
#' \code{\link{me_mode_bin}}, to determine the mode with binned data
#' 
#' @examples 
#' data = c("a", "a", "a", "b", "b", "b", "c")
#' me_mode(data)
#' me_mode(data, allEq="all")
#' 
#' data = c(1, 1, 1, 2, 2, 2, 3)
#' me_mode(data)
#' 
#' data = c(1, 2, 1, 2, 3, 3)
#' me_mode(data)
#' 
#' @references 
#' Larson, R., & Farber, E. (2014). *Elementary statistics: Picturing the world* (6th ed.). Pearson.
#' 
#' Pearson, K. (1895). Contributions to the mathematical theory of evolution. II. Skew variation in homogeneous material. *Philosophical Transactions of the Royal Society of London. (A.), 186*, 343–414. https://doi.org/10.1098/rsta.1895.0010
#' 
#' Spiegel, M. R., & Stephens, L. J. (2008). *Schaum’s outline of theory and problems of statistics* (4th ed.). McGraw-Hill.
#' 
#' Weisstein, E. W. (2002). *CRC concise encyclopedia of mathematics* (2nd ed.). Chapman & Hall/CRC.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
me_mode <- function(data, allEq = c("none", "all")){
  
  if (length(allEq)>1) {
    allEq="none"
  }
  
  freq = table(data)
  fMode = max(freq)
  modes = names(freq)[freq == fMode]
  
  if (is.numeric(data)) {
    modes = as.numeric(modes)  
  }
  
  if (length(modes) == length(freq)) {
    if (allEq=="none") {
      modes = NA
      fMode = NA
    }
  }
  
  return(list("mode"=modes, "modeFreq"=fMode))
  
}