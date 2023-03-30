#' Mode for Binned Data
#' 
#' @param data the binned data
#' @param allEq indicator on what to do if maximum frequency is equal for more than one category (see details)
#' @param value which value to show in the output (see details) 
#' @returns 
#' A list with:
#' \item{modes}{the mode(s)}
#' \item{modeFD}{frequency density of the mode}
#' 
#' @description 
#' The mode is a measure of central tendency and defined as “the abscissa corresponding to the ordinate of maximum frequency” (Pearson, 1895, p. 345). 
#' A more modern definition would be “the most common value obtained in a set of observations” (Weisstein, 2002). 
#' For binned data the mode is the bin with the highest frequency density. This will have
#' the same result as using the highest frequency if all bins are of equal size.
#' A frequency density is the frequency divided by the bin size (Zedeck, 2014, pp. 144-145).
#' Different methods exist to narrow this down to a single value. See the details for more info on this.
#' 
#' To create the bins, R's *cut()* function could be used.
#' 
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
#' ## Value to return
#' 
#' If *value="midpoint"* is used the modal bin(s) midpoints are shown, using:
#' \deqn{MP_m = \frac{UB_m + LB_m}{2}}
#' Where \eqn{UB_m} is the upper bound of the modal bin, and \eqn{LB_m} the lower bound.
#' 
#' If *value="quadratic"* is used a quadratic curve is made from the midpoint of the
#' bin prior to the modal bin, to the midpoint of the bin after the modal bin.
#' This is done using:
#' \deqn{M = \LB_{m} + \frac{d_1}{d_1 + d_2}\times\left(UB_m - LB_m\right)}
#' With:
#' \deqn{d_1 = FD_m - FD_{m -1}}
#' \deqn{d_2 = FD_m - FD_{m + 1}}
#' Where \eqn{FD_m} is the frequency density of the modal category.
#' 
#' ## Multimode
#' 
#' One small controversy exists if all categories have the same frequency. 
#' In this case none of them has a higher occurence than the others, so none of them would be the mode 
#' (see for example Spiegel & Stephens, 2008, p. 64, Larson & Farber, 2014, p. 69). 
#' This is used when *allEq="none"* and the default.
#' 
#' On a rare occasion someone might argue that if all categories have the same frequency, 
#' then all categories are part of the mode since they all have the highest frequency. 
#' This is used when *allEq="all"*.
#' 
#' The function can return the bins that are the modal bins, by setting *value="none"*.
#' 
#' @examples 
#' data = sample(100, size=100)
#' binData = cut(data, breaks=3)
#' me_mode_bin(binData)
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
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
me_mode_bin <- function(binData, 
                        allEq=c("none", "all"), 
                        value=c("none", "midpoint", "quadratic")){
  
  if (length(allEq)>1) {
    allEq="none"
  }
  
  if (length(value)>1) {
    value="none"
  }
  
  freq = table(binData)
  k = length(freq)
  
  #taken from https://stackoverflow.com/a/61757432
  pattern = "(\\(|\\[)(-*[0-9]+\\.*[0-9]*),(-*[0-9]+\\.*[0-9]*)(\\)|\\])"
  
  #lower, upper bounds and frequency density
  LB = rep(0, k)
  UB = rep(0, k)
  FD = rep(0, k)
  for (i in 1:k) {
    LB[i] = as.numeric(gsub(pattern, "\\2", names(freq)[i]))
    UB[i] = as.numeric(gsub(pattern, "\\3", names(freq)[i]))
    FD[i] = freq[i]/(UB[i] - LB[i])
  }
  
  #Frequency Density of Modal bin
  FDm = max(FD)
  
  #Which index(es)
  iFDm = which(FD==FDm)
  
  #number of modes
  nModes = length(iFDm)
  
  if (nModes == length(freq)) {
    if (allEq=="none") {
      modes = NA
      fMode = NA
    }
  }
  else{
    
    if (value=="none") {
      modes = names(freq)[FD==FDm]
    }
    
    if (value=="midpoint") {
      modes = rep(0, nModes)
      for (i in 1:nModes) {
        modes[i] = (UB[i] + LB[i])/2
      }
      
    }
    
    else if (value=="quadratic") {
      modes = rep(0, nModes)
      
      for (i in 1:nModes) {
        if (i == 1) {d1 = FDm}
        else{d1 = FDm - FD[i-1]}
        
        if (i == k) {d2 = FDm}
        else{d2 = FDm - FD[i+1]}
        
        modes[i] = LB[i] + d1/(d1 + d2)*(UB[i] - LB[i])
      }
    }
  }

  return(list("mode"=modes, "modeFD"=FDm))
  
}