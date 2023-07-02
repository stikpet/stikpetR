#' Mode for Binned Data
#' 
#' @description 
#' The mode is a measure of central tendency and defined as “the abscissa corresponding to the ordinate of maximum frequency” (Pearson, 1895, p. 345). A more modern definition would be “the most common value obtained in a set of observations” (Weisstein, 2002).
#' 
#' For binned data the mode is the bin with the highest frequency density. This will have the same result as using the highest frequency if all bins are of equal size. A frequency density is the frequency divided by the bin size (Zedeck, 2014, pp. 144-145). Different methods exist to narrow this down to a single value. See the notes for more info on this.
#' 
#' The word mode might even come from the French word 'mode' which means fashion. Fashion is what most people wear, so the mode is the option most people chose.
#' 
#' If one category has the highest frequency this category will be the modal category and if two or more categories have the same highest frequency each of them will be the mode. If there is only one mode the set is sometimes called unimodal, if there are two it is called bimodal, with three trimodal, etc. For two or more, thse term multimodal can also be used.
#' 
#' An advantage of the mode over many other measures of central tendency (like the median and mean), is that it can be determined for already nominal data types.
#' 
#' @param data list or dataframe
#' @param bins dataframe with the lower and upper bounds
#' @param allEq optional indicator on what to do if maximum frequency is equal for more than one category. Either "none" (default) or "all"
#' @param value optional which value to show in the output. Either "none" (default), "midpoint", or "quadratic"
#' 
#' @returns 
#' A dataframe with
#' \item{mode}{the mode(s)}
#' \item{mode fd}{frequency density of the mode}
#' 
#' @details 
#' **Value to return**
#' 
#' If *value="midpoint"* is used the modal bin(s) midpoints are shown, using:
#' \deqn{MP_m = \frac{UB_m + LB_m}{2}}
#' Where \eqn{UB_m} is the upper bound of the modal bin, and \eqn{LB_m} the lower bound.
#' 
#' If *value="quadratic"* is used a quadratic curve is made from the midpoint of the bin prior to the modal bin, to the midpoint of the bin after the modal bin. This is done using:
#' \deqn{M = LB_{m} + \frac{d_1}{d_1 + d_2}\times\left(UB_m - LB_m\right)}
#' 
#' With:
#' \deqn{d_1 = FD_m - FD_{m -1}}
#' \deqn{d_2 = FD_m - FD_{m + 1}}
#' 
#' Where \eqn{FD_m} is the frequency density of the modal category.
#' 
#' **Multimode**
#' 
#' One small controversy exists if all categories have the same frequency. In this case none of them has a higher occurence than the others, so none of them would be the mode (see for example Spiegel & Stephens, 2008, p. 64, Larson & Farber, 2014, p. 69). This is used when *allEq="none"* and the default.
#' 
#' On a rare occasion someone might argue that if all categories have the same frequency, then all categories are part of the mode since they all have the highest frequency. This is used when *allEq="all"*.
#' 
#' The function can return the bins that are the modal bins, by setting *value="none"*.
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
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
me_mode_bin <- function(data, bins, allEq="none", value="none"){
  
  k = nrow(bins)
  
  tab = data.frame(matrix(ncol=2, nrow=0))
  colnames(tab) = c("f", "fd")
  for (i in 1:k ){
    f = sum(ex1 >= bins[i,1]) - sum(ex1 >= bins[i,2])
    fd = f/(bins[i,2] - bins[i,1])
    tab[nrow(tab) + 1, ] = c(f, fd)
  }
  
  modeFD = max(tab[,'fd'])
  nModes = sum(tab[,'fd'] == modeFD)
  
  if ((nModes==k) & (allEq=="none")){
    mode = "none"
    modeFD = "none"
  }
  else {
    if (value=="midpoint"){
      ff = 0
      for (i in 1:k){
        if (tab[i, "fd"]==modeFD){
          newMode = (bins[i,1] + bins[i,2])/2
          if (ff==0){
            mode=newMode
            ff = ff + 1}
          else {
            mode = paste(mode, ",", newMode)
          }
        }
      }
    }
    else if (value=="quadratic"){
      ff = 0
      for (i in 1:k){
        if (tab[i, "fd"]==modeFD){
          if (i==1){
            d1 = modeFD}
          else if (i==k){
            d2 = modeFD}
          else{
            d1 = modeFD - tab[(i - 1), "fd"]
            d2 = modeFD - tab[(i + 1), "fd"]}
          
          newMode = bins[i, 1] + d1/(d1 + d2) * (bins[i, 2] - bins[i,1])
          
          if (ff==0){
            mode=newMode
            ff = ff + 1}
          else {
            mode = paste(mode, ",", newMode)
          }
        }
      }
    }
    else if (value=="none"){
      ff = 0
      for (i in 1:k){
        if (tab[i, "fd"]==modeFD){
          newMode = paste(bins[i,1], "<", bins[i,2])
          if (ff==0){
            mode=newMode
            ff = ff + 1}
          else {
            mode = paste(mode, ",", newMode)
          }
        }
      }
    }
  }
  
  res <- data.frame(mode, modeFD)
  colnames(res)<-c("mode", "mode fd.")
  
  return (res)
}