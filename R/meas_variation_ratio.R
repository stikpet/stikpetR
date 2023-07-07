#' Variation Ratio
#' 
#' @description 
#' the Variation Ratio (VR) (Freeman, 1965) is simply the proportion that does not belong to the modal category (Zedeck, 2014, p.406). 
#' 
#' It is a measure of dispersion for categorical data. There are many other measures of dispersion for categorical data. A good start for more info on other measures could be an article from Kader and Perry (2007).
#' 
#' @param data the scores from which to determine the variation ratio
#' 
#' @return VR the variation ratio value
#' 
#' @details 
#' The formula used is (Freeman, 1965, p. 41):
#' \deqn{VR = 1 - p_{m}}
#' With:
#' \deqn{p_{m} = \frac{F_{m}}{n}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} the total sample size
#' \item \eqn{k_m} the number of categories with a frequency equal to \eqn{F_{mode}}
#' \item \eqn{F_{m}} the frequency (count) of the modal category (categories)
#' }
#' 
#' @references 
#' Freeman, L. C. (1965). *Elementary applied statistics: For students in behavioral science*. Wiley.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' #Example 1: dataframe
#' ex1 = df1['mar1']
#' me_variation_ratio(ex1)
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", 
#' "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", 
#' "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' me_variation_ratio(ex2)
#'  
#' @export
me_variation_ratio <- function(data){
  freq = table(data)
  maxFreq = max(freq)
  maxCount = sum(freq==maxFreq)
  
  if (maxCount==length(freq)) {
    warning("no mode in data, so also no variation ratio")
    return()
  }
  VR = 1-maxCount*maxFreq/sum(freq)
  
  return(VR)
}
