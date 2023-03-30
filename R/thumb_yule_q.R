#' Rules of thumb for Yule Q
#' 
#' @param q the Yule Q value
#' @param qual optional for which rule-of-thumb to use. Currently only "glen"
#' @return the qualification/classification
#' 
#' @details 
#' 
#' Glen rule of thumb for Yule Q (2017):
#' 
#' |\|Q\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.30 | negligible |
#' |0.30 < 0.50 | moderate |
#' |0.50 < 0.70 | substantial |
#' |0.70 or more | very strong |
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references
#' Glen, S. (2017, August 16). Gamma Coefficient (Goodman and Kruskal’s Gamma) & Yule’s Q. Statistics How To. https://www.statisticshowto.com/gamma-coefficient-goodman-kruskal/
#' 
#' @examples 
#' q = 0.4285714
#' th_yule_q(q)
#' 
#' @export
th_yule_q <- function(q, qual="glen"){
  
  if (qual=="glen") {
    
    if (abs(q)<0.29) {
      qual = "very small"}
    else if (abs(q)<0.49) {
      qual = "moderate"}
    else if (abs(q)<0.69) {
      qual = "substantion"}
    else{
      qual = "very strong"}
  }
  
  return(qual)
  
}