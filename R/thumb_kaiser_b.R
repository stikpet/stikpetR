#' Rule-of-Thumb for Kaiser b
#' 
#' @description 
#' Simple function to use a rule-of-thumb for the Kaiser b variation measure.
#' 
#' @param b the Cohen b value
#' @param qual optional setting for which rule of thumb to use. Currently only "kaiser"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' Kaiser's rule of thumb for Kaiser b (1968, p. 212):
#' 
#' |\|b\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.70 | terrible |
#' |0.70 < 0.80 | poor  |
#' |0.80 < 0.90 | fair |
#' |0.90 < 0.95 | good |
#' |0.95 < 1.00 | excellent |
#' 
#' @seealso 
#' \code{\link{me_qv}}, to determine Kaiser b
#' 
#' @references 
#' Kaiser, H. F. (1968). A measure of the population quality of legislative apportionment. *American Political Science Review, 62*(1), 208–215. doi:10.2307/1953335
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#'  
#' @export
th_kaiser_b <- function(b, qual="kaiser"){
  
  if (qual=="kaiser") {
    
    ref = "Kaiser (1968, p. 212)"
    
    if (abs(b)<0.70) {
      qual = "terrible"}
    else if (abs(b)<0.80) {
      qual = "poor"}
    else if (abs(b)<0.90) {
      qual = "fair"}
    else if (abs(b)<0.95) {
      qual = "good"}
    else{
      qual = "excellent"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}