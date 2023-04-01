#' Rule-of-Thumb for Cohen g
#' 
#' @param g the Cohen g value
#' @param qual optional setting for which rule of thumb to use. Currently only "cohen"
#' @return dataframe with the qualification and source
#' 
#' @examples 
#' g <- -0.2391304
#' th_cohen_g(g)
#' 
#' @details 
#' Cohen's rule of thumb for Cohen g (1988, pp. 147-149):
#' 
#' |\|g\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.05 | negligible |
#' |0.05 < 0.15 | small |
#' |0.15 < 0.25 | medium |
#' |0.25 or more | large |
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'  
#' @export
th_cohen_g <- function(g, qual="cohen"){
  
  if (qual=="cohen") {
    
    ref = "Cohen (1988, pp. 147-149)"
    
    if (abs(g)<0.05) {
      qual = "negligible"}
    else if (abs(g)<0.15) {
      qual = "small"}
    else if (abs(g)<0.25) {
      qual = "medium"}
    else{
      qual = "large"}
  }
  
  results = data.frame(qual, ref)
  
  return(results)
  
}