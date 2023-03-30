#' Rule-of-Thumb for Cohen h
#' 
#' @param h the Cohen h value
#' @param qual optional setting for which rule of thumb to use. Currently only 'cohen'
#' @return dataframe with the qualification and source
#' 
#' @examples 
#' h <- -0.2391304
#' th_cohen_h(h)
#' 
#' @details 
#' Cohen's rule of thumb for Cohen g (1988, p. 198):
#' 
#' |\|h\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | negligible |
#' |0.20 < 0.50 | small |
#' |0.50 < 0.80 | medium |
#' |0.80 or more | large |
#' 
#' Note that Cohen actually just lists small = 0.20, medium = 0.50, and large = 0.80.
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
th_cohen_h <- function(h, qual="cohen"){
  
  #Cohen (1988, pp. 184-185)
  if (qual=="cohen") {
    
    ref = "Cohen (1988, p. 198)"
    
    if (abs(h)<0.2) {
      qual = "negligible"}
    else if (abs(h)<0.5) {
      qual = "small"}
    else if (abs(h)<0.8) {
      qual = "medium"}
    else{
      qual = "large"}
  }
  
  results = data.frame(qual, ref)
  
  return(results)
  
}