#' Rule-of-Thumb for Cohen h
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Cohen h effect size.
#' 
#' @param h the Cohen h value
#' @param qual optional setting for which rule of thumb to use. Currently only 'cohen'
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
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
#' @seealso 
#' \code{\link{es_cohen_h}}, to determine Cohen h
#' 
#' \code{\link{es_cohen_h_os}}, to determine Cohen h', then use \code{\link{es_convert}}
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' es = 0.6
#' th_cohen_d(es)
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
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}