#' Rule-of-Thumb for Cohen g
#' 
#' @description 
#' Simple function to use a rule-of-thumb for the Cohen g effect size.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/3DEmngmws2U) and the effect size is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CohenG.html)
#' 
#' 
#' @param g the Cohen g value
#' @param qual optional setting for which rule of thumb to use. Currently only "cohen"
#' 
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
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
#' 
#' @section Before, After and Alternatives:
#' \code{\link{es_cohen_g}}, to determine Cohen g
#' 
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples
#' es = 0.6
#' th_cohen_g(es)
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
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}