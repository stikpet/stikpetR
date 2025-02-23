#' Rule-of-Thumb for Cohen w
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Cohen w effect size.
#' 
#' The measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CohenW.html)
#' 
#' @param w the Cohen w value
#' @param qual optional setting for which rule of thumb to use. Currently only 'cohen'
#' 
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' 
#' @details 
#' Cohen's rule of thumb for Cohen w (1988, p. 227):
#' 
#' |\|w\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.30 | small |
#' |0.30 < 0.50 | medium |
#' |0.50 or more | large |
#' 
#' 
#' @section Before, After and Alternatives:
#' Before using this function you need to obtain a Cohen w value:
#' \code{\link{es_cohen_w}}, to determine Cohen w
#' 
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' es = 0.6
#' th_cohen_w(es)
#' 
#' @export
th_cohen_w <- function(w, qual="cohen"){
  
  #Use Cohen (1988, p. 227)
  ref = "Cohen (1988, p. 227)"
  if (abs(w) < 0.1){
    qual = "negligible"}
  else if (abs(w) < 0.3){
    qual = "small"}
  else if (abs(w) < 0.5){
    qual = "medium"}
  else{
    qual = "large"}
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}


