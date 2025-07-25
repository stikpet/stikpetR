#' Rule-of-Thumb for Cohen f
#' 
#' @description 
#' Simple function to use a rule-of-thumb for the Cohen f effect size.
#' 
#' @param f the Cohen f value
#' @param qual optional setting for which rule of thumb to use. Currently only "cohen"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' Cohen's rule of thumb for Cohen f (1988, pp. 285-287):
#' 
#' |\|f\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.25 | small |
#' |0.25 < 0.40 | medium |
#' |0.40 or more | large |
#' 
#' @seealso 
#' \code{\link{es_cohen_f}}, to determine Cohen f
#' 
#' \code{\link{r_rosenthal}}, to determine the Rosenthal correlation, which Cohen called also Cohen f
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_cohen_f <- function(f, qual="cohen"){
  
  if (qual=="cohen") {
    
    ref = "Cohen (1988, pp. 285-287)"
    
    if (abs(f)<0.10) {
      qual = "negligible"}
    else if (abs(f)<0.25) {
      qual = "small"}
    else if (abs(f)<0.40) {
      qual = "medium"}
    else{
      qual = "large"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}



