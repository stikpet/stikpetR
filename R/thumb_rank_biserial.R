#' Rule-of-Thumb for Rank Biserial Correlation
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Rank Biserial Correlation.
#' 
#' @param rb the rank-biserial correlation value
#' @param qual optional setting for which rule of thumb to use. Currently only 'cohen'
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' Cohen's rule of thumb for rank-biserial correlation (1988, p. 82):
#' 
#' |\|r_b\|| Interpretation|
#' |---|----------|
#' |0.000 < 0.125 | negligible |
#' |0.125 < 0.304 | small |
#' |0.304 < 0.465 | medium |
#' |0.465 or more | large |
#' 
#' @seealso 
#' \code{\link{es_cohen_w}}, to determine Cohen w
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' es = 0.6
#' th_rank_biserial(es)
#' 
#' @export
th_rank_biserial <- function(rb, qual="cohen"){
  
  #Use Cohen (1988, p. 82)
  ref = "Cohen (1988, p. 82)"
  if (abs(rb) < 0.125){
    qual = "negligible"}
  else if (abs(rb) < 0.304){
    qual = "small"}
  else if (abs(rb) < 0.465){
    qual = "medium"}
  else{
    qual = "large"}
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}