#' Rules of thumb for Yule Q
#' @description 
#' Simple function to use a rule-of-thumb for Yule Q effect size.
#' 
#' @param q the Yule Q value
#' @param qual optional for which rule-of-thumb to use. Currently only "glen"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
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
#' @seealso 
#' \code{\link{es_bin_bin}}, to determine Yule Q
#' 
#' @references
#' Glen, S. (2017, August 16). Gamma Coefficient (Goodman and Kruskal’s Gamma) & Yule’s Q. Statistics How To. https://www.statisticshowto.com/gamma-coefficient-goodman-kruskal/
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' q = 0.6
#' th_yule_q(q)
#' 
#' @export
th_yule_q <- function(q, qual="glen"){
  
  if (qual=="glen") {
    ref = "Glen (2017)"
    if (abs(q)<0.29) {qual = "very small"}
    else if (abs(q)<0.49) {qual = "moderate"}
    else if (abs(q)<0.69) {qual = "substantion"}
    else{qual = "very strong"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}