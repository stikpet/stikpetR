#' Rule of thumb for Point Biserial Correlation
#' 
#' @description 
#' 
#' Simple function to use a rule-of-thumb for the Point Biserial Correlation.
#' 
#' @param rp the correlation coefficient
#' @param qual optional the rule of thumb to be used. Currently only "cohen"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' 
#' Cohen's rule of thumb for rank biserial (1988, p. 82):
#' 
#' |\|r_p\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.100 | negligible |
#' |0.100 < 0.243 | small |
#' |0.243 < 0.371 | medium |
#' |0.371 < 1 | large |
#' 
#' @references
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' es = 0.6
#' th_point_biserial(es)
#' 
#' @export
th_point_biserial <- function(rp, qual="cohen"){
  # Agnes (2011).
  if (qual=="cohen"){
    ref = "Cohen (1988, p. 82)"
    if (abs(rp) < 0.1){qual = "negligible"}
    else if (abs(rp) < 0.243){qual = "small"}
    else if (abs(rp) < 0.371){qual = "medium"}
    else{qual = "large"}
  }

  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}



