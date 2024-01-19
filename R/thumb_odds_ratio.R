#' Rules of thumb for Odds Ratio
#' 
#' @description 
#' This function will give a qualification (classification) for a given Odds Ratio
#' 
#' @param or the odds ratio
#' @param qual optional rule of thumb to use. Either "chen" (default), "wuensch", "jones1", "jones2", or "hopkins"
#' 
#' @returns 
#' Dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' If the OR is less than 1, the alternative is used, i.e. 1/OR.
#' 
#' "chen" => Chen et al. (2010, p. 864)
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.68 | negligible |
#' |1.68 < 3.47 | weak |
#' |3.47 < 6.71 | moderate |
#' |6.71 or more | strong |
#' 
#' "hopkins" => Hopkins (1997, tbl. 1)
#' 
#' |\\(OR^{\\ast}\\)| Interpretation|
#' |---|----------|
#' |1.00 < 1.50 | trivial |
#' |1.50 < 3.50 | small |
#' |3.50 < 9.00 | moderate |
#' |9.00 < 32.0 | large |
#' |32.0 < 360 | very large |
#' |360 or more | nearly perfect |
#' 
#' "jones1" => Jones (2014)
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.5 | negligible |
#' |1.5 < 2.5 | small |
#' |2.5 < 4.3 | medium |
#' |4.3 or more | large |
#' 
#' "jones2" => Jones (2014)
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.5 | negligible |
#' |1.5 < 3.5 | small |
#' |2.5 < 9.0 | medium |
#' |9.0 or more | large |
#' 
#' "wuensch" => Wuensch (2009, p. 2)
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.49 | negligible |
#' |1.49 < 3.45 | small |
#' |3.45 < 9 | medium |
#' |9 or more | large |
#' 
#' 
#' @references
#' Chen, H., Cohen, P., & Chen, S. (2010). How big is a big Odds Ratio? Interpreting the magnitudes of Odds Ratios in epidemiological studies. *Communications in Statistics - Simulation and Computation, 39*(4), 860–864. doi:10.1080/03610911003650383
#' 
#' Hopkins, W. G. (2006, August 7). New view of statistics: Effect magnitudes. http://www.sportsci.org/resource/stats/effectmag.html
#' 
#' Jones, K. (2014, June 5). How do you interpret the odds ratio (OR)? ResearchGate. https://www.researchgate.net/post/How_do_you_interpret_the_odds_ratio_OR
#' 
#' Wuensch, K. (2009). Cohen’s conventions for small, medium, and large effects. https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize?action=AttachFile&do=get&target=esize.doc
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' th_odds_ratio(5.23)
#' th_odds_ratio(5.23, qual="wuensch")
#' 
#' @export
th_odds_ratio <- function(or, qual="chen"){
  
  if (or < 1) {
    or = 1/or
  }
  
  #Chen et al. (2010, p. 862).
  if (qual=="chen") {
    src = "Chen et al. (2010, p. 862)"
    if (abs(or) < 1.68) {qual = "negligible"}
    else if (abs(or) < 3.47){qual = "weak"}
    else if (abs(or) < 6.71){qual = "moderate"}
    else {qual = "strong"}
  }
  
  # Hopkins (2006, tbl. 1)
  else if (qual=="hopkins") {
    src = "Hopkins (2006, tbl. 1)"
    if (abs(or) < 1.5) {qual = "trivial"}
    else if (abs(or) < 3.5){qual = "small"}
    else if (abs(or) < 9){qual = "moderate"}
    else if (abs(or) < 32){qual = "large"}
    else if (abs(or) < 360){qual = "very large"}
    else {qual = "nearly perfect"}
  }
  
  # Jones (2014)
  else if (qual=="jones1") {
    src = "Jones (2014)"
    if (abs(or) < 1.5) {qual = "negligible"}
    else if (abs(or) < 2.5){qual = "small"}
    else if (abs(or) < 4.3){qual = "medium"}
    else {qual = "large"}
  }
  
  # Jones (2014)
  else if (qual=="jones2") {
    src = "Jones (2014)"
    if (abs(or) < 1.5) {qual = "negligible"}
    else if (abs(or) < 3.5){qual = "small"}
    else if (abs(or) < 9){qual = "medium"}
    else {qual = "large"}
  }
  
  # Wuensch (2009, p. 2).
  else if (qual=="wuensch") {
    src = "Wuensch (2009, p. 2)"
    if (abs(or) < 1.49) {qual = "negligible"}
    else if (abs(or) < 3.45){qual = "small"}
    else if (abs(or) < 9){qual = "medium"}
    else {qual = "large"}
  }
  
  # the results
  results = data.frame(qual, src)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}