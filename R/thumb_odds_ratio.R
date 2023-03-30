#' Rules of thumb for Odds Ratio
#' 
#' @param or the odds ratio
#' @param qual c("chen", "wuensch", "jones") optional rule of thumb to use (default is "chen")
#' @return data frame with the qualification/classification  and source
#' 
#' @details 
#' If the OR is less than 1, the alternative is used, i.e. 1/OR.
#' 
#' 
#' Chen et al. (2010, p. 862):
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.68 | negligible |
#' |1.68 < 3.47 | weak |
#' |3.47 < 6.71 | moderate |
#' |6.71 or more | strong |
#' 
#' 
#' Wuensch (2009, p. 2):
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.49 | negligible |
#' |1.49 < 3.45 | small |
#' |3.45 < 9 | medium |
#' |9 or more | large |
#' 
#' Jones (2014):
#' 
#' |OR| Interpretation|
#' |---|----------|
#' |1.00 < 1.5 | negligible |
#' |1.5 < 2.5 | small |
#' |2.5 < 4.3 | medium |
#' |4.3 or more | large |
#' 
#' 
#' @examples 
#' th_odds_ratio(5.23)
#' th_odds_ratio(5.23, qual="wuensch")
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references
#' Chen, H., Cohen, P., & Chen, S. (2010). How big is a big Odds Ratio? Interpreting the magnitudes of Odds Ratios in epidemiological studies. *Communications in Statistics - Simulation and Computation, 39*(4), 860–864. https://doi.org/10.1080/03610911003650383
#' 
#' Jones, K. (2014, June 5). How do you interpret the odds ratio (OR)? ResearchGate. https://www.researchgate.net/post/How_do_you_interpret_the_odds_ratio_OR
#' 
#' Wuensch, K. (2009). Cohen’s conventions for small, medium, and large effects. https://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize?action=AttachFile&do=get&target=esize.doc  
#' 
#' @export
th_odds_ratio <- function(or, qual="chen"){
  
  if (or < 1) {
    or = 1/or
  }
  
  #Chen et al. (2010, p. 862).
  if (qual=="chen") {
    if (abs(or) < 1.68) {
      qual = "negligible"}
    else if (abs(or) < 3.47){
      qual = "weak"}
    else if (abs(or) < 6.71){
      qual = "moderate"}
    else {
      qual = "strong"}
  }
  
  #Wuensch (2009, p. 2).
  if (qual=="wuensch") {
    if (abs(or) < 1.49) {
      qual = "negligible"}
    else if (abs(or) < 3.45){
      qual = "small"}
    else if (abs(or) < 9){
      qual = "medium"}
    else {
      qual = "large"}
  }
  
  #Jones (2014)
  if (qual=="jones") {
    if (abs(or) < 1.5) {
      qual = "negligible"}
    else if (abs(or) < 2.5){
      qual = "small"}
    else if (abs(or) < 4.3){
      qual = "medium"}
    else {
      qual = "large"}
  }
  
  return(qual)
}

th_odds_ratio(5.23)
th_odds_ratio(5.23, qual="chen")