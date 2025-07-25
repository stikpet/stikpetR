#' Rule-of-Thumb for Cramer V
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Cramer V effect size. Note however that many will actually use the rule-of-thumb for Cohen w and convert Cramer V to Cohen w first.
#' 
#' The measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CramerV.html)
#' 
#' @param v the Cramer V value
#' @param qual optional setting for which rule of thumb to use. Either "rea-parker" (default), "akoglu", "calamba-rustico"
#' 
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' 
#' @details 
#' *"rea-parker"* => Uses Rea and Parker (1992, p. 203):
#' 
#' |\\|v\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.20 | weak |
#' |0.20 < 0.40 | moderate |
#' |0.40 < 0.60 | relatively strong |
#' |0.60 < 0.80 | strong |
#' |0.80 or more | very strong |
#' 
#' *"akoglu"* => Uses Akoglu (2018, p. 92):
#' 
#' |\\|v\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.05 | very weak |
#' |0.05 < 0.10 | weak |
#' |0.10 < 0.15 | moderate |
#' |0.15 < 0.25 | strong |
#' |0.25 or more | very strong |
#' 
#' *"calamba-rustico"* => Uses Calamba and Rustico (2019, p. 7):
#' 
#' |\\|v\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.15 | very weak |
#' |0.15 < 0.20 | weak |
#' |0.20 < 0.25 | moderate |
#' |0.25 < 0.30 | moderately strong |
#' |0.30 < 0.35 | strong |
#' |0.35 < 0.50 | worrisomely strong |    
#' |0.50 or more | redundant |
#' 
#' Note that the original source has a gap from 0.40 < 0.50, I added this to the 'worrisomely strong' category.
#' 
#' 
#' @section Before, After and Alternatives:
#' Before using this function you need to obtain a Cramer v value:
#' \code{\link{es_cramer_v_gof}}, to determine Cramer V for a Goodness-of-Fit test.
#' \code{\link{es_cramer_v_ind}}, to determine Cramer V for a test of independence.
#' 
#' 
#' @references 
#' Akoglu, H. (2018). User's guide to correlation coefficients. *Turkish Journal of Emergency Medicine, 18*(3), 91-93. doi:10.1016/j.tjem.2018.08.001
#' 
#' Calamba, S. S., & Rustico, E. M. P. (2019). Usefulness of code of ethics for professional accountants in resolving ethical conflicts in the Philippines.
#' 
#' Rea, L. M., & Parker, R. A. (1992). *Designing and conducting survey research: A comprehensive guide*. Jossey-Bass Publishers.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples
#' es = 0.6
#' th_cramer_v(es)
#' 
#' 
#' @export
th_cramer_v <- function(v, qual="rea-parker"){
  
  #Rea and Parker (1992, p. 203).
  if (qual=="rea-parker") {
    ref = "Rea and Parker (1992, p. 203)"
    if (abs(v) < 0.10) {
      qual = "negligible"}
    else if (abs(v) < 0.20){
      qual = "weak"}
    else if (abs(v) < 0.40){
      qual = "moderate"}
    else if (abs(v) < 0.60){
      qual = "relatively strong"}
    else if (abs(v) < 0.80){
      qual = "strong"}
    else {
      qual = "very strong"}
  }
  
  #Akoglu (2018, p. 92)
  else if (qual=="akoglu") {
    ref = "Akoglu (2018, p. 92)"
    if (abs(v) < 0.05) {
      qual = "very weak"}
    else if (abs(v) < 0.10){
      qual = "weak"}
    else if (abs(v) < 0.15){
      qual = "moderate"}
    else if (abs(v) < 0.25){
      qual = "strong"}
    else {
      qual = "very strong"}
  }
  
  # Calamba and Rustico (2019, p. 7)
  else if (qual=="calamba-rustico") {
    ref = "Calamba and Rustico (2019, p. 7)"
    if (abs(v) < 0.15) {
      qual = "very weak"}
    else if (abs(v) < 0.20){
      qual = "weak"}
    else if (abs(v) < 0.25){
      qual = "moderate"}
    else if (abs(v) < 0.30){
      qual = "moderately strong"}
    else if (abs(v) < 0.35){
      qual = "strong"}
    else if (abs(v) < 0.50){
      qual = "worrisomely strong"}
    else {
      qual = "redundant"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}



