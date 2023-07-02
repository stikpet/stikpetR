#' Rules of Thumb for Cohen d
#' 
#' @param d the Cohen d value
#' @param qual c("cohen", "sawilowsky", "lovakov", "rosenthal") optional the rule of thumb to be used
#' @return dataframe with the qualification and source
#' 
#' @examples  
#' th_cohen_d(0.34)
#' 
#' @details 
#' The following rules-of-thumb can be used:
#' 
#' "cohen" => Cohen (1988, p. 40)
#' 
#' |\|d\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | negligible |
#' |0.20 < 0.50 | small |
#' |0.50 < 0.80 | medium |
#' |0.80 or more | large |
#' 
#' "sawilowsky" => Sawilowsky (2009, p. 599)
#' 
#' |\|d\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.20 | very small |
#' |0.20 < 0.50 | small |
#' |0.50 < 0.80 | medium |
#' |0.80 < 1.20 | large |
#' |1.20 < 2.00 | very large |
#' |2.00 or more | huge |
#' 
#' 
#' "lovakov" => Lovakov and Agadullina (2021, p. 501)
#' 
#' |\|d\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.15 | negligible |
#' |0.15 < 0.36 | small |
#' |0.36 < 0.65 | medium |
#' |0.65 or more | large |
#' 
#' "rosenthal" => Rosenthal (1996, p. 45)
#' 
#' |\|d\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | negligible |
#' |0.20 < 0.50 | small |
#' |0.50 < 0.80 | medium |
#' |0.80 < 1.30 | large |
#' |0.80 or more | very large |
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Lovakov, A., & Agadullina, E. R. (2021). Empirically derived guidelines for effect size interpretation in social psychology. *European Journal of Social Psychology, 51*(3), 485–504. https://doi.org/10.1002/ejsp.2752
#' 
#' Rosenthal, J. A. (1996). Qualitative descriptors of strength of association and effect size. *Journal of Social Service Research, 21*(4), 37–59. https://doi.org/10.1300/J079v21n04_02
#' 
#' Sawilowsky, S. (2009). New effect size rules of thumb. *Journal of Modern Applied Statistical Methods, 8*(2). https://doi.org/10.22237/jmasm/1257035100
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
th_cohen_d <- function(d, qual="sawilowsky"){
  
  #Cohen (1988, p. 40)
  if (qual=="cohen") {
    ref = "Cohen (1988, p. 40)"
    if (abs(d) < 0.2) {
      qual = "negligible"}
    else if (abs(d) < 0.50){
      qual = "small"}
    else if (abs(d) < 0.80){
      qual = "medium"}
    else {
      qual = "large"}
  }
  
  #Lovakov and Agadullina (2021, p. 501)
  else if (qual=="lovakov") {
    ref = "Lovakov and Agadullina (2021, p. 501)"
    if (abs(d) < 0.15) {
      qual = "negligible"}
    else if (abs(d) < 0.35){
      qual = "small"}
    else if (abs(d) < 0.65){
      qual = "medium"}
    else {
      qual = "large"}
  }
  
  #Rosenthal (1996, p. 45)
  else if (qual=="rosenthal") {
    ref = "Rosenthal (1996, p. 45)"
    if (abs(d) < 0.20) {
      qual = "negligible"}
    else if (abs(d) < 0.5){
      qual = "small"}
    else if (abs(d) < 0.8){
      qual = "medium"}
    else if (abs(d) < 1.3){
      qual = "large"}
    else {
      qual = "very large"}
  }
  
  #Sawilowsky (2009, p. 599)
  else if (qual=="sawilowsky") {
    ref = "Sawilowsky (2009, p. 599)"
    if (abs(d) < 0.10) {
      qual = "negligible"}
    else if (abs(d) < 0.20){
      qual = "very small"}
    else if (abs(d) < 0.5){
      qual = "small"}
    else if (abs(d) < 0.8){
      qual = "medium"}
    else if (abs(d) < 1.2){
      qual = "large"}
    else if (abs(d) < 2){
      qual = "very large"}
    else {
      qual = "huge"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}