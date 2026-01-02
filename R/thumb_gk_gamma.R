#' Rule-of-Thumb for Goodman-Kruskal gamma
#' 
#' @description
#' Simple function to use a rule-of-thumb for Goodman-Kruskal gamma
#' 
#' @param g the Goodman-Kruskal gamma value
#' @param qual optional setting for which rule of thumb to use. Either "blaikie" (default), "rea-parker", "metsamuuronen" 
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' 
#' @details 
#' The following rules-of-thumb can be used:
#' 
#' *"blaikie"* => Uses Blaikie (2003, p. 100):
#' 
#' |\\|g\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.30 | weak |
#' |0.30 < 0.60 | moderate |
#' |0.60 < 0.75 | strong |
#' |0.75 or more | very strong |
#' 
#' *"rea-parker"* => Rea and Parker (2014, p. 229):
#' 
#' |\\|g\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.30 | low |
#' |0.30 < 0.60 | moderate |
#' |0.60 < 0.75 | strong |
#' |0.75 or more | very strong |
#' 
#' *"metsamuuronen"* => Metsämuuronen (2023, p. 17):
#' 
#' |\\|g\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.14 | negligible |
#' |0.14 < 0.31 | small |
#' |0.31 < 0.45 | medium |
#' |0.45 < 0.62 | large |
#' |0.62 < 0.84 | very large |
#' |0.84 or more | huge |
#' 
#' @references 
#' Blaikie, N. W. H. (2003). *Analyzing quantitative data: From description to explanation*. Sage Publications Ltd.
#' 
#' Metsämuuronen, J. (2023). Somers’ delta as a basis for nonparametric effect sizes: Grissom-Kim PS, Cliff’s d, and Vargha-Delaney A as specific cases of Somers delta. doi:10.13140/RG.2.2.36002.09925
#' 
#' Rea, L. M., & Parker, R. A. (2014). *Designing and conducting survey research: A comprehensive guide* (4th ed.). Jossey-Bass, a Wiley brand.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_gk_gamma <- function(g, qual="blaikie"){
  #Blaikie (2003, p. 100)
  if (qual=="blaikie"){
    ref = "Blaikie (2003, p. 100)"
    if (abs(g) < 0.1){
      qual = "negligible"}
    else if (abs(g) < 0.3){
        qual = "weak"}
    else if (abs(g) < 0.6){
      qual = "moderate"}
    else if (abs(g) < 0.75){
      qual = "strong"}
    else{
      qual = "very strong"}
  }
  
  else if (qual=="rea-parker"){
    ref = "Rea and Parker (2014, p. 229)"
    if (abs(g) < 0.1){
        qual = "negligible"}
    else if (abs(g) < 0.3){
        qual = "low"}
    else if (abs(g) < 0.6){
        qual = "moderate"}
    else if (abs(g) < 0.75){
        qual = "strong"}
    else{
        qual = "very strong"}
  }
  
  #Metsämuuronen (2023, p. 17).
  else if (qual=="metsamuuronen"){
    ref = "Metsämuuronen (2023, p. 17)"
    if (abs(g) < 0.14){
        qual = "negligible"}
    else if (abs(g) < 0.31){
        qual = "small"}
    else if (abs(g) < 0.45){
        qual = "medium"}
    else if (abs(g) < 0.62){
        qual = "large"}
    else if (abs(g) < 0.84){
        qual = "very large"}
    else{
        qual = "huge"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}
