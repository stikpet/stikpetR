#' Rule-of-Thumb for Somers d
#' 
#' @description
#' Simple function to use a rule-of-thumb for Somers d
#' 
#' @param d the Somers d value
#' @param qual optional setting for which rule of thumb to use. Currently only "vargha"
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
#' *"metsamuuronen"* => Uses Metsämuuronen (2023, p. 17):
#' 
#' |\\|d\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.13 | negligible |
#' |0.13 < 0.29 | small |
#' |0.29 < 0.43 | medium |
#' |0.43 < 0.59 | large |
#' |0.59 < 0.81 | very large |
#' |0.81 or more | huge |
#' 
#' @references 
#' Metsämuuronen, J. (2023). Somers’ delta as a basis for nonparametric effect sizes: Grissom-Kim PS, Cliff’s d, and Vargha-Delaney A as specific cases of Somers delta. doi:10.13140/RG.2.2.36002.09925
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_somers_d <- function(d, qual="metsamuuronen"){
  #Metsämuuronen (2023, p. 17)
  if (qual=="metsamuuronen"){
    ref = "Metsämuuronen (2023, p. 17)"
    if (abs(d) < 0.13){
      qual = "negligible"}
    else if (abs(d) < 0.29){
        qual = "small"}
    else if (abs(d) < 0.43){
        qual = "medium"}
    else if (abs(d) < 0.59){
        qual = "large"}
    else if (abs(d) < 0.81){
        qual = "very large"}
    else{
        qual = "huge"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}
