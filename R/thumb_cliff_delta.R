#' Rule-of-Thumb for Cliff Delta
#' 
#' @description
#' Simple function to use a rule-of-thumb for Cliff Delta
#' 
#' @param d the Cliff Delta value
#' @param qual optional setting for which rule of thumb to use. Either "romano" (default), "metsamuuronen" 
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
#' *"romano"* => Uses Romano et al. (2006, p. 14):
#' 
#' |\\|d\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.15 | negligible |
#' |0.15 < 0.33 | small |
#' |0.33 < 0.47 | medium |
#' |0.47 or more | large |
#' 
#' *"metsamuuronen"* => Metsämuuronen (2023, p. 17):
#' 
#' |\\|d\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.11 | negligible |
#' |0.11 < 0.28 | small |
#' |0.28 < 0.43 | medium |
#' |0.43 or more | large |
#' 
#' @section Before, After and Alternatives:
#' Cliff delta could be converted to Cohen d (Marfo & Okyere, 2019, p. 4) or Vargha-Delaney A (use \code{\link{es_convert}} to convert the effect size measure)
#' 
#' @references 
#' Marfo, P., & Okyere, G. A. (2019). The accuracy of effect-size estimates under normals and contaminated normals in meta-analysis. *Heliyon, 5*(6), e01838. doi:10.1016/j.heliyon.2019.e01838
#' 
#' Metsämuuronen, J. (2023). Somers’ delta as a basis for nonparametric effect sizes: Grissom-Kim PS, Cliff’s d, and Vargha-Delaney A as specific cases of Somers delta. doi:10.13140/RG.2.2.36002.09925
#' 
#' Romano, J., Kromrey, J. D., Coraggio, J., Skowronek, J., & Devine, L. (2006). Exploring methods for evaluating group differences on the NSSE and other surveys: Are the t-test and Cohen’s d indices the most appropriate choices? 1–51.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_cliff_delta <- function(d, qual="romano"){
  #Romano et al. (2006, p. 14).
  if (qual=="romano"){
    ref = "Romano et al. (2006, p. 14)."
    if (abs(d) < 0.15){
      qual = "negligible"}
    else if (abs(d) < 0.33){
      qual = "small"}
    else if (abs(d) < 0.47){
      qual = "medium"}
    else{
      qual = "large"}
    }
  
  #Metsämuuronen (2023, p. 17).
  else if (qual=="metsamuuronen"){
    ref = "Metsämuuronen (2023, p. 17)"
    if (abs(d) < 0.11){
      qual = "negligible"}
    else if (abs(d) < 0.28){
      qual = "small"}
    else if (abs(d) < 0.43){
      qual = "medium"}
    else{
      qual = "large"}
    }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}
