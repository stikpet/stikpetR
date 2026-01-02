#' Rule-of-Thumb for Vargha-Delaney A
#' 
#' @description
#' Simple function to use a rule-of-thumb for Vargha-Delaney A
#' 
#' @param a the Vargha-Delaney A value
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
#' *"vargha"* => Uses Vargha and Delaney (2000, p. 106):
#' 
#' |\\|0.5 - A\\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.06 | negligible |
#' |0.06 < 0.14 | small |
#' |0.14 < 0.21 | medium |
#' |0.21 or more | large |
#' 
#' @references 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101–132. doi:10.3102/10769986025002101
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_vda <- function(a, qual="vargha"){
  #Vargha and Delaney (2000, p. 106)
  if (qual=="vargha"){
    ref = "Vargha and Delaney (2000, p. 106)"
    if (abs(0.5 - a) < 0.06){
      qual = "negligible"}
    else if (abs(0.5 - a) < 0.14){
        qual = "small"}
    else if (abs(0.5 - a) < 0.21){
        qual = "medium"}
    else{
        qual = "large"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}
