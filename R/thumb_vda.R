#' Rule-of-Thumb for Vargha-Delaney A
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Vargha-Delaney A
#' 
#' @param a the Vargha-Delaney A value
#' @param qual optional setting for which rule of thumb to use. Either "vd" (default), "sawilowsky", "cohen", "lovakov", "rosenthal"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' Vargha and Delaney (2000, p. 106):
#' 
#' |\|0.5 - A\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.06 | negligible |
#' |0.11 < 0.14 | small |
#' |0.28 < 0.21 | medium |
#' |0.21 or more | large |
#' 
#' "sawilowsky", "cohen", "lovakov", and "rosenthal" will use the rule-of-thumb from Cohen d, by converting VD-A first to a (Glass) Rank Biserial (Cliff delta) then to Cohen d, and use the rule of thumb from Cohen d.
#' 
#' @seealso 
#' \code{\link{es_vargha_delaney_a}}, to determine a Vargha-Delaney A
#' 
#' \code{\link{es_convert}}, to convert this to Cohen d
#' 
#' \code{\link{th_cohen_d}}, rules of thumb for Cohen d
#' 
#' @references 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101–132. doi:10.3102/10769986025002101
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_vda <- function(a, qual="vd"){
  
  if (qual=="vd"){
    #Use Vargha and Delaney (2000, p. 106)
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
  else if (qual %in% c("sawilowsky", "cohen-conv", "lovakov", "rosenthal")){
    if (qual=="cohen-conv"){
      qual="cohen"}
    
    #convert to Cohen's d
    rb = es_convert(a, fr="vda", to="rb")
    d = es_convert(rb, fr="rb", to="cohend")
    
    res = th_cohen_d(d, qual)
    qual = res$classification
    ref = res$reference
    
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}