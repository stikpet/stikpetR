#' Rule-of-Thumb for Rank Biserial Correlation
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Rank Biserial Correlation.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/31jr7jUCni4) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Correlations/RankBiserialCorrelation.html)
#' 
#' 
#' @param rb the rank-biserial correlation value
#' @param qual optional setting for which rule of thumb to use. Either "cohen" (default), "vd", "sawilowsky", "cohen-conv", "lovakov", "rosenthal", "brydges"
#' 
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' 
#' @details 
#' Cohen's rule of thumb for rank-biserial correlation (1988, p. 82):
#' 
#' |\|r_b\|| Interpretation|
#' |---|----------|
#' |0.000 < 0.125 | negligible |
#' |0.125 < 0.304 | small |
#' |0.304 < 0.465 | medium |
#' |0.465 or more | large |
#' 
#' Vargha and Delaney (2000, p. 106):
#' 
#' |\|r_b\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.11 | negligible |
#' |0.11 < 0.28 | small |
#' |0.28 < 0.43 | medium |
#' |0.43 or more | large |
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to obtain the measure:
#' \code{\link{r_rank_biserial_is}}, to determine a the rank biserial for independent samples.
#' \code{\link{r_rank_biserial_os}}, to determine a the rank biserial for one-sample.
#' 
#' The function uses the convert function and corresponding rules of thumb:
#' \code{\link{es_convert}}, to convert this to Cohen d.
#' \code{\link{th_cohen_d}}, rules of thumb for Cohen d.
#' 
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples
#' # Example 1: using Cohen's rules:
#' rb = 0.6
#' th_rank_biserial(rb)
#' 
#' # Example 2: Convert to Cohen d, then use Cohen d rules:
#' rb= 0.23
#' th_rank_biserial(rb, qual="cohen-conv")
#' 
#' @export
th_rank_biserial <- function(rb, qual="cohen"){
  
  if (qual=="cohen"){
    #Use Cohen (1988, p. 82)
    ref = "Cohen (1988, p. 82)"
    if (abs(rb) < 0.125){
      qual = "negligible"}
    else if (abs(rb) < 0.304){
      qual = "small"}
    else if (abs(rb) < 0.465){
      qual = "medium"}
    else{
      qual = "large"}
  }
  else if (qual=="vd"){
    #Use Vargha and Delaney (2000, p. 106)
    ref = "Vargha and Delaney (2000, p. 106)"
    if (abs(rb) < 0.11){
      qual = "negligible"}
    else if (abs(rb) < 0.28){
      qual = "small"}
    else if (abs(rb) < 0.43){
      qual = "medium"}
    else{
      qual = "large"}
  }
  else if (qual %in% c("sawilowsky", "cohen-conv", "lovakov", "rosenthal", "brydges")){
    if (qual=="cohen-conv"){
      qual="cohen"}
    
    #convert to Cohen's d
    d = es_convert(rb, fr="rb", to="cohend")
    
    res = th_cohen_d(d, qual)
    qual = res$classification
    ref = res$reference
    
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}



