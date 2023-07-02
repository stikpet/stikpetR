#' Rule-of-Thumb for Cohen w
#' 
#' @param w the Cohen w value
#' @param qual optional setting for which rule of thumb to use. Currently only 'cohen'
#' @return dataframe with the qualification and source
#' 
#' @examples 
#' w <- -0.2391304
#' th_cohen_w(w)
#' 
#' @details 
#' Cohen's rule of thumb for Cohen w (1988, p. 227):
#' 
#' |\|w\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.30 | small |
#' |0.30 < 0.50 | medium |
#' |0.50 or more | large |
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'  
#' @export
th_cohen_w <- function(w, qual="cohen"){
  
  #Use Cohen (1988, p. 227)
  ref = "Cohen (1988, p. 227)"
  if (abs(w) < 0.1){
    qual = "negligible"}
  else if (abs(w) < 0.3){
    qual = "small"}
  else if (abs(w) < 0.5){
    qual = "medium"}
  else{
    qual = "large"}
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
  
}


