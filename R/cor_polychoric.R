#' Polychoric Correlation Coefficient
#' 
#' @param dataVar A vector with the scores data
#' @param groupVar A vector with the group data
#' @return Polychoric Correlation Coefficient value
#' 
#' @details
#' This simply uses the *polychor()* function from the *polycor* library 
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @examples 
#' scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
#' groups = c("A","A","A","B","B","B","B", NA, "C")
#' r_polychoric(scores, groups)
#' 
#' @export
r_polychoric <- function(dataVar, groupVar){
  
  #make sure data is numeric
  scores = as.numeric(dataVar)
  
  #remove rows with missing values
  df = data.frame(dataVar, groupVar)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  rp = polycor::polychor(table(df))
  
  return(rp)
  
}