#' Yule's Y / Coefficient of Colligation
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Yule's Y
#' 
#' @details 
#' The formula used (Yule, 1912, p. 592):
#' \deqn{Y = \frac{\sqrt{a\times d} - \sqrt{b\times c}}{\sqrt{a\times d} + \sqrt{b\times c}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' Yule Y can be converted to Yule Q using: *es_convert(q, "yuley", "yuleq")*
#' 
#' Rules-of-thumb for Yule Q can then be used for classification/qualification: *th_yule_q(q)*
#' 
#' Yule Y can be converted to Odds Ratio using: *es_convert(q, "yuley", "or")*
#' 
#' Rules-of-thumb for Odds Ratio can then be used for classification/qualification: *th_odds_ratio(or)*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Yule, G. U. (1912). On the methods of measuring association between two attributes. *Journal of the Royal Statistical Society, 75*(6), 579–652. https://doi.org/10.2307/2340126
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_yule_y(bin1, bin2)
#' 
#' @export
es_yule_y <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #Yule's Y (1912, p. 592) is a further adaptation of Yule's Q into:
  y = ((a*d)**0.5 - (b*c)**0.5) / ((a*d)**0.5 + (b*c)**0.5)
  
  return(y)
}