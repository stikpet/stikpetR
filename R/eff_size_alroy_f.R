#' Alroy's Forbes Adjustment
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Alroy F
#' 
#' @details 
#' The formula used (Alroy, 2015):
#' \deqn{F' = \frac{a\times\left(n' + \sqrt{n'}\right)}{a\times\left(n' + \sqrt{n'}\right) + \frac{3}{2}\times b\times c}}
#' With:
#' \deqn{n' = a + b + c}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' Note Alroy refers to the Forbes coefficient as a measure of similarity.  
#' Alroy then sets out to improve the measure by disregarding the lower left value (d).
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Alroy, J. (2015). A new twist on a very old binary similarity coefficient. *Ecology, 96*(2), 575–586. https://doi.org/10.1890/14-0471.1
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", 
#' "female","female","female", "male", "male", "male", "male", "male", "male", "male", 
#' "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#' "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#' "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other",
#' "nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", 
#' "other", "other", "other", "other", "other", "other", "other", "other", "other", 
#' "other", "other", "other", "other", "other", "other")
#' es_alroy_f(bin1, bin2)
#' 
#' @export
es_alroy_f <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  nAdj = a + b + c
  
  fadj = a*(nAdj + sqrt(nAdj))/(a*(nAdj + sqrt(nAdj)) + 3/2*b*c)
  
  return(fadj)
}