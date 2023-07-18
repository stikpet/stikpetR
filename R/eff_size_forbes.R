#' Forbes Coefficient
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Forbes Coefficient
#' 
#' @details 
#' 
#' The formula used (Forbes, 1907, p. 279):
#' \deqn{F = \frac{n\times\min\left(a, d\right)}{C_1\times R_1}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_1} the sum of counts in the 1st row 
#' \item \eqn{C_1} the sum of counts in the 1st column 
#' }
#' 
#' The coefficient has a value of 1 if there is no association, while it has a value of 0 or 2 when there is a perfect one.
#' To adjust to the more traditional range of -1 to 1, Cole 1 simply subtracts one from the Forbes coefficient.
#' 
#' Note that Alroy added an adjustment to Forbes F, see *es_alroy_f()* for this adjusted version.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Forbes, S. A. (1907). On the local distribution of certain Illinois fishes: An essay in statistical ecology. *Illinois Natural History Survey Bulletin, 7*(8), 273–303.
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_forbes(bin1, bin2)
#' 
#' @export
es_forbes <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #the row totals
  rowTots <- margin.table(ct, 1)
  R1 <- unname(rowTots[1])
  
  #the column totals
  colTots <- margin.table(ct, 2)
  C1 <- unname(colTots[1])
  
  #grand total
  n <- sum(colTots)
  
  forb = n*a/(C1 * R1)
  
  return(forb)
}