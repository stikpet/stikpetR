#' Cole C5
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return Cole C5
#' 
#' @details 
#' The formula used (Cole, 1949, p. 416):
#' \deqn{C_5 = \frac{\sqrt{2}\times\left(a\times d - b\times c\right)}{\sqrt{\left(a\times d - b\times c\right)^2 + R1\times R2\times C1\times C2}}}
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' This is actually an adjustment of the contingency coefficient, by dividing it over \eqn{\sqrt{\frac{1}{2}}}, which is the maximum value in 2x2 tables.
#' 
#' Note that Cole C2 is the phi coefficient, Cole C3 is McEwen-Michael coefficient, Cole C4 is Yule Q, and Cole C6 is Yule r.
#' 
#' @author 
#' P. Stikker
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_cole_c5(bin1, bin2)
#' 
#' @export
es_cole_c5 <- function(var1, var2){
  
  data = data.frame(var1, var2)
  
  #remove missing values
  data = na.omit(data)
  
  #Create a cross table first
  ct = table(data)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #the row totals
  rowTots <- margin.table(ct, 1)
  R1 <- unname(rowTots[1])
  R2 <- unname(rowTots[2])
  
  #the column totals
  colTots <- margin.table(ct, 2)
  C1 <- unname(colTots[1])
  C2 <- unname(colTots[2])
  
  c5 = sqrt(2) * (a*d-b*c)/sqrt((a*d - b*c)**2 + R1*R2*C1*C2)
  
  return(c5)
}