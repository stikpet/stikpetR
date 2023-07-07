#' Bonett and Price Y*
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return Bonett and Price Y*
#' 
#' @details 
#' The formula used (Bonett and Price, 2007, pp. 433-434):
#' \deqn{Y* = \frac{\hat{\omega}^x-1}{\hat{\omega}^x+1}}
#' With:
#' \deqn{x = \frac{1}{2}-\left(\frac{1}{2}-p_{min}\right)^2}
#' \deqn{p_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)}{n}}
#' \deqn{\hat{\omega} = \frac{\left(a+0.1\right)\times\left(d+0.1\right)}{\left(b+0.1\right)\times\left(c+0.1\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' Note that \eqn{\hat{\omega}} is a biased corrected version of the Odds Ratio
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Bonett, D. G., & Price, R. M. (2007). Statistical inference for generalized yule coefficients in 2 × 2 contingency tables. *Sociological Methods & Research, 35*(3), 429–446. https://doi.org/10.1177/0049124106292358
#'
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female",
#' "female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#' "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#' "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl",
#' "nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", 
#' "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", 
#' "other", "other")
#' es_bonett_price_y(bin1, bin2)
#'
#' @export
es_bonett_price_y <- function(var1, var2){
  
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
  rowTots <- rowSums(ct)
  
  #the column totals
  colTots <- colSums(ct)

  #grand total
  n <- sum(colTots)
  
  #smallest proportion of rows and columns
  pMin = min(min(rowTots),  min(colTots))/n
  
  #determine power to use
  pwr = 1/2 - (1/2 - pMin)^2
  
  #odds ratio
  OR = a*d/(b*c)
  #adjust the odds ratio
  ORadj = (a+0.1)*(d+0.1)/((b+0.1)*(c+0.1))
  
  ybp = (ORadj^pwr - 1)/(ORadj^pwr + 1)
  
  return(ybp)
}