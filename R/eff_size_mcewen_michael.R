#' McEwen-Michael Coefficient / Cole C3
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return McEwen-Michael Coefficient
#' 
#' @details 
#' The formula used (Michael, 1920, p. 57)):
#' \deqn{\frac{a\times d - b\times c}{\left(\frac{a + d}{2}\right)^2 + \left(\frac{b + c}{2}\right)^2}}
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' Note that Cole (1949, p. 415) refers to this as C3.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' Michael, E. L. (1920). Marine Ecology and the coefficient of association: A plea in behalf of quantitative biology. *Journal of Ecology, 8*(1), 54–59. https://doi.org/10.2307/2255213
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_mcewen_michael(bin1, bin2)
#' 
#' @export
es_mcewen_michael <- function(var1, var2){
  
  data = data.frame(var1, var2)
  
  #remove missing values
  data = na.omit(data)
  
  #Create a cross table first
  dataTable = table(data)
  
  #store the individual cells
  a = dataTable[1,1]
  b = dataTable[1,2]
  c = dataTable[2,1]
  d = dataTable[2,2]
  
  mm = (a*d - b*c)/(((a+d)/2)**2 + ((b+c)/2)**2)
  
  return(mm)
}