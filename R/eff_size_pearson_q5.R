#' Pearson Q5
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return Pearson Q4
#' 
#' @details
#' This is an approximation for a tetrachoric correlation coefficient.
#' 
#' The formula used (Pearson, 1900, p. 16):
#' \deqn{Q_5 = \sin\left(\frac{\pi}{2} \times \frac{1}{\sqrt{1+ k^2}}\right)}
#' with:
#' \deqn{k^2 = \frac{4\times a\times b\times c\times d\times n^2}{\left(a\times d-b\times c\right)^2\times\left(a+d\right)\times\left(b+c\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{n} the sum of all counts
#' }
#' 
#' Note that Pearson Q2 is the same as Yule Q and Pearson Q3 is the same as Yule r.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London, 195*, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_pearson_q5(bin1, bin2)
#' 
#' @export
es_pearson_q5 <- function(var1, var2){
  
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
  
  n = a+b+c+d
  
  k = 4*a*b*c*d*n**2 / ((a*d-b*c)**2 * (a+d)*(b+c))
  q5 = sin(pi/2 * 1/(sqrt(1+k)))
  
  return(q5)
  
}