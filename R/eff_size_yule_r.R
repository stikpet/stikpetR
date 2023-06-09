#' Yule r
#' 
#' @param var1 A vector with the data from the first variable
#' @param var2 A vector with the data from the second variable
#' @return Yule r
#' 
#' @details 
#' This is an approximation for the tetrachoric correlation coefficient.
#' 
#' The formula used is the one from Cole (1949, p. 416):
#' \deqn{C_6 = \cos\left(\frac{\pi\times\sqrt{b\times c}}{\sqrt{a\times d} + \sqrt{b\times c}}\right)}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' This version can also be found in Walker and Lev (1953, p. 274)
#' 
#' Yule original formula (Yule, 1900, p. 276):
#' \deqn{r_{Yule} = \cos\left(\frac{\sqrt{k}}{1 + \sqrt{k}}\times\pi\right)}
#' With:
#' \deqn{k = \frac{1 - Q}{1 + Q}}
#' \deqn{Q = \frac{a\times d - b\times c}{a\times d + b\times c}}
#' 
#' Pearson original formula (Pearson, 1900, p. 16):
#' \deqn{Q_3 = \sin\left(\frac{\pi}{2}\times\frac{\sqrt{a\times d} - \sqrt{b\times c}}{\sqrt{a\times d} - \sqrt{b\times c}}\right)}
#' Although the original article has a small error in the equation (it divides the numerator and denominator by the same value)
#' 
#' Note that \eqn{r_{Yule} = Q_3 = C_6}
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
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London, 195*, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' Walker, H. M., & Lev, J. (1953). Statistical inference. Holt. https://catalog.hathitrust.org/Record/001306434
#' 
#' Yule, G. U. (1900). On the association of attributes in statistics: With illustrations from the material of the childhood society, &c. *Philosophical Transactions of the Royal Society of London, 194*, 257–319. https://doi.org/10.1098/rsta.1900.0019
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_yule_r(bin1, bin2)
#' 
#' @export
es_yule_r <- function(var1, var2){
  
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
  
  c6 = cos(pi*sqrt(b*c)/(sqrt(a*d)+sqrt(b*c)))
  
  return(c6)
}