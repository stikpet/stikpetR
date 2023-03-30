#' Yule Q
#' 
#' @param var1 A vector with the data from the first variable
#' @param var2 A vector with the data from the second variable
#' @return Yule Q
#' 
#' @details
#' Yule's Q (1900) is also Cole's C4 (1949, p. 415) and Pearson's Q2 (1900, p. 15).
#' 
#' The formula used (Yule, 1900, p. 272):
#' \deqn{Q = \frac{a\times d - b\times c}{a\times d + b\times c}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' }
#' 
#' For a rule-of-thumb on the classification use *th_yule_q()*, 
#' or convert the Yule Q to an Odds Ratio using *es_convert(q, "yuleq", "or")* and 
#' use a rule of thumb for odds ratios (*th_odds_ratio()*).
#' 
#' Yule Q can also be converted to Yule Y using: *es_convert(q, "yuleq", "yuley")*
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
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London*, 195, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' Yule, G. U. (1900). On the association of attributes in statistics: With illustrations from the material of the childhood society, &c. *Philosophical Transactions of the Royal Society of London*, 194, 257–319. https://doi.org/10.1098/rsta.1900.0019
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_yule_q(bin1, bin2)
#' 
#' @export
es_yule_q <- function(var1, var2){
  
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
  
  q = (a*d - b*c)/(a*d + b*c)
  
  return(q)
  
}