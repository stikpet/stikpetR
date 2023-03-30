#' Edward Q
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return Edward Q
#' 
#' @details 
#' The formula used (Edwards 1957, as cited in Becker & Clogg, 1988, p. 408):
#' \deqn{Q_E = \frac{OR^{\frac{\pi}{4}} - 1}{OR^{\frac{\pi}{4}} + 1}}
#' With:
#' \deqn{OR = \frac{a\times d}{b \times c}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' \item \eqn{OR} the Odds Ratio
#' }
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Becker, M. P., & Clogg, C. C. (1988). A note on approximating correlations from Odds Ratios. *Sociological Methods & Research, 16*(3), 407–424. https://doi.org/10.1177/0049124188016003003
#' 
#' Edwards, J. H. (1957). A note on the practical interpretation of 2 x 2 tables. *Journal of Epidemiology & Community Health, 11*(2), 73–78. https://doi.org/10.1136/jech.11.2.73
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_edward_q(bin1, bin2)
#' 
#' @export
es_edward_q <- function(var1, var2){
  
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
  
  or = a*d/(b*c)
  q = ((or)**(pi/4)-1)/((or)**(pi/4)+1)
  
  return(q)
}