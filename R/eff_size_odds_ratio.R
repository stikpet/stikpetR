#' Odds Ratio
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return dataframe with the odds ratio value, z-statistic, and 2-sided p-value
#' 
#' @details 
#' 
#' The p-value is for the null-hypothesis that the population OR is 1.
#' 
#' \deqn{OR = \frac{a\times d}{b \times c}}
#' \deqn{sig. = 2\times\left(1 - Z\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}}
#' \deqn{z = \frac{\ln{\left(OR\right)}}{SE}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' \item \eqn{Z\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' I could not find an original source for this, but one source is McHugh (2009).
#' The formula for the standard error can be found on page 123.
#' 
#' For a classification/qualification of the odds ratio use: *th_odds_ratio()*.
#' 
#' Or convert it to Cohen's d using either: *es_convert(or, from="or", to="cohend", ex1="chinn")* or *es_convert(or, from="or", to="cohend", ex1="borenstein")*
#' 
#' Then use classification/qualification of Cohen d: *th_cohen_d()*
#' 
#' Or convert it to Yule Q using: *es_convert(or, from="or", to="yuleq")*
#' 
#' Then use classification/qualification of Yule Q: *th_yule_q()*
#' 
#' Or convert it to Yule Y using: *es_convert(or, from="or", to="yuley")*
#' 
#' **Alternative**
#' 
#' R's *stats* library has a function that also shows an odds ratio: *fisher.test()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' McHugh, M. (2009). The odds ratio: Calculation, usage, and interpretation. *Biochemia Medica, 19*(2), 120–126. https://doi.org/10.11613/BM.2009.011
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_odds_ratio(bin1, bin2)
#' 
#' @export
es_odds_ratio <- function(var1, var2){
  
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
  
  #The odds ratio:
  or = (a/c) / (b/d)
  
  #Significance
  L = log(or)
  SE = sqrt(sum(1/ct))
  z = L/SE
  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic = z
  results = data.frame(or, statistic, pValue)
  
  return (results)
}