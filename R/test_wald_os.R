#' One-Sample Wald Test
#' 
#' @param data A vector with the data
#' @param codes Optional vector with the two codes to use
#' @param p0 The hypothesized proportion for the first category (default is 0.5)
#' @param cc c(NULL, "yates") use of continuity correction (default is NULL)
#' @return dataframe with test-value, two-sided p-value and test used
#' 
#' @examples 
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Female","Male", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Female")
#' ts_wald_os(data, c("Female", "Male"), p0 = 0.5)
#' 
#' @details 
#' This test differs from the one-sample score test in the calculation of the standard error. 
#' For the ‘regular’ version this is based on the expected proportion, 
#' while for the Wald version it is done with the observed proportion.
#' 
#' The formula used (Wald, 1943):
#' \deqn{z=\frac{x - \mu}{SE}}
#' With:
#' \deqn{\mu = n\times p_0}
#' \deqn{SE = \sqrt{x\times\left(1 - \frac{x}{n}\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x} is the number of successes in the sample
#' \item \eqn{p_0} the expected proportion (i.e. the proportion according to the null hypothesis)
#' }
#' 
#' If the Yates continuity correction is used the formula changes to (Yates, 1934, p. 222):
#' \deqn{z_{Yates} = \frac{\left|x - \mu\right| - 0.5}{SE}}
#' 
#' The formula used in the calculation is the one from IBM (2021, p. 997).
#' IBM refers to Agresti, most likely Agresti (2013, p. 10), who in turn
#' refer to Wald (1943)
#' 
#' **Alternative**
#' 
#' No alternative library that has this function exist to my knowledge.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @references 
#' Agresti, A. (2013). *Categorical data analysis* (3rd ed.). Wiley.
#' 
#' IBM SPSS Statistics Algorithms. (2021). IBM.
#' 
#' Wald, A. (1943). Tests of statistical hypotheses concerning several parameters when the number of observations is large. *Transactions of the American Mathematical Society, 54*(3), 426–482. https://doi.org/10.2307/1990256
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' @export
ts_wald_os <- function(data, codes=NULL, p0=0.5, cc=NULL){
  
  #one-sample Wald test
  #an approximation using the normal distribution for a one-sample binomial test
  
  #if no codes provided use first found
  if (is.null(codes)) {
    n1 = unname(table(data)[1])
    n2 = sum(table(data)) - n1
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
  }
  n = n1 + n2
  
  minCount = n1
  ExpProp = p0
  if (n2 < n1){
    minCount = n2
    ExpProp = 1 - ExpProp}
  
  #Wald approximation
  if (is.null(cc)){
    p = minCount / n
    q = 1 - p
    se = (p * q / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-sample Wald"}
  
  else if (cc == "yates"){
    #Wald approximation with continuity correction
    p = (minCount + 0.5) / n
    q = 1 - p
    se = (p * q / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-sample Wald with Yates continuity correction"}
  
  testResults <- data.frame(statistic, pValue, testUsed)
  
  return (testResults)
}