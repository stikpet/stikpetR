#' One-Sample Score Test
#' 
#' @param data A vector with the data
#' @param codes optional vector with the two codes to use
#' @param p0 optional the hypothesized proportion for the first category (default is 0.5)
#' @param cc Use of continuity correction (default is "none")
#' @returns 
#' Dataframe with:
#' \item{statistic}{the test value}
#' \item{pValue}{two-sided p-value}
#' \item{testUsed}{a description of the test used}
#' 
#' @description 
#' A one-sample score test could be used with binary data, to test if the two categories
#' have a significantly different proportion. It is an approximation of a binomial test, 
#' by using a standard normal distribution. Since the binomial distribution is discrete 
#' while the normal is continuous, a so-called continuity correction can (should?) be 
#' applied.
#' 
#' The null hypothesis is usually that the proportions of the two categories in the 
#' population are equal (i.e. 0.5 for each). If the p-value of the test is below the 
#' pre-defined alpha level (usually 5% = 0.05) the null hypothesis is rejected and 
#' the two categories differ in proportion significantly.
#' 
#' The input for the function doesn't have to be a binary variable. 
#' A nominal variable can also be used and the two categories to compare indicated.
#' 
#' A significance in general is the probability of a result as in the sample, 
#' or more extreme, if the null hypothesis is true. 
#' 
#' Some info on the different tests can be found in \href{https://youtu.be/jQ-nSPTGOgE}{video}.
#' 
#' @details 
#' Also sometimes called a 'proportion' test.
#' The formula used is (Wilson, 1927):
#' \deqn{z=\frac{x - \mu}{SE}}
#' With:
#' \deqn{\mu = n\times p_0}
#' \deqn{SE=\sqrt{\mu\times\left(1-p_0\right)}}
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
#' The formula used and naming comes from IBM (2021, p. 997) who refer to Agresti, most likeli Agresti (2013, p. 10)
#' 
#' @section Alternatives:
#' 
#' R's *stats* library: prop.test().
#' 
#' @examples 
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male")
#' ts_score_os(data, c("Female", "Male"), p0 = 0.5)
#' 
#' @seealso 
#' Effect size measures that could go with the test are Cohen g (\code{\link{es_cohen_g}}), 
#' Cohen h (\code{\link{es_cohen_h_os}}), or the alternative ratio (\code{\link{es_alt_ratio}})
#' 
#' Other tests for a binary variable are the binomial test (\code{\link{ts_binomial_os}}), 
#' and Wald test (\code{\link{ts_wald_os}})
#' 
#' @references 
#' Agresti, A. (2013). *Categorical data analysis* (3rd ed.). Wiley.
#' 
#' IBM SPSS Statistics Algorithms. (2021). IBM.
#' 
#' Wilson, E. B. (1927). Probable Inference, the Law of Succession, and Statistical Inference. *Journal of the American Statistical Association, 22*(158), 209–212. https://doi.org/10.2307/2276774
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_score_os <- function(data, codes=NULL, p0 = 0.5, cc = c("none", "yates")){
  
  if (length(cc)>1) {cc = "none"}
  
  #if no codes provided use first found
  if (is.null(codes)) {
    n1 = unname(table(data)[1])
    n2 = sum(table(data)) - n1
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
  }
  
  #total sample size
  n = n1 + n2
  
  #use minimum of the two
  minCount = n1
  ExpProp = p0
  if (n2 < n1){
    minCount = n2
    ExpProp = 1 - ExpProp}
    
  #Normal approximation
  if (cc=="none"){
    p = minCount / n
    q = 1 - p
    se = (p0 * (1 - p0) / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-Sample Score"}
  
  else if (cc == "yates"){
    #Normal approximation with continuity correction
    p = (minCount + 0.5) / n
    q = 1 - p
    se = (p0 * (1 - p0) / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-Sample Score with Yates continuity correction"}
  
  testResults <- data.frame(statistic, pValue, testUsed)
  
  return (testResults)
}
