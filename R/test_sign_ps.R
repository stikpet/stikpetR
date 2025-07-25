#' Paired Samples Sign Test
#' @description 
#' This test compares the number of pairs that have a difference above the hypothesized difference, with those below the difference. It can be considered an alternative for the paired samples t-test.
#' 
#' @param field1 the numeric scores of the first variable
#' @param field2 the numeric scores of the second variable
#' @param levels vector, optional. the levels from field1 and field2
#' @param dmu float, optional. The difference according to the null hypothesis (default is 0)
#' @param method string, optional. Test to be used. Either "exact" (default), "appr".
#' 
#' @returns 
#' A dataframe with:
#' \item{n pos}{the number of scores with a positive difference}
#' \item{n neg}{the number of scores with a negative difference}
#' \item{statistic}{the test statistic (only applicable if method="appr")}
#' \item{p-Value}{the significance (p-value)}
#' 
#' @details 
#' If method="exact" the binomial distribution will be used. The formula used is (Dixon & Mood, 1946):
#' \deqn{sig. = 2\times \text{Bin}\left(n, \min\left(n_{pos}, n_{neg}\right), \frac{1}{2}\right)}
#' 
#' When using the approximation, the standard normal distribution is used (SPSS, 2006, p. 483):
#' \deqn{z = \frac{\max\left(n_{pos}, n_{neg}\right)-0.5\times\left(n_{pos} + n_{neg}\right)-0.5}{0.5\times\sqrt{n_{pos}+n_{neg}}}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \deqn{n_{pos} = \sum_{i=1}^n \begin{cases} 1 & \text{ if } d_i>d_{H0} \\ 0 & \text{ if } d_i \le d_{H0} \end{cases}}
#' \deqn{n_{neg} = \sum_{i=1}^n \begin{cases} 0 & \text{ if } d_i \ge d_{H0} \\ 1 & \text{ if } d_i < d_{H0} \end{cases}}
#' \deqn{d_i = x_i - y_i}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{n} is the number of pairs with a difference unequal to zero
#' \item \eqn{n_{pos}} the number of pairs with a positive difference
#' \item \eqn{n_{neg}} the number of pairs with a negative difference
#' \item \eqn{d_{H0}} the difference according to the null hypothesis, usually 0
#' \item \eqn{x_i} the i-th score from the first variable
#' \item \eqn{y_i} the i-th score from the second variable
#' \item \eqn{\text{Bin}\left(\dots,\dots\right)} the cumulative probability mass function of the binomial distribution
#' }
#' 
#' The test was described by Arbuthnott (1710) but with more modern notation see Dixon and Mood (1946).
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' SignTest(ord1, ord2)
#' 
#' *library(EnvStats)*
#' 
#' signTest(ord1, ord2, paired=TRUE)
#' 
#' *library(BSDA)*
#' 
#' SIGN.test(ord1, ord2)
#' 
#' @references 
#' Arbuthnott, J. (1710). An argument for divine providence, taken from the constant regularity observed in the births of both sexes. *Philosophical Transactions of the Royal Society of London, 27*(328), 186-190. doi:10.1098/rstl.1710.0011
#' 
#' Dixon, W. J., & Mood, A. M. (1946). The statistical sign test. *Journal of the American Statistical Association, 41*(236), 557-566. doi:10.1080/01621459.1946.10501898
#' 
#' SPSS. (2006). SPSS 15.0 algorithms.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_sign_ps <- function(field1, field2, levels=NULL, dmu=0, method="exact"){
  
  if (!is.null(levels)){
    myFieldOrd = factor(field1, ordered = TRUE, levels = levels)
    field1 = as.numeric(myFieldOrd)
  }
  
  if (!is.null(levels)){
    myFieldOrd = factor(field2, ordered = TRUE, levels = levels)
    field2 = as.numeric(myFieldOrd)
  }
  
  datFrame = na.omit(data.frame(field1,field2))
  
  d = datFrame$field1 - datFrame$field2
  
  pos = sum(d>dmu)
  neg = sum(d<dmu)
  ties = sum(d==dmu)
  
  nAdj = pos + neg
  
  k = min(pos, neg)
  
  if (method=="exact"){
    z = NA
    pValue = 2*pbinom(k,nAdj,0.5)} 
  else if (method=="appr"){
    if (k==neg){
      nMax = pos}
    else{
      nMax = neg}
    
    z = (nMax - 0.5 * nAdj - 0.5) / (0.5 * (nAdj)**0.5)
    pValue = 2 * (1 - pnorm(abs(z)))
  }
  
  
  if (pValue>1) {
    pValue = 1
  }
  
  results = data.frame(pos, neg, z, pValue)
  colnames(results) = c("n pos", "n neg", "statistic", "p-value")
  return(results)
  
}



