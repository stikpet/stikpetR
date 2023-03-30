#' Paired Samples Sign Test
#' 
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @param mu the difference according to the null hypothesis (default is 0)
#' @returns 
#' A dataframe with:
#' \item{pos}{the number of scores with a positive difference}
#' \item{neg}{the number of scores with a negative difference}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details
#' 
#' The formula used:
#' \deqn{sig. = 2\times \text{Bin}\left(n, \min\left(n_{pos}, n_{neg}\right), \frac{1}{2}\right)}
#' With:
#' \deqn{n_{pos} = \sum_{i=1}^n \begin{cases} 1 & \text{ if } d_i>d_{H0} \\ 0 & \text{ if } d_i \le d_{H0} \end{cases}}
#' \deqn{n_{neg} = \sum_{i=1}^n \begin{cases} 0 & \text{ if } d_i \ge d_{H0} \\ 1 & \text{ if } d_i < d_{H0} \end{cases}}
#' \deqn{d_i = x_i - y_i}
#' 
#' *Symbols used:*
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
#' @examples 
#' ord1 = c(5, 3, 3, 4, 3, 4, 3, 4, 4, 4, 5, 3, 1, 3, 2)
#' ord2 = c(1, 3, 3, 3, 3, 3, 2, 1, 3, 4, 5, 3, 2, 5, 2)
#' 
#' ts_sign_ps(ord1, ord2)
#' ts_sign_ps(ord1, ord2, mu=1)
#' 
#' @references 
#' Arbuthnott, J. (1710). An argument for divine providence, taken from the constant regularity observ’d in the births of both sexes. *Philosophical Transactions of the Royal Society of London, 27*(328), 186–190. https://doi.org/10.1098/rstl.1710.0011
#' 
#' Dixon, W. J., & Mood, A. M. (1946). The statistical sign test. *Journal of the American Statistical Association, 41*(236), 557–566. https://doi.org/10.1080/01621459.1946.10501898
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' 
#' @export
ts_sign_ps <- function(ord1, ord2, mu=0){

  datFrame = na.omit(data.frame(ord1,ord2))
  
  d = datFrame$ord1 - datFrame$ord2
  
  pos = sum(d>mu)
  neg = sum(d<mu)
  ties = sum(d==mu)
  
  nAdj = pos + neg
  
  k = min(pos, neg)
  
  pValue = 2*pbinom(k,nAdj,0.5) 
  
  if (pValue>1) {
    pValue = 1
  }

  results = data.frame(pos, neg, pValue)
  
  return(results)
  
}