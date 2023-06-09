#' Spearman Rho / Rank Correlation Coefficient
#' 
#' @description 
#' The Spearman Rank Correlation Coefficient is the Pearson Correlation Coefficient, 
#' after the scores first have been converted to ranks.
#' 
#' This function makes use of *di_spearman()* for the test of this correlation,
#' which requires the *pspearman* library for exact computations.
#' 
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @param method the test to be used
#' @param cc boolean to indicate the use of a continuity correction
#' @param iters the number of iterations to use, only applicable if Iman-Conover is used
#' @returns 
#' A dataframe with:
#' \item{rs}{the correlation coefficient}
#' \item{statistic}{the statistic from the test (only if applicable)}
#' \item{df}{the degrees of freedom (only if applicable)}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used is (Spearman, 1904, p. 77):
#' \deqn{r_s = \frac{SS_{r_x, r_y}}{SS_{r_x}\times SS_{r_y}}}
#' With:
#' \deqn{SS_{r_x} = \sum_{i=1}^n \left(r_{x_i} - \bar{r}_x\right)^2}
#' \deqn{SS_{r_y} = \sum_{i=1}^n \left(r_{y_i} - \bar{r}_y\right)^2}
#' \deqn{SS_{r_x, r_y} = \sum_{i=1}^n \left(r_{x_i} - \bar{r}_x\right) \times \left(r_{y_i} - \bar{r}_y\right)}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{r_{x_i}} the i-th rank of the scores of the first variable
#' \item \eqn{r_{y_i}} the i-th rank of the scores of the second variable
#' \item \eqn{n} the total sample size (number of ranks)
#' }
#' 
#' If all the ranks are distinct (i.e. no ties) the formula can also be written as:
#' 
#' \deqn{r_s = 1 - \frac{6}{n\times\left(n^2 - 1\right)}\times S}
#' With:
#' \deqn{S = \sum_{i=1}^n d_i^2}
#' \deqn{d_i^2 = \left(r_{x_i} - r_{y_i}\right)^2}
#' 
#' The test can be performed in different ways. Options to choose from are:
#' \itemize{
#' \item "t" uses a Student t distribution approximation
#' \item "z-fieller" uses a standard normal approximation from Fieller
#' \item "z-old" uses standard normal approximation from Old
#' \item "iman-conover" a combination of z and t distribution from Iman and Conover
#' \item "AS89" uses the AS 89 algorithm
#' \item "exact" uses an exact distribution
#' }
#' 
#' See for the details of each the *di_spearman()* function
#' 
#' A continuity correction can be applied (Zar, 1972, p. 579):
#' \deqn{r_s^{cc} = \left|r_s\right| - \frac{6}{n^3 - n}}
#' 
#' **Alternatives**
#' 
#' *R's stats*
#' 
#' Using the t-approximation:
#' 
#' cor.test(ord1, ord2, method="spearman")
#' 
#' Using AS89
#' 
#' cor.test(ord1, ord2, method="spearman", exact=TRUE)
#' 
#' *library(pspearman)*
#' 
#' spearman.test(ord1, ord2, approximation="t-distribution")
#' 
#' spearman.test(ord1, ord2, approximation="AS89")
#' 
#' spearman.test(ord1, ord2, approximation="exact")
#' 
#' @references 
#' Spearman, C. (1904). The proof and measurement of association between two things. *The American Journal of Psychology, 15*(1), 72–101.
#' 
#' Zar, J. H. (1972). Significance testing of the Spearman rank correlation coefficient. *Journal of the American Statistical Association, 67*(339), 578–580. https://doi.org/10.1080/01621459.1972.10481251
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' ord1 = c(5, 3, 3, 4, 3, 4, 3, 4, 4, 4, 5, 3, 1, 3, 2)
#' ord2 = c(5, 3, 3, 3, 3, 3, 5, 4, 3, 4, 5, 3, 2, 5, 2)
#' 
#' r_spearman_rho(ord1, ord2, method="t")
#' r_spearman_rho(ord1, ord2, method="z-fieller")
#' r_spearman_rho(ord1, ord2, method="z-olds")
#' r_spearman_rho(ord1, ord2, method="iman-conover")
#' r_spearman_rho(ord1, ord2, method="AS89")
#' r_spearman_rho(ord1, ord2, method="exact")
#' 
#' r_spearman_rho(ord1, ord2, method="t", cc=TRUE)
#' r_spearman_rho(ord1, ord2, method="z-fieller", cc=TRUE)
#' r_spearman_rho(ord1, ord2, method="z-olds", cc=TRUE)
#' r_spearman_rho(ord1, ord2, method="iman-conover", cc=TRUE)
#' r_spearman_rho(ord1, ord2, method="AS89", cc=TRUE)
#' 
#' @export
r_spearman_rho <- function(ord1, ord2, 
                           method=c("t", "z-fieller", "z-olds", "iman-conover", "AS89", "exact"), 
                           cc=FALSE, iters=500){
  
  if (length(method)>1) {
    method = "t"
  }
  
  dFr = na.omit(data.frame(ord1, ord2))
  
  rx = rank(dFr$ord1)
  ry = rank(dFr$ord2)
  
  mrx = mean(rx)
  mry = mean(ry)
  
  n = length(rx)
  
  SSrx = var(rx)*(n - 1)
  SSry = var(ry)*(n - 1)
  
  SSxy = sum((rx - mrx)*(ry - mry))
  
  rs = SSxy / sqrt(SSrx*SSry)
  
  if (cc && method!="exact") {
    #(Zar, 1972, p. 579)
    rs = abs(rs) - 6/(n**3 - n)
  }
  
  df = n - 2
  
  if (method=="iman-conover") {
    distRes = di_spearman(n, rs, method=method, iters=iters)
  }
  else{
    distRes = di_spearman(n, rs, method=method)
  }
  
  pValue = distRes$pValue
  
  if (method=="t" || method=="iman-conover" || method=="AS89") {
    statistic = distRes$statistic
    results = data.frame(rs, statistic, df, pValue)
  }
  else if (method=="z-fieller" || method=="z-olds"){
    statistic = distRes$statistic
    results = data.frame(rs, statistic, pValue)
  }
  
  else if (method=="exact") {
    results = data.frame(pValue)
  }
  return(results)
}