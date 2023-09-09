#' (Pearson) Contingency Coefficient
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param adj c(NULL, "sakoda") adjustment to use (see details)
#' @param r the number of rows, required if adj="sakado"
#' @param c the number of columns, required if adj="sakado"
#' @return value of (Pearson) Contingency Coefficient
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' es_cont_coeff(chi2Value, n)
#' 
#' @details 
#' The formula used is (Pearson, 1904, p. 9):
#' \deqn{C = \sqrt\frac{\chi^{2}}{n + \chi^{2}}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{\chi^{2}} the chi-square test statistic
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies
#' }
#' 
#' The maximum value for C would be (Sakoda, 1977, p. 778):
#' \deqn{C_{max} = \sqrt{\frac{m - 1}{m}}}
#' Where \eqn{m} is the minimum of the number of rows, or number of columns.
#' 
#' Sakoda propses to divide the contingency coefficient by this maximum, i.e.
#' \deqn{\frac{C}{C_{max}}}
#' 
#' For a 2x2 table Cole C1 will also divide by the maximum and produce the same result.
#' 
#' **Alternative**
#' 
#' The *'DescTools'* library also has a function for this: *ContCoef()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Pearson, K. (1904). *Contributions to the Mathematical Theory of Evolution. XIII. On the theory of contingency and its relation to association and normal correlation*. Dulau and Co.
#' 
#' Sakoda, J. M. (1977). Measures of Association for Multivariate Contingency Tables. *In Proceedings of the Social Statistics Section of the American Statistical Association: Vol. Part III* (pp. 777–780).
#' 
#' @export
es_cont_coeff <- function(chi2, n, adj=NULL, r=NULL, c=NULL){
  
  es = sqrt(chi2/(n + chi2))
  
  if (!is.null(adj) && adj=="sakoda") {
    m = r
    if (c < r){m = c}
    cmax = sqrt((m - 1)/m)
    es = es/cmax
  }
  
  return (es)
  
}