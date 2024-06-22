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
#' Blaikie-Roberts suggest to use as \eqn{C_{max}} (Blaikie, 1969, p.19):
#' \deqn{C_{max} = \sqrt[4]{\frac{r - 1}{r}\times\frac{c - 1}{c}}}
#' 
#' Blaikie refers to his mentor Roberts for this (Blaikie, 2003, p. 115)
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
#' Blaikie, N. W. H. (1969). Religion, social status, and community involvement: A study in Christchurch. *The Australian and New Zealand Journal of Sociology, 5*(1), 14–31. doi:10.1177/144078336900500102
#' 
#' Blaikie, N. W. H. (2003). *Analyzing quantitative data: From description to explanation*. Sage Publications Ltd.
#' 
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
  else if (!is.null(adj) && adj=="br") {
    cmax = ((r - 1)/r * (c - 1)/c)**0.25
    es = es/cmax
  }
  
  return (es)
  
}