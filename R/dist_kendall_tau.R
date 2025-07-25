#' Kendall Tau Distribution
#' 
#' @param n the sample size (number of pairs)
#' @param tau Kendall tau value
#' @param method algorithm to use
#' @returns 
#' \item{pValue}{the two-tailed significance (p-value)}
#' 
#' @details 
#' If *method="AS71"* Algorithm AS 71 (Best & Gipps, 1974) will be used, 
#' by running the helper function *he_AS71(S, n)*.
#' The test statistic is:
#' \deqn{S = \binom{n}{2}\times\left|\tau\right| = \frac{n\times\left(n - 1\right)}{2}\times\left|\tau\right|}
#' AS 71 returns upper values only, so they get doubled for a two-sided test.
#' 
#' If *method="kendall"* the algorithm found at https://github.com/scipy/scipy/blob/v1.10.1/scipy/stats/_mstats_basic.py#L774-L898
#' was adapted. This refers to Kendall (1970), and uses the helper function *he_kendall(n, C)*.
#' Where \eqn{C = n_c}, i.e. the number of concordant pairs.
#' This algorithm already returns a two-tailed result.
#' 
#' @references 
#' Best, D. J., & Gipps, P. G. (1974). Algorithm AS 71: The upper tail probabilities of Kendall's tau. *Applied Statistics, 23*(1), 98-100. https://doi.org/10.2307/2347062
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
di_kendall_tau <- function(n, tau, method=c("kendall", "AS71")){
  
  if (length(method)>1) {
    method="kendall"
  }
  
  if (method == "kendall") {
    C = (tau + 1)*(n*(n - 1))/4
    pValue = he_kendall(n, C)
  }
  
  else if (method == "AS71") {
    S = choose(n, 2)*abs(tau)
    pValue = 2*he_AS71(S, n)
  }
  
  return(pValue)
  
}



