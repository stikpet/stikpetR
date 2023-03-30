#' Hedges g (for one-sample)
#' 
#' @param data pandas series with the numeric scores
#' @param mu optional parameter to set the hypothesized mean. If not used the midrange is used
#' @param appr c(NULL, 'hedges', 'durlak', 'xue') approximation to use, NULL will use exact.
#' @return dataframe with the hypothesized mean used, the effect size value, and method used.
#' 
#' @examples
#' data <- c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, 90, 25, 40, 70, 65, 70, 98, 40, 65, 60)
#' es_hedges_g_os(data, mu=50)
#'  
#' @details 
#' Hedges g is a correction for Cohen's d'. Actually Hedges (1981) didn't seem to have
#' a one-sample version for Hedges g, and this correction is the one for Hedges g used for the
#' independent samples.
#' 
#' The formula used for the exact method (appr=NULL) (Hedges, 1981, p. 111):
#' \deqn{g = d' \times\frac{\Gamma\left(m\right)}{\Gamma\left(m - \frac{1}{2}\right)\times\sqrt{m}}}
#' With:
#' \deqn{m = \frac{df}{2}}
#' \deqn{df = n - 1}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{d'} Cohen’s d for one-sample
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{\Gamma\left(\dots\right)} the gamma function
#' }
#' 
#' The formula used for the approximation from Hedges (1981, p. 114) (appr="hedges"):
#' \deqn{g = d' \times\left(1 - \frac{3}{4\times df - 1}\right)}
#' 
#' The formula used for the approximation from Durlak (2009, p. 927) (appr="durlak"):
#' \deqn{g = d' \times\frac{n - 3}{n - 2.25} \times\sqrt{\frac{n - 2}{n}}}
#' 
#' The formula used for the approximation from Xue (2020, p. 3) (appr="xue"):
#' \deqn{g = d' \times \sqrt[12]{1 - \frac{9}{df} + \frac{69}{2\times df^2} - \frac{72}{df^3} + \frac{687}{8\times df^4} - \frac{441}{8\times df^5} + \frac{247}{16\times df^6}}}
#' 
#' **Alternative**
#' 
#' The *effectsize* library has a similar function: *hedges_g()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references
#' Durlak, J. A. (2009). How to select, calculate, and interpret effect sizes. *Journal of Pediatric Psychology, 34*(9), 917–928. https://doi.org/10.1093/jpepsy/jsp004
#' 
#' Hedges, L. V. (1981). Distribution Theory for Glass’s Estimator of Effect Size and Related Estimators. *Journal of Educational Statistics, 6*(2), 107–128. https://doi.org/10.2307/1164588
#' 
#' Xue, X. (2020). Improved approximations of Hedges’ g*. https://doi.org/10.48550/arXiv.2003.06675
#' 
#' @export
es_hedges_g_os <- function(data, mu=NULL, appr=NULL){
  # Hedges g for a one-sample test
  # Numeric data vector as input
  # mu: the hypothesized mean (default is "none")
  # appr: approximation to use. Can be set to 'hedges', 'durlak', or 'xue'. NULL will give the exact result.
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2}
  
  n = length(data)
  df = n - 1
  avg = mean(data)
  s = sd(data)
  
  d = (avg - mu)/s
  
  m = df/2
  
  if (is.null(appr) && m <172){
    # Hedges g (exact)
    g = d*gamma(m)/(gamma(m-0.5)*sqrt(m))
    comment = "exact"
  }
  else if(appr=="hedges"){
    # Hedges approximation
    g = d*(1-3/(4*df-1))
    comment = "Hedges approximation"
  }
  else if(appr=="durlak"){
    # Durlak (2009, p. 927) approximation:
    g = d*(n-3)/(n-2.25)*sqrt((n-2)/n)
    comment = "Durlak approximation"
  }
  else{
    # Xue (2020, p. 3) approximation:
    g = d*(1 - 9/df + 69/(2*df^2) - 72/(df^3) + 687/(8*df^4) - 441/(8*df^5) + 247/(16*df^6))^(1/12)
    comment = "Xue approximation"
  }
  
  if (is.null(appr) && m > 171){
    print("WARNING: exact method could not be computed due to large sample size, Xue approximation used instead")
  }
  
  #prepare results
  res <- data.frame(mu, g, comment)
  colnames(res) <- c('mu', 'g', 'version')
  
  return(res)
  
}