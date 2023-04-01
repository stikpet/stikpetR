#' Cramer's V for Goodness-of-Fit
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param k the number of categories
#' @param bergsma optional boolean to indicate the use of the Bergsma correction (default is False)
#' @return Cramer's V value
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' k <- 4
#' es_cramer_v_gof(chi2Value, n, k)
#' es_cramer_v_gof(chi2Value, n, k, TRUE)
#' 
#' @details
#' The formula used is:
#' \deqn{V=\sqrt\frac{\chi_{GoF}^{2}}{n\times \left(k - 1\right)}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{k} the number of categories
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies
#' \item \eqn{\chi_{GoF}^{2}} the chi-square value of a Goodness-of-Fit test
#' }
#' 
#' The Bergsma correction uses a different formula.
#' \deqn{\tilde{V} = \sqrt{\frac{\tilde{\varphi}^2}{\tilde{k} - 1}}}
#' With:
#' \deqn{\tilde{\varphi}^2 = max\left(0,\varphi^2 - \frac{k - 1}{n - 1}\right)}
#' \deqn{\tilde{k} = k - \frac{\left(k - 1\right)^2}{n - 1}}
#' \deqn{\varphi^2 = \frac{\chi_{GoF}^{2}}{n}}
#' 
#' Cramér described V (1946, p. 282) for use with a test of independence.
#' Others (e.g. K. Kelley & Preacher, 2012, p. 145; Mangiafico, 2016a, p. 474) 
#' added that this can also be use for goodness-of-fit tests.
#' For the Bergsma (2013, pp. 324-325) correction the same thing applies
#' 
#' Cramér's V can be converted to Cohen's w using *es_convert(from="cramervgof", to = "cohenw", ex1 = df)*
#' 
#' Rules-of-thumb for the interpretation can then be used, using *th_cohen_w(w)*
#' 
#' **Alternatives**
#' The *lsr* library has a similar function: *cramersV()*
#' 
#' The *DescTools* library has a similar function: *CramerV()*
#' 
#' @references 
#' Bergsma, W. (2013). A bias-correction for Cramér’s and Tschuprow’s. Journal of the Korean Statistical Society, 42(3), 323–328. https://doi.org/10.1016/j.jkss.2012.10.002
#' 
#' Cramér, H. (1946). Mathematical methods of statistics. Princeton University Press.
#' 
#' Kelley, K., & Preacher, K. J. (2012). On effect size. Psychological Methods, 17(2), 137–152. https://doi.org/10.1037/a0028086
#' 
#' Mangiafico, S. S. (2016). Summary and analysis of extension program evaluation in R (1.13.5). Rutger Cooperative Extension.
#' 
#' #' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
es_cramer_v_gof <- function(chi2, n, k, bergsma=FALSE){
  
  df <- k - 1
  
  if (bergsma){
    kAvg <- k - (k - 1)^2/(n - 1)
    phi2 <- chi2/n
    phi2Avg <- max(0, phi2 - (k - 1)/(n - 1))
    v <- sqrt(phi2Avg/(kAvg - 1))}
  else{
    v <- sqrt(chi2/(n * df))}

    
    return(v)
}