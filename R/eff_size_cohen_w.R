#' Cohen's w
#' 
#' @description 
#' An effect size measure that could be used with a chi-square test. It has no upper limit, but can be compared to 
#' Cohen's rules-of-thumb. 
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' 
#' @return value of Cohen's w
#' 
#' @details 
#' The formula used is (Cohen, 1988, p. 216):
#' \deqn{w = \sqrt\frac{\chi_{GoF}^{2}}{n}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{\chi_{GoF}^{2}} the Pearson chi-square goodness-of-fit value 
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies
#' }
#' 
#' @section Alternative:
#' 
#' The *'rcompanion'* library also has a function for this: *cohenW()*
#' 
#' @seealso 
#' \code{\link{ts_pearson_gof}}, to obtain a chi-square value
#' 
#' \code{\link{th_cohen_w}}, rules-of-thumb for Cohen w
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' es_cohen_w(chi2Value, n)
#' 
#' @export
es_cohen_w <-function(chi2, n){
  #this function calculates Cohens w
  w = sqrt(chi2 / n)
  
  return (w)
}
