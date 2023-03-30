#' Cohen's w
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @return value of Cohen's w
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' es_cohen_w(chi2Value, n)
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
#' The chi-square value could for example be obtained using the *ts_pearson_gof()* function.
#' 
#' For a classification using some rule-of-thumb run *th_cohen_w(w)*
#' 
#' **Alternative**
#' 
#' The *'rcompanion'* library also has a function for this: *cohenW()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @export
es_cohen_w <-function(chi2, n){
  #this function calculates Cohens w
  w = sqrt(chi2 / n)
  
  return (w)
}
