#' Cohen's w
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @return value of Cohen's w
#' 
#' @description 
#' An effect size measure that could be used with a chi-square test. It has no upper limit, but can be compared to 
#' Cohen's rules-of-thumb. 
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
#' The chi-square value could for example be obtained using the \code{\link{ts_pearson_gof}} function.
#' 
#' For a classification using some rule-of-thumb run \code{\link{th_cohen_w}}
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' es_cohen_w(chi2Value, n)
#' 
#' @section Alternative:
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
#' @seealso 
#' 
#' Alternative effect sizes that might be of interest:
#' \itemize{
#' \item \code{\link{es_cramer_v_gof}} Cramér V
#' \item \code{\link{es_jbm_e}} Johnston-Berry-Mielke E
#' }
#' 
#' Tests with a nominal variable where this effect size might be used:
#' \itemize{
#' \item \code{\link{ts_pearson_gof}} Pearson chi-square test of goodness-of-fit
#' \item \code{\link{ts_g_gof}} G / Likelihood Ratio / Wilks test of goodness-of-fit
#' \item \code{\link{ts_freeman_tukey_gof}} Freeman-Tukey test of goodness-of-fit
#' \item \code{\link{ts_neyman_gof}} Neyman test of goodness-of-fit
#' \item \code{\link{ts_mod_log_likelihood_gof}} mod-log likelihood test of goodness-of-fit
#' \item \code{\link{ts_cressie_read_gof}} Cressie-Read / Power Divergence test of goodness-of-fit
#' }
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
es_cohen_w <-function(chi2, n){
  #this function calculates Cohens w
  w = sqrt(chi2 / n)
  
  return (w)
}
