#' Cohen's w
#' 
#' @description 
#' An effect size measure that could be used with a chi-square test. It has no upper limit, but can be compared to Cohen's rules-of-thumb. 
#' 
#' This function is shown in this [YouTube video](https://youtu.be/0l7SmTlmuYs) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CohenW.html)
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
#' @section Before, After and Alternatives:
#' Before this you will need a chi-square value. From either:
#' \code{\link{ts_freeman_tukey_gof}}, for Freeman-Tukey Test of Goodness-of-Fit. 
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_mod_log_likelihood_gof}}, for Mod-Log Likelihood Test of Goodness-of-Fit. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_pearson_gof}}, for Pearson Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test. 
#' \code{\link{ph_pairwise_gof}} for Pairwise Goodness-of-Fit Tests.
#' \code{\link{ph_residual_gof_gof}} for Residuals Using Goodness-of-Fit Tests
#' 
#' After this you might want to use some rule-of-thumb for the interpretation:
#' \code{\link{th_cohen_w}} for various rules-of-thumb for Cohen w.
#' 
#' Alternative effect sizes that use a chi-square value:
#' \code{\link{es_cramer_v_gof}},  for Cramer's V for Goodness-of-Fit.
#' \code{\link{es_fei}}, for Fei.
#' \code{\link{es_jbm_e}}, for Johnston-Berry-Mielke E.
#' 
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' chi2Value <- 3.106
#' n <- 19
#' es_cohen_w(chi2Value, n)
#' 
#' @export
es_cohen_w <-function(chi2, n){
  #this function calculates Cohens w
  w = sqrt(chi2 / n)
  
  return (w)
}
