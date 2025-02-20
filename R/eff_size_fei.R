#' Fei
#' 
#' @description 
#' An effect size measure that could be used with a chi-square test or g-test. 
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param minExp the minimum expected count
#' 
#' @return the value of Fei
#' 
#' @details 
#' The formula used (Ben-Shachar et al., 2023, p. 6):
#' \deqn{Fei = \\sqrt{\\frac{\\chi_{GoF}^2}{n\\times\\left(\\frac{1}{\\min\\left(p_E\\right)}-1\\right)}}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{\\chi_{GoF}^2}, the chi-square value of the goodness-of-fit chi-square test
#' \item \eqn{n}, the sample size
#' \item \eqn{p_E}, the expected proportions
#' }
#' 
#' *Classification*
#' A qualification rule-of-thumb could be obtained by converting this to Cohen's w (use **es_convert(Fei, fr="fei", to="cohenw", ex1=minExp/n)**)
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
#' After this you might want to use some rule-of-thumb for the interpretation by converting it to Cohen w:
#' \code{\link{es_convert}} to convert Fei to Cohen w (using fr="fei", to="cohenw", ex1=minExp/n).
#' \code{\link{th_cohen_w}} for various rules-of-thumb for Cohen w.
#' 
#' Alternative effect sizes that use a chi-square value:
#' \code{\link{es_cohen_w}}, for Cohen w.
#' \code{\link{es_cramer_v_gof}},  for Cramer's V for Goodness-of-Fit.
#' \code{\link{es_jbm_e}}, for Johnston-Berry-Mielke E.
#' 
#' @references 
#' Ben-Shachar, M. S., Patil, I., Thériault, R., Wiernik, B. M., & Lüdecke, D. (2023). Phi, fei, fo, fum: Effect sizes for categorical data that use the chi-squared statistic. *Mathematics, 11*(1982), 1–10. doi:10.3390/math11091982
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' chi2 = 23.5
#' n = 53
#' minExp = 14
#' es_fei(chi2=chi2, n=n, minExp=minExp)
#' 
#' @export
es_fei <-function(chi2, n, minExp){
  
  pe = minExp/n
  f = (chi2/(n*(1/pe - 1)))**0.5
  
  return (f)
}