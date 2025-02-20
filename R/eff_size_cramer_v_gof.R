#' Cramer's V for Goodness-of-Fit
#' 
#' @description 
#' Cramér's V is one possible effect size when using a chi-square test. This measure is actually designed for the chi-square test for independence but can be adjusted for the goodness-of-fit test (Kelley & Preacher, 2012, p. 145; Mangiafico, 2016, p. 474). 
#' 
#' It gives an estimate of how well the data then fits the expected values, where 0 would indicate that they are exactly equal. If you use the equal distributed expected values the maximum value would be 1, otherwise it could actually also exceed 1.
#' 
#' As for a classification Cramér's V can be converted to Cohen w, for which Cohen provides rules of thumb.
#' 
#' A Bergsma correction is also possible.
#' 
#' A general explanation can also be found in this [YouTube video](https://youtu.be/FZcnk4EYpek). The function is shown in this [YouTube video]() and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CramerV.html)
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param k the number of categories
#' @param bergsma optional boolean to indicate the use of the Bergsma correction: `FALSE`(the default), or `TRUE`
#' 
#' @return Cramer's V value
#' 
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
#' Cramér described V (1946, p. 282) for use with a test of independence. Others (e.g. K. Kelley & Preacher, 2012, p. 145; Mangiafico, 2016a, p. 474) added that this can also be use for goodness-of-fit tests.
#' 
#' For the Bergsma (2013, pp. 324-325) correction the same thing applies
#' 
#' Cramér's V can be converted to Cohen's w using *es_convert(from="cramervgof", to = "cohenw", ex1 = df)*
#' 
#' Rules-of-thumb for the interpretation can then be used, using *th_cohen_w(w)*
#' 
#' @section Alternatives:
#' 
#' The *lsr* library has a similar function: *cramersV()*
#' 
#' The *DescTools* library has a similar function: *CramerV()*
#' 
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
#' \code{\link{th_cramer_v}} for various rules-of-thumb for Cramér V.
#' 
#' or convert to Cohen w:
#' \code{\link{es_convert}} to convert Cramér's V to Cohen w (using fr="cramervgof", to="cohenw", ex1=k).
#' \code{\link{th_cohen_w}} for various rules-of-thumb for Cohen w.
#' 
#' Alternative effect sizes that use a chi-square value:
#' \code{\link{es_cohen_w}}, for Cohen w.
#' \code{\link{es_fei}}, for Fei.
#' \code{\link{es_jbm_e}}, for Johnston-Berry-Mielke E.
#' 
#' 
#' @references 
#' Bergsma, W. (2013). A bias-correction for Cramér’s and Tschuprow’s. *Journal of the Korean Statistical Society, 42*(3), 323–328. https://doi.org/10.1016/j.jkss.2012.10.002
#' 
#' Cramér, H. (1946). *Mathematical methods of statistics*. Princeton University Press.
#' 
#' Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological Methods, 17*(2), 137–152. https://doi.org/10.1037/a0028086
#' 
#' Mangiafico, S. S. (2016). Summary and analysis of extension program evaluation in R (1.13.5). Rutger Cooperative Extension.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' chi2Value <- 3.106
#' n <- 19
#' k <- 3
#' es_cramer_v_gof(chi2Value, n, k)
#' es_cramer_v_gof(chi2Value, n, k, bergsma=TRUE)
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