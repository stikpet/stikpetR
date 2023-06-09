#' Johnston-Berry-Mielke E
#' 
#' @description 
#' An effect size measure that could be used with a chi-square test or g-test. 
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param minExp the minimum expected count
#' @param test optional to indicate if a chi-square tests, or a g (likelihood ratio) test was used. Either `"chi"` (default), or `"g"`.
#' 
#' @return JBM's E value
#' 
#' @details 
#' Two versions of this effect size. The formula for a chi-square test is:
#' \deqn{E_{\chi^2}=\frac{q}{1-q}\times \left(\sum_{i=1}^{k}\frac{p_{i}^{2}}{q_{i}}-1\right) = \frac{\chi_{GoF}^2\times E_{min}}{n\times\left(n - E_{min}\right)}}
#' For a Likelihood Ratio (G) test:
#' \deqn{E_{L}=-\frac{1}{\textup{ln}(q)}\times\sum_{i=1}^{k}\left(p_{i}\times\textup{ln}\left(\frac{p_{i}}{q_{i}}\right)\right) = -\frac{1}{\textup{ln}\left(q\right)\times\frac{\chi_L^2}{2\times n}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{q} the minimum of all \eqn{q_i}
#' \item \eqn{q_i} the expected proportion in category i
#' \item \eqn{p_i} the observed proportion in category i
#' \item \eqn{n} the total sample size
#' \item \eqn{E_min} the minimum expected count
#' \item \eqn{\chi_{GoF}^2} the chi-square test statistic of a Pearson chi-square test of goodness-of-fit
#' \item \eqn{\chi_L^2} the chi-square test statistic of a likelihood ratio test of goodness-of-fit
#' }
#' 
#' Both formulas are from Johnston et al. (2006, p. 413)
#' 
#' A qualification rule-of-thumb could be obtained by converting this to Cohen's w
#' 
#' @seealso 
#' \code{\link{es_convert}} to convert JBM-E to Cohen w, use from="jbme", to="cohenw", and ex1=minExp/n
#' \code{\link{th_cohen_w}} rules-of-thumb for Cohen w
#' 
#' @references 
#' Johnston, J. E., Berry, K. J., & Mielke, P. W. (2006). Measures of effect size for chi-squared and likelihood-ratio goodness-of-fit tests. *Perceptual and Motor Skills, 103*(2), 412–414. https://doi.org/10.2466/pms.103.2.412-414
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' minExp <- n/4
#' es_jbm_e(chi2Value, n, minExp)
#' es_jbm_e(chi2Value, n, minExp, test="g")
#' 
#' @export
es_jbm_e <-function(chi2, n, minExp, test=c("chi", "g")){
  
  if (length(test)>1) {test="chi"}
  
  if (test=="chi"){
    E = chi2*minExp /(n*(n - minExp))
    
    }
  else if (test=="g"){
    E = -1/log(minExp/n)*chi2/(2*n)
  }
  
  return (E)
}