#' Johnston-Berry-Mielke E
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param minExp the minimum expected count
#' @param test c("chi", "g"). "chi" (default) for chi-square tests, or "g" for likelihood ratio (g) tests
#' @return JBM's E value
#' 
#' @examples 
#' chi2Value <- 3.105263
#' n <- 19
#' minExp <- n/4
#' es_jbm_e(chi2Value, n, minExp)
#' es_jbm_e(chi2Value, n, minExp, test="g")
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
#' A qualification rule-of-thumb could be obtained by converting this to Cohen's w, using *es_convert(from=="jbme", to=="cohenw", ex1 = minExp/n)*
#' The rule-of-thumb can then be obtained using *th_cohen_w(w)*
#' 
#' **Alternative**
#' 
#' I'm not aware of any other library with a similar function.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Johnston, J. E., Berry, K. J., & Mielke, P. W. (2006). Measures of effect size for chi-squared and likelihood-ratio goodness-of-fit tests. *Perceptual and Motor Skills, 103*(2), 412–414. https://doi.org/10.2466/pms.103.2.412-414
#'  
#' @export
es_jbm_e <-function(chi2, n, minExp, test="chi"){
  if (test=="chi"){
    E = chi2*minExp /(n*(n - minExp))
    
    }
  else{
    E = -1/log(minExp/n)*chi2/(2*n)
  }
  
  return (E)
}