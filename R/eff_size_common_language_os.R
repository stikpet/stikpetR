#' Common Language Effect Size (One-Sample)
#' @description 
#' The Common Language Effect Size is most often used for independent samples or paired samples, but some have adapted the concept for one-sample as well.
#' 
#' It is the probability of taking a random score and the probability it is higher than the selected value: 
#' $$P(X > \\mu_{H_0})$$
#' 
#' Some will also argue to count ties equally, which makes the definition:
#' $$P(X > \\mu_{H_0}) + \\frac{P(X = \\mu_{H_0})}{2}$$
#' 
#' This version is implemented in MatLab (see https://nl.mathworks.com/matlabcentral/fileexchange/113020-cles) based on a Python version from Tulimieri (2021) 
#' 
#' For scale data, an approximation using the standard normal distribution is also available using Cohen's d, alternatively a conversion via the rank-biserial coefficient can be done. These two are used in R’s *effectsize* library from Ben-Shachar et al. (2020).
#' 
#' @param scores list with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param mu optional hypothesized statistic, otherwise the midrange will be used
#' @param version optional {"brute", "brute-it", "rb", "normal"} method to use. see details
#' 
#' @returns 
#' CLES : float, the Common Language Effect Size
#' 
#' @details
#' For "brute" simply counts all scores above the test statistic and half of the ones that are equal (Tulimieri, 2021):
#' $$CL = P(X > \\mu_{H_0}) + \\frac{P(X = \\mu_{H_0})}{2}$$
#' With:
#' $$P\\left(x \\gt \\mu\\right) = \\frac{\\sum_{i=1}^n \\begin{cases} 1, & \\text{if } x_i \\gt \\mu \\\\ 0, & \\text{otherwise}\\end{cases}}{n}$$
#' $$P\\left(x = \\mu\\right) = \\frac{\\sum_{i=1}^n \\begin{cases} 1, & \\text{if } x_i = \\mu \\\\ 0, & \\text{otherwise}\\end{cases}}{n}$$
#' 
#' With "brute-it" the ties are ignored (it = ignore ties):
#' $$CL = P(X > \\mu_{H_0})$$
#' 
#' The "normal", uses Cohen's d and a normal approximation (Ben-Shachar et al., 2020):
#' $$CL = \\Phi\\left(\\frac{d'}{\\sqrt{2}}\\right)$$
#' Where \\(d'\\) is Cohen's d for one-sample, and \\(\\Phi\\left(\\dots\\right)\\) the cumulative density function of the normal distribution
#' This is like a one-sample version of the McGraw and Wong (1992, p. 361) version with the independent samples.
#' 
#' The "rb", uses the rank-biserial correlation coefficient (Ben-Shachar et al., 2020):
#' $$CL = \\frac{1+r_b}{2}$$
#' The CLE can be converted to a Rank Biserial (= Cliff delta) using the **es_convert()** function. This can then be converted to a Cohen d, and then the rules-of-thumb for Cohen d could be used (**th_cohen_d()**)
#' 
#' @seealso 
#' \code{\link{es_cohen_d_os}}, Cohen d' for one-sample
#' \code{\link{r_rank_biserial_os}}, Rank Biserial Correlation Coefficient (one-sample)
#' \code{\link{es_convert}}, to convert an CLE to a rank biserial use `fr="cle", to="rb"`. To convert the result to Cohen d, use `fr="rb", to="cohend"`.
#' \code{\link{th_cohen_d}}, rules of thumb for Cohen d
#' 
#' @references
#' Ben-Shachar, M., Lüdecke, D., & Makowski, D. (2020). effectsize: Estimation of Effect Size Indices and Standardized Parameters. *Journal of Open Source Software, 5*(56), 1–7. doi:10.21105/joss.02815
#'  
#' Grissom, R. J. (1994). Statistical analysis of ordinal categorical status after therapies. *Journal of Consulting and Clinical Psychology, 62*(2), 281–284. doi:10.1037/0022-006X.62.2.281
#' 
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. *Psychological Bulletin, 111*(2), 361–365. doi:10.1037/0033-2909.111.2.361
#'  
#' Tulimieri, D. (2021). CLES/CLES. https://github.com/tulimid1/CLES/tree/main
#'  
#' Wolfe, D. A., & Hogg, R. V. (1971). On constructing statistics and reporting data. *The American Statistician, 25*(4), 27–30. doi:10.1080/00031305.1971.10477278
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_common_language_os <- function(scores, levels=NULL, mu=NULL, version="brute"){
  if (!is.null(levels)){
    scores = factor(na.omit(scores), ordered = TRUE, levels = levels)
    scores = as.numeric(scores)
  }
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(scores) + max(scores)) / 2
  }
  n = length(scores)
  
  if (version=="brute-it"){        
    n_gt = sum(sapply(scores, function(i) i > mu))
    cles = n_gt / n}
  else if (version=="brute"){
    n_gt = sum(sapply(scores, function(i) i > mu))
    n_eq = sum(sapply(scores, function(i) i == mu))
    cles = n_gt / n + 1/2*(n_eq/n)}
  else if (version=="rb"){
    rb = r_rank_biserial_os(scores, mu=mu)[1,2]
    cles = (1 + rb)/2}
  else if (version=="normal"){
    d_os = es_cohen_d_os(scores, mu=mu)
    cles = pnorm(d_os/(2**0.5))}
  return (cles)
}