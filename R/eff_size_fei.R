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
#' @seealso 
#' \code{\link{es_convert}} to convert Fei to Cohen w, use from="fei", to="cohenw", and ex1=minExp/n
#' \code{\link{th_cohen_w}} rules-of-thumb for Cohen w
#' 
#' @references 
#' Ben-Shachar, M. S., Patil, I., Thériault, R., Wiernik, B. M., & Lüdecke, D. (2023). Phi, fei, fo, fum: Effect sizes for categorical data that use the chi-squared statistic. *Mathematics, 11*(1982), 1–10. doi:10.3390/math11091982
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_fei <-function(chi2, n, minExp){
  
  pe = minExp/n
  f = (chi2/(n*(1/pe - 1)))**0.5
  
  return (f)
}