#' Rosenthal Correlation Coefficient
#' 
#' This function will calculate Rosenthal Correlation Coefficient. A simple correlation coefficient that divides a z-score by the square root of the sample size.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/ycip6vF-7_U) and the effect size is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Correlations/RosenthalCorrelationCoefficient.html)
#' 
#' @param zVal z-value of test
#' @param n total sample size
#' 
#' @return r the effect size measure
#' 
#' @details 
#' The formula used (Rosenthal, 1991, p. 19):
#' \deqn{r = \frac{z}{\sqrt{n}}}
#' 
#' *Symbols used:*
#' \itemize{
#' 	\item \eqn{n} the sample size
#' 	\item \eqn{z} the calculated z-statistic value
#' }
#' 
#' Rosenthal (1991) is the oldest reference I could find for this correlation coefficient. However, Cohen (1988, p. 275) actually has a measure 'f' that has the same equation.
#' 
#' For a classification the same as for Pearson correlation use *th_pearson_r()*
#' 
#' @section Before, After and Alternatives:
#' Before this measure you might want to perform the test:
#' \code{\link{ts_sign_os}}, for One-Sample Sign Test.
#' \code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
#' \code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
#' 
#' After this you might want to use the rules-of-thumb for a Pearson Correlation or Cohen f:
#' \code{\link{th_pearson_r}}, for rules of thumb for a Pearson correlation coefficient
#' \code{\link{th_cohen_f}}, for rules of thumb for a Cohen f
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Rosenthal, R. (1991). *Meta-analytic procedures for social research* (Rev. ed). Sage Publications.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' z = 1.143943
#' n = 20
#' r_rosenthal(z, n)
#'  
#' @export
r_rosenthal <- function(zVal, n){
  
  r = zVal/sqrt(n)
  
  return (r)
  
}