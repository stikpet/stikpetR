#' Rosenthal Correlation Coefficient
#' 
#' @param zVal z-value of test
#' @param n total sample size
#' @return the effect size measure
#' 
#' @examples  
#' z = 1.143943
#' n = 20
#' r_rosenthal(z, n)
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
#' Rosenthal (1991) is the oldest reference I could find for this correlation coefficient.
#' However, Cohen (1988, p. 275) actually has a measure 'f' that has the same equation.
#' 
#' For a classification the same as for Pearson correlation use *th_pearson_r()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Rosenthal, R. (1991). *Meta-analytic procedures for social research* (Rev. ed). Sage Publications.
#'  
#' @export
r_rosenthal <- function(zVal, n){
  
  r = zVal/sqrt(n)
  
  return (r[1,1])
  
}