#' Point Biserial Correlation Coefficient
#' @description
#' This can be seen as coding a binary variable with the groups into 0 and 1, and then calculates a Pearson correlation coefficient between the those values and the scores.
#' 
#' This gives the same result as the formula used and as input the Student t-test statistic and corresponding degrees of freedom.
#' 
#' @param t the test statistic value
#' @param df the degrees of freedom
#' 
#' @return Point Biserial Correlation Coefficient
#' 
#' @details
#' The formula used is (Friedman, 1968, p. 245):
#' \deqn{r_{pb} = \sqrt{\frac{t^2}{t^2 + df}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{t} the test statistic of the independent samples Student t-test
#' \item \eqn{df} the degrees of freedom of the independent samples Student t-test
#' }
#' 
#' @seealso 
#' \code{\link{ts_student_t_is}}, Student t-test
#' 
#' @references 
#' Friedman, H. (1968). Magnitude of experimental effect and a table for its rapid estimation. *Psychological Bulletin, 70*(4), 245–251. https://doi.org/10.1037/h0026258
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' r_point_biserial(0.9984, 1967)
#' 
#' @export
r_point_biserial <- function(t, df){
  
  r = sqrt(t**2 / (t**2 + df))
  
  return(r)
  
}