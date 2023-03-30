#' Epsilon Squared for Kruskal-Wallis
#' 
#' @param H the H value of the Kruskal-Wallis test
#' @param n the sample size
#' @return the effect size value
#' 
#' @details 
#' The formula used is (Tomczak & Tomczak, 2014, p. 24):
#' \deqn{\epsilon_{KW}^2 = \frac{H}{n - 1}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{H} the H value of the Kruskal-Wallis test
#' \item \eqn{n} the sample size
#' }
#' 
#' Epsilon squared is mentioned and in general described by Kelley (1935).
#' 
#' @references 
#' Kelley, T. L. (1935). An unbiased correlation ratio measure. *Proceedings of the National Academy of Sciences, 21*(9), 554–559. https://doi.org/10.1073/pnas.21.9.554
#' 
#' Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates revisited. An overview of some recommended measures of effect size. *Trends in Sport Sciences, 1*(21), 19–25.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' H = 21.32807
#' n = 54
#' es_epsilon_sq_kw(H, n)
#' 
#' @export
es_epsilon_sq_kw <- function(H, n){
  
  es = H/(n - 1)
  
  return(es)
}