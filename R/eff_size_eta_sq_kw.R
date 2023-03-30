#' Eta Squared for Kruskal-Wallis
#' 
#' @param H the H value of the Kruskal-Wallis test
#' @param n the sample size
#' @param k the number of categories
#' @return the effect size value
#' 
#' @details 
#' The formula used is (Tomczak & Tomczak, 2014, p. 24):
#' \deqn{\eta_{KW}^2 = \frac{H - k + 1}{n - k}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{H} the H value of the Kruskal-Wallis test
#' \item \eqn{n} the sample size
#' \item \eqn{k} the number of categories
#' }
#' 
#' Eta squared in general is much older, for example Pearson (1911, p. 254).
#' 
#' @references 
#' Pearson, K. (1911). On a correction to be made to the correlation ratio η. *Biometrika, 8*(1/2), 254. https://doi.org/10.2307/2331454
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
#' k = 3
#' es_eta_sq_kw(H, n)
#' 
#' @export
es_eta_sq_kw <- function(H, n, k){
  
  es = (H - k + 1)/(n - k)
  
  return(es)
}