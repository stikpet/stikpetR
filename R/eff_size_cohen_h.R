#' Cohen h
#' 
#' @param p1 the first proportion
#' @param p2 the second proportion
#' @return Cohen h
#' 
#' @details 
#' Formula used (Cohen, 1988, p. 181):
#' \deqn{h=\phi_{1}-\phi_{2}}
#' With:
#' \deqn{\phi_{i}=2\times\textup{arcsin}\sqrt{p_{i}}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{p_i} the proportion of cases in category i
#' }
#' 
#' For classification rule-of-thumb use: *th_cohen_h()*
#' 
#' @examples 
#' es_cohen_h(0.2, 0.4)
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
#' @export
es_cohen_h <- function(p1, p2){
  
  phi1 = 2 * asin(sqrt(p1))
  phic = 2 * asin(sqrt(p2))
  
  h = phi1 - phic
  
  return (h)

}