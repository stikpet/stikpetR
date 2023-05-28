#' Consensus
#' 
#' @param ordData ordinal data
#' @return cns the consensus score
#' 
#' @description 
#' The Consensus is a measure of agreement or dispersion for ordinal data. If there 
#' is no agreement the value is 0, and with full agreement 1. 
#' 
#' @details 
#' The formula used (Tastle et al., 2005, p. 98):
#' \deqn{\text{Cns}\left(X\right) = 1 + \sum_{i=1}^k p_i \log_2\left(1 - \frac{\left|X_i - \mu_X\right|}{d_X}\right)}
#' With:
#' \deqn{\mu_X = \frac{\sum_{i=1}^k X_i\times F_i}{n}}
#' \deqn{d_X = \max\left(X_i\right) - \min\left(X_i\right)}
#' \deqn{p_i = \frac{F_i}{n}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{X_i} the rank of category i
#' \item \eqn{F_i} the frequency (count) of the i-th category (after they have been sorted)
#' \item \eqn{n} the sample size
#' \item \eqn{k} the number of categories.
#' }
#' 
#' @examples 
#' ordData <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' me_consensus(ordData)
#' 
#' labels = c("fully disagree", "disagree", "neutral", "agree", "fully agree")
#' factData = as.factor(ordData)
#' levels(factData) = labels
#' me_consensus(factData)
#' 
#' @section Alternatives:
#' 
#' The *agrmt* library has a function *consensus(table(ordData))*
#' 
#' @references 
#' Tastle, W. J., & Wierman, M. J. (2007). Consensus and dissention: A measure of ordinal dispersion. *International Journal of Approximate Reasoning, 45*(3), 531–545. https://doi.org/10.1016/j.ijar.2006.06.024
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
me_consensus <- function(ordData){
  freq = table(ordData)
  n = sum(freq)
  p = freq/n
  r = seq(length(freq))
  m = sum(freq*r/n)
  d = max(r) - min(r)
  
  cns = 1 + sum(p*log(1 - abs(r - m)/d, base=2))
  
  return(cns)  
}