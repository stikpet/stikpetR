#' Dominance and a Vargha-Delaney A like effect size measure
#' 
#' @param data dataframe with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param mu optional parameter to set the hypothesized median. If not used the midrange is used
#' @param out c("dominance","vda"). optional to either show the dominance score (default), or a VDA like measure
#' @return dataframe with the hypothesized median (mu) and the effect size value
#' 
#' @examples  
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' es_dominance(data)
#' es_dominance(data, mu=2)
#' 
#' @details 
#' The formula used is (Mangiafico, 2016, p. 223-224):
#' \deqn{D = p_{pos} - p_{neg}}
#' Where:
#' \deqn{p_i = \frac{n_i}{n}}
#' *Symbols used:*
#' \itemize{
#' \item \eqn{p_{pos}} the proportion of cases above the hypothesized median
#' \item \eqn{p_{neg}} the proportion of cases below the hypothesized median
#' \item \eqn{n_{pos}} the number of cases above the hypothesized median
#' \item \eqn{n_{neg}} the number of cases below the hypothesized median
#' \item \eqn{n} the total number of cases
#' }
#' 
#' The dominance score will range from -1 to 1.
#' 
#' A Vargha-Delaney A (VDA) style effect size is calculated with (Mangiafico, 2016, p. 223-224):
#' \deqn{VDA_{like} = \frac{D + 1}{2}}
#' 
#' This will range from 0 to 1, with 0.5 being the same as a dominance score of 0.
#' 
#' **Alternative**
#' The *rcompanion* library has a similar function: *oneSampleDominance()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Mangiafico, S. S. (2016). Summary and analysis of extension program evaluation in R (1.20.01). Rutger Cooperative Extension.
#'  
#' @export
es_dominance <- function(data, levels=NULL, mu=NULL, out="dominance"){
  
  if (!is.null(levels)){
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    data = as.numeric(myFieldOrd)
  }
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }
  
  #total sample size
  n = length(data)
  
  pPlus = sum(data > mu)/n
  pMin = sum(data < mu)/n
  
  dominance = pPlus - pMin
  es = dominance
  
  if (out=="vda") {
    VDA = (dominance + 1)/2
    es = VDA
  }
  
  return(es)
}