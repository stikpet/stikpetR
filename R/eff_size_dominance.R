#' Dominance and a Vargha-Delaney A like effect size measure
#' 
#' @param data dataframe with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param mu optional parameter to set the hypothesized median. If not used the midrange is used
#' @param out c("dominance","vda"). optional to either show the dominance score or a VDA like measure: `"dominance"` (default), `"vda"`
#' 
#' @return dataframe with the hypothesized median (mu) and the effect size value
#' 
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
#' @section Alternatives:
#' The *rcompanion* library has a similar function: *oneSampleDominance()*
#' 
#' @references 
#' Mangiafico, S. S. (2016). Summary and analysis of extension program evaluation in R (1.20.01). Rutger Cooperative Extension.
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples  
#' #Example 1: Text dataframe
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' es_dominance(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' es_dominance(ex2)
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
  title = "dominance"
  
  if (out=="vda") {
    VDA = (dominance + 1)/2
    es = VDA
    title = "VDA-like"
  }
  
  
  results <- data.frame(mu, es)
  colnames(results) = c("mu", title)
  
  return(results)
}