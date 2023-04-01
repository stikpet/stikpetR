#' Vargha and Delaney A
#' 
#' @param dataVar A vector with the scores data
#' @param groupVar A vector with the group data
#' @return Vargha and Delaney A value
#' 
#' @details
#' 
#' The formula used is (Vargha & Delaney, 2000, p. 107):
#' \deqn{A = \frac{1}{n_j}\times\left(\frac{R_i}{n_i} - \frac{n_i + 1}{2}\right)}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_i} the number of scores in category i
#' \item \eqn{R_i} the sum of the ranks in category i
#' }
#' 
#' This effect size is an adaptation of the Common Language effect size, 
#' adapted for ordinal data.
#' 
#' It could also be calculated from the Mann-Whitney U value:
#' \deqn{A = \frac{U}{n_1\times n_2}}
#' 
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101–132. https://doi.org/10.3102/10769986025002101
#' 
#' @examples 
#' scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
#' groups = c("A","A","A","B","B","B","B", NA, "C")
#' es_vargha_delaney_a(scores, groups)
#' 
#' @export
es_vargha_delaney_a <- function(dataVar, groupVar){
  
  #make sure data is numeric
  scores = as.numeric(dataVar)
  
  #remove rows with missing values
  df = data.frame(scores, groupVar)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  n = length(df$group)
  n1 = sum(df$group==df$group[1])
  n2 = n - n1

  rankScores = rank(df$score)
  
  R1 = sum(rankScores[df$group==df$group[1]])
  
  A = 1/n2 * (R1/n1 - (n1 + 1)/2)
  
  return(A)
  
}

var1 = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
var2 = c("A","A","A","B","B","B","B", NA, "C")

es_vargha_delaney_a(var1, var2)
