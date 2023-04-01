#' (Glass) Rank Biserial Correlation / Cliff Delta
#' 
#' @param dataVar A vector with the scores data
#' @param groupVar A vector with the group data
#' @return Rank Biserial Correlation value
#' 
#' @details
#' 
#' The formula used is (Glass, 1966, p. 626):
#' \deqn{r_b = \frac{2\times\left(\bar{R}_1 - \bar{R}_2\right)}{n}}
#' With:
#' \deqn{\bar{R}_i=\frac{R_i}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\bar{R}_i} the average of ranks in category i
#' \item \eqn{R_i} the sum of ranks in category i
#' \item \eqn{n} the total sample size
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' Glass (1966) showed that the formula was the same as that of the 
#' rank biserial from Cureton (1956). Cliff's delta (Cliff, 1993, p. 495) 
#' is actually also the same.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. *Psychological Bulletin, 114*(3), 494–509. https://doi.org/10.1037/0033-2909.114.3.494
#' 
#' Cureton, E. E. (1956). Rank-biserial correlation. *Psychometrika, 21*(3), 287–290. https://doi.org/10.1007/BF02289138
#' 
#' Glass, G. V. (1966). Note on rank biserial correlation. *Educational and Psychological Measurement, 26*(3), 623–631. https://doi.org/10.1177/001316446602600307
#' 
#' @examples 
#' scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
#' groups = c("A","A","A","B","B","B","B", NA, "C")
#' r_rank_biserial_is(scores, groups)
#' 
#' @export
r_rank_biserial_is <- function(dataVar, groupVar){
  
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
  R1a = R1 / n1
  R2 = n*(n+1)/2 - R1
  R2a = R2 / n2
  
  rb = 2*(R1a - R2a)/n
  
  return(rb)
  
}
