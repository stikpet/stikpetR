#' Common Language (CL/CLES) (Independent Samples)
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @return The Common Language Effect Size value
#' 
#' @details
#' 
#' The formula used is (McGraw & Wong, 1992, p. 361):
#' \deqn{CL = \Phi\left(z\right)}
#' 
#' With:
#' \deqn{z = \frac{\left|\bar{x}_1 - \bar{x}_2\right|}{\sqrt{s_1^2 + s_2^2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' \item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. *Psychological Bulletin, 111*(2), 361–365. https://doi.org/10.1037/0033-2909.111.2.361
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' es_common_language_is(scores, groups)
#' 
#' @export
es_common_language_is <- function(scores, groups){

  #make sure data is numeric
  scores = as.numeric(scores)
  
  #remove rows with missing values
  df = data.frame(scores, groups)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  X1 = df$score[df$group == df$group[1]]
  X2 = df$score[df$group != df$group[1]]
  
  var1 = var(X1)
  var2 = var(X2)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  z = abs(m1 - m2)/sqrt(var1 + var2)
  
  cl = pnorm(z)
  
  return(cl)
}