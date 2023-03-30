#' Point Biserial Correlation Coefficient
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @return Point Biserial Correlation Coefficient
#' 
#' @details
#' 
#' The formula used is (Friedman, 1968, p. 245):
#' \deqn{r_{pb} = \sqrt{\frac{t^2}{t^2 + df}}}
#' 
#' With:
#' \deqn{df = n_1 + n_2 - 2}
#' \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{SE}}
#' \deqn{SE = s_p\times\sqrt{\frac{1}{n_1} + \frac{1}{n_2}}}
#' \deqn{s_p = \sqrt{\frac{\left(n_1 - 1\right)\times s_1^2 + \left(n_2 - 1\right)\times s_2^2}{df}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{t} the test statistic of the independent samples Student t-test
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' This is codes the binary variable with the groups into 0 and 1, and then
#' calculates a Pearson correlation coefficient between the those values and
#' the scores.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Friedman, H. (1968). Magnitude of experimental effect and a table for its rapid estimation. *Psychological Bulletin, 70*(4), 245–251. https://doi.org/10.1037/h0026258
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' r_point_biserial(scores, groups)
#' 
#' @export
r_point_biserial <- function(scores, groups){
  #make sure data is numeric
  scores = as.numeric(scores)
  
  #remove rows with missing values
  df = data.frame(scores, groups)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  X1 = df$score[df$group == df$group[1]]
  X2 = df$score[df$group != df$group[1]]
  
  n1 = length(X1)
  n2 = length(X2)
  
  var1 = var(X1)
  var2 = var(X2)
  
  sp = sqrt(((n1 - 1)*var1 + (n2 - 1)*var2)/(n1 + n2 - 2))
  
  se = sp*sqrt(1/n1 + 1/n2)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  t = (m1 - m2)/se
  
  df = n1 + n2 - 2
  
  rpb = sqrt(t**2/(t**2 + df))
  
  if (t < 0) {
    rpb = -rpb
  }
  
  return(rpb)
  
}