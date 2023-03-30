#' Welch t Test (Independent Samples)
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @return dataframe with the test statistic, degrees of freedom and p-value
#' 
#' @details
#' 
#' The formula used is:
#' \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{SE}}
#' \deqn{df = \frac{SE^4}{\frac{\left(s_1^2\right)^2}{n_1^2\times\left(n_1 - 1\right)} + \frac{\left(s_2^2\right)^2}{n_2^2\times\left(n_2 - 1\right)}}}
#' \deqn{sig. = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
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
#' Welch, B. L. (1947). The generalization of `Student’s’ problem when several different population variances are involved. *Biometrika, 34*(1/2), 28–35. https://doi.org/10.2307/2332510
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' ts_welch_t_is(scores, groups)
#' 
#' @export
ts_welch_t_is <- function(scores, groups){
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
  
  sse = var1/n1 + var2/n2
  se = sqrt(sse)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  t = (m1 - m2)/se
  
  df = sse**2/(var1**2/(n1**2*(n1 - 1)) + var2**2/(n2**2*(n2 - 1)))
  
  pValue = 2*(1-pt(abs(t), df))
  
  statistic = t
  
  testResults <- data.frame(statistic, df, pValue)
  
  return(testResults)
  
}