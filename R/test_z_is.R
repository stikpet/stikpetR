#' Independent Samples Z Test
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @param dmu difference according to null hypothesis (default is 0)
#' @param sigma1 population standard deviation of the first group, if NULL sample results will be used
#' @param sigma2 population standard deviation of the second group, if NULL sample results will be used
#' @returns 
#' A dataframe with:
#' \item{statistic}{the test statistic}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details
#' 
#' The formula used is:
#' \deqn{z = \frac{\bar{x}_1 - \bar{x}_2}{SE}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}}}
#' \deqn{\sigma_i^2 \approx s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' ts_z_is(scores, groups)
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @export
ts_z_is <- function(scores, groups, dmu=0, sigma1=NULL, sigma2=NULL){
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
  
  if (is.null(sigma1)) {
    var1 = var(X1)
  }
  else{
    var1 = sigma1**2
  }
  
  if (is.null(sigma2)) {
    var2 = var(X2)
  }
  else{
    var2 = sigma2**2
  }
  
  sse = var1/n1 + var2/n2
  se = sqrt(sse)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  z = (m1 - m2 - dmu)/se
  
  pValue = 2*(1-pnorm(abs(z)))
  
  statistic = z
  
  testResults <- data.frame(statistic, pValue)
  
  return(testResults)
  
}
