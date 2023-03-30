#' Independent Samples Trimmed/Yuen Mean Test
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized trimmed mean, otherwise the midrange will be used
#' @param trim optional proportion to trim from each side (so in total twice this will be trimmed)
#' @return dataframe test statistic, degrees of freedom, p-value (sig.) and name of test used
#'  
#' @details 
#' The formula used is (Yuen & Dixon, 1973, p. 394):
#' \deqn{t = \frac{\bar{x}_{t,1} - \bar{x}_{t,2}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{SSD_{w,1} + SSD_{w,2}}{m_1 + m_2 - 2}\times\left(\frac{1}{m_1} + \frac{1}{m_2}\right)}}
#' \deqn{df = m_1 + m_2 - 2}
#' \deqn{\bar{x}_{t,i} = \frac{\sum_{j=g_i+1}^{n_i - g_i}y_{i,j}}{}}
#' \deqn{g_i = \lfloor n_i\times p_t\rfloor}
#' \deqn{m_i = n_ - 2\times g_i}
#' \deqn{SSD_{w,i} = g_i\times\left(y_{i,g_i+1} - \bar{x}_{w,i}\right)^2 + g_i\times\left(y_{i,n_i-g_i} - \bar{x}_{w,i}\right)^2 + \sum_{j=g+1}^{n_i - g_i} \left(y_{i,j} - \bar{x}_{w,i}\right)^2}
#' \deqn{\bar{x}_{w,i} = \frac{\bar{x}_{t,i}\times m_i + g_i\times\left(y_{i, g_i+1} + y_{i, n_i-g_i}\right)}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{ti}} the trimmed mean of the scores in category i
#' \item \eqn{x_{wi}} The Winsorized mean of the scores in category i
#' \item \eqn{SSD_{wi}} the sum of squared deviations from the Winsorized mean of category i
#' \item \eqn{m_i} the number of scores in the trimmed data set from category i
#' \item \eqn{y_{i,j}} the j-th score after the scores in category i, after they are sorted from low to high
#' \item \eqn{p_t} the proportion of trimming on each side, we can define
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
#' Yuen, K. K., & Dixon, W. J. (1973). The approximate behaviour and performance of the two-sample trimmed t. *Biometrika, 60*(2), 369–374. https://doi.org/10.2307/2334550
#'  
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' ts_trimmed_means_t_is(scores, groups)
#' 
#' @export
ts_trimmed_means_t_is <- function(scores, groups, trimProp = 0.1){
  df = data.frame(scores, groups)
  df = na.omit(df)
  
  scores1 = df$scores[df$groups==df$groups[1]]
  n1 = length(scores1)
  n1t = n1 - 2*round(n1*trimProp)
  m1t = mean(scores1, trim=trimProp)
  #note that R base mean(..., trim) will round down
  
  scores2 = df$scores[df$groups!=df$groups[1]]
  n2 = length(scores2)
  n2t = n2 - 2*round(n2*trimProp)
  m2t = mean(scores2, trim=trimProp)
  
  #Winsorize the data
  sc1Sort = sort(scores1)
  minReplace1 = sc1Sort[round(n1*trimProp)+1]
  maxReplace1 = sc1Sort[n1 - round(n1*trimProp)]
  scores1 = replace(scores1, scores1<minReplace1, minReplace1)
  scores1 = replace(scores1, scores1>maxReplace1, maxReplace1)
  
  sc2Sort = sort(scores2)
  minReplace2 = sc2Sort[round(n2*trimProp)+1]
  maxReplace2 = sc2Sort[n2 - round(n2*trimProp)]
  scores2 = replace(scores2, scores2<minReplace2, minReplace2)
  scores2 = replace(scores2, scores2>maxReplace2, maxReplace2)
  var1 = var(scores1)
  var2 = var(scores2)
  
  s2 = ((n1 - 1)*var1 + (n2 - 1)*var2)/((n1t - 1) + (n2t - 1))
  
  tValue = (m1t - m2t) / (sqrt(s2 * (1/n1t + 1/n2t)))
  
  df = n1t + n2t - 2
  
  pValue = 2*(1 - pt(abs(tValue), df))
  
  statistic=tValue
  testResults <- data.frame(statistic, df, pValue)
  
  return(testResults)
  
}