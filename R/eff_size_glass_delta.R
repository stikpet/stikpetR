#' Glass Delta
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @param conrol Optional to indicate which group is the 'control' group, otherwise the 1st found will be used.
#' @return Hedges g value
#' 
#' @details
#' The formula used is (Glass, 1976, p. 7):
#' \deqn{\delta = \frac{\bar{x}_1 - \bar{x}_2}{s_2}}
#' 
#' With:
#' \deqn{s_2 = \sqrt{\frac{\sum_{i=1}^{n_2} \left(x_{2,i} - \bar{x}_2\right)^2}{n_2 - 1}}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' Glass actually uses a ‘control group’ and \eqn{s_2} is then the standard deviation of the control group. 
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. *Educational Researcher, 5*(10), 3–8. https://doi.org/10.3102/0013189X005010003
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' es_glass_delta(scores, groups)
#' es_glass_delta(scores, groups, control="international")
#' 
#' @export
es_glass_delta <- function(scores, groups, control=NULL){
  #make sure data is numeric
  scores = as.numeric(scores)
  
  #remove rows with missing values
  df = data.frame(scores, groups)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  if (is.null(control)){
    control = df$group[1]
  }
  
  X1 = df$score[df$group != control]
  X2 = df$score[df$group == control]

  var2 = var(X2)
  
  mean1 = mean(X1)
  mean2 = mean(X2)
  
  n = length(df$group)
  n1 = sum(df$group==df$group[1])
  n2 = n - n1
  
  sd2 = sqrt(var2)
  
  gd = (mean1 - mean2)/sd2
  
  return(gd)
  
}