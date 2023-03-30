#' Fligner-Policello Test
#' 
#' @param var1 A vector with the scores data
#' @param var2 A vector with the group data
#' @param ties boolean to indicate the use of a ties correction
#' @param corr boolean to indicate the use of a continuity correction
#' @return dataframe with the test statistic, p-value, and the test used
#' 
#' @details
#' The formula used is:
#' \deqn{z = \frac{N_Y - N_X}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
#' With:
#' \deqn{SS_X = \sum_{x\in X} \left(N_X - M_X\right)^2, SS_Y = \sum_{y\in Y} \left(N_Y - M_Y\right)^2}
#' \deqn{M_X = \frac{N_X}{n_x}, M_Y = \frac{N_Y}{n_y}}
#' \deqn{N_X = \sum_{x \in X} N\left(x\right), N_Y = \sum_{y \in Y} N\left(y\right)}
#' \deqn{N\left(y\right) = \sum_{x\in X} f\left(y, x\right)}
#' \deqn{N\left(x\right) = \sum_{y\in Y} f\left(x, y\right)}
#' \deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0 & \text{ if } a\leq b \end{cases}}
#' 
#' In case of a tie correction (Hollander et al., 2014, p. 146):
#' \deqn{z = \frac{\left|N_Y - N_X\right| - 0.5}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
#' \deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0.5 & \text{ if } a = b \\ 0 & \text{ if } a\leq b \end{cases}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{X} the scores in the first category
#' \item \eqn{Y} the scores in the second category
#' \item \eqn{n_i} the number of scores in the i category
#' }
#' 
#' The test is described by Fligner and Policello (1981), and can also be found in Kloke and McKean (2015, p. 68)
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Fligner, M. A., & Policello, G. E. (1981). Robust rank procedures for the Behrens-Fisher problem. *Journal of the American Statistical Association, 76*(373), 162–168. https://doi.org/10.1080/01621459.1981.10477623
#' 
#' Hollander, M., Wolfe, D. A., & Chicken, E. (2014). *Nonparametric statistical methods* (3rd ed.). John Wiley & Sons, Inc.
#' 
#' Kloke, J., & McKean, J. W. (2015). *Nonparametric statistical methods using R*. CRC Press, Taylor & Francis.
#' 
#' @examples 
#' scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
#' groups = c("A","A","A","B","B","B","B", NA, "C")
#' ts_fligner_policello(scores, groups)
#' ts_fligner_policello(scores, groups, ties=FALSE, cont=FALSE)
#' ts_fligner_policello(scores, groups, ties=FALSE, cont=TRUE)
#' ts_fligner_policello(scores, groups, ties=TRUE, cont=FALSE)
#' ts_fligner_policello(scores, groups, ties=TRUE, cont=TRUE)
#' 
#' @export
ts_fligner_policello <- function(dataVar, groupVar, ties= TRUE, cont=FALSE){
  testUsed = "Fligner-Policello test"
  
  #remove rows with missing values
  df = data.frame(dataVar, groupVar)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  data1 <- unlist(unname(split(df$score, df$group)[1]))
  data2 <- unlist(unname(split(df$score, df$group)[2]))
  
  n1 = c()
  for (i in data1) {
    s = 0
    for (j in data2) {
      if (i == j && ties) {
        s = s + 0.5
      } else if(i > j){
        s = s + 1
      }
    }
    n1 = c(n1, s)
  }
  
  n2 = c()
  for (i in data2) {
    s = 0
    for (j in data1) {
      if (i == j && ties) {
        s = s + 0.5
      } else if(i > j){
        s = s + 1
      }
    }
    n2 = c(n2, s)
  }
  
  s1 = sum(n1)
  s2 = sum(n2)
  
  m1 = s1/length(n1)
  m2 = s2/length(n2)
  
  ss1 = sum((n1 - m1)^2)
  ss2 = sum((n2 - m2)^2)
  
  se = sqrt(ss1 + ss2 + m1*m2)
  num = (s2 - s1)/2
  
  if (cont) {
    num = abs(num)-0.5
  }
  z = num/se
  pValue = 2*(1 - pnorm(abs(z)))
  
  
  if (cont && ties){
    testUsed = paste0(testUsed, ", with continuity and ties correction")
  }
  else if (cont){
    testUsed = paste0(testUsed, ", with continuity correction")
  }
  else if (ties){
    testUsed = paste0(testUsed, ", with ties correction")
  }
  
  statistic = z  
  testResults <- data.frame(statistic, pValue, testUsed)
  
  return (testResults)
  
}