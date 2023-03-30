#' Mehrotra Test
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns 
#' A dataframe with:
#' \item{statistic}{the F-statistic from the test}
#' \item{df1}{the first degrees of freedom}
#' \item{df2}{the second degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Mehrotra, 1997, p. 11141):
#' \deqn{F_{M} = \frac{\sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}{\sum_{j=1}^k\left(1 - \frac{n_j}{n}\right)\times s_j^2}}
#' \deqn{df_1 = \frac{\left(\sum_{j=1}^k s_j^2 - \frac{n_j\times s_j^2}{n}\right)^2}{\sum_{j=1}^k s_j^4 + \left(\frac{\sum_{j=1}^k n_j\times s_j^2}{n}\right)^2 - 2\times\frac{\sum_{j=1}^k n_j\times s_j^4}{n}}}
#' \deqn{df_2 =\frac{\left(\sum_{j=1}^k\left(1 - \frac{n_j}{n}\right)\times s_j^2\right)^2}{\sum_{j=1}^k\frac{\left(1 - \frac{n_j}{n}\right)\times s_j^4}{n_j - 1}}} 
#' \deqn{sig. = 1 - F\left(F_{BF}, df_1, df_2\right)}
#' With:
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^k n_j\times\bar{x}_j}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{F\left(\dots,\dots,\dots\right)} the cumulative distribution function of the F distribution.
#' }
#' 
#' The same as the Brown-Forsythe test for means, except for \eqn{df_1}.
#' 
#' @references 
#' Mehrotra, D. V. (1997). Improving the Brown-Forsythe solution to the generalized Behrens-Fisher problem. *Communications in Statistics - Simulation and Computation, 26*(3), 1139–1145. https://doi.org/10.1080/03610919708813431
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' scores = c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, NA, 90, 25, 40, 70, 65, NA, 70, 98, 40, 65, 60, 35, NA, 50, 40, 75, NA, 65, 70, NA, 20, 80, 35, NA, 68, 70, 60, 70, NA, 80, 98, 10, 40, 63, 75, 80, 40, 90, 100, 33, 36, 65, 78, 50)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' ts_mehrota(scores, groups)
#' 
#' @export
ts_mehrota <- function(scores, groups){
  
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  n <- sum(myRes$n)
  xbar <- sum(myRes$n*myRes$mean)/n
  Fstat <- sum(myRes$n*(myRes$mean - xbar)^2) / sum((1 - myRes$n/n)*myRes$var) 
  k <- dim(myRes)[1]
  df1 <- (sum(myRes$var) - (sum(myRes$n*myRes$var)/n))^2 / (sum(myRes$var^2) + (sum(myRes$n*myRes$var)/n)^2 - 2*sum(myRes$n*myRes$var^2)/n)
  df2 <- sum((1 - myRes$n/n)*myRes$var)^2 / sum((1 - myRes$n/n)^2*myRes$var^2/(myRes$n - 1))
  pVal <- pf(Fstat, df1, df2, lower.tail = FALSE)
  
  testResults <- data.frame(Fstat, df1, df2, pVal)
  colnames(testResults)<-c("statistic", "df1", "df2", "pValue")
  
  return (testResults)
  
}