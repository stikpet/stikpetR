#' Box F-Test
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns 
#' A dataframe with:
#' \item{statistic}{the chi-square-statistic from the test}
#' \item{df1}{the first degrees of freedom}
#' \item{df2}{the second degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Box, 1954, p. 299):
#' \deqn{F_{B} = \frac{F_F}{c}}
#' \deqn{df_1 = \frac{\left(\sum_{j=1}^k\left(n - n_j\right)\times s_j^2\right)^2}{\left(\sum_{j=1} n_j\times s_j^2\right)^2 + n\times\sum_{j=1}^k\left(n - 2\times n_j\right)\times s_j^4}}
#' \deqn{df_2 =\frac{\left(\sum_{j=1}^k\left(n_j-1\right)\times s_j^2\right)^2}{\sum_{j=1}^k\left(n_j-1\right)\times s_j^4}} 
#' \deqn{sig. = 1 - F\left(F_B, df_1, df_2\right)}
#' With:
#' \deqn{c = \frac{n - k}{n\times\left(k - 1\right)}\times \frac{\sum_{j=1}^k\left(n-n_j\right)\times s_j^2}{\sum_{j=1}^k\left(n_j-1\right)\times s_j^2}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{F-F} the F statistic of the classic/Fisher one-way ANOVA. See *ts_fisher_owa()* for details.
#' \item \eqn{n} the total sample size
#' \item \eqn{k} the number of categories
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{F\left(\dots,\dots,\dots\right)} the cumulative distribution function of the F distribution.
#' }
#' 
#' This also appears to give the same results for the test statistic, df_1 as the Brown-Forsythe test for means
#' but a different df_2
#' 
#' The *doex* and *onewaytests* libraries used to have a different method for calculating df_1 but after
#' personal communication with the creators of those packages they mentioned to fix it in an update.
#' 
#' @references 
#' Box, G. E. P. (1954). Some theorems on quadratic forms applied in the study of analysis of variance problems, I: Effect of inequality of variance in the one-way classification. *The Annals of Mathematical Statistics, 25*(2), 290–302. https://doi.org/10.1214/aoms/1177728786
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
#' ts_box(scores, groups)
#' 
#' @export
ts_box <- function(scores, groups){
  
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  k <- dim(myRes)[1]
  n = sum(myRes$n)
  xbar = mean(datFrame$scores)
  
  b = (n - k)/(n*(k - 1)) * sum((n - myRes$n)*myRes$var)/sum((myRes$n - 1)*myRes$var)
  Fstat = ((n - k)/(k - 1))*((sum(myRes$n*(myRes$mean-xbar)^2))/(sum((myRes$n-1)*myRes$var)))
  Fadj = Fstat / b
  
  df1 = sum((n - myRes$n)*myRes$var)^2/(sum(myRes$n*myRes$var)^2 + n*sum((n - 2*myRes$n)*myRes$var^2))
  
  df2 = sum((myRes$n - 1)*myRes$var)^2/sum((myRes$n-1)*myRes$var^2)
  
  pVal <- pf(Fadj, df1, df2, lower.tail = FALSE)
  
  testResults <- matrix(c(Fadj, df1, df2, pVal), ncol=4)
  colnames(testResults)<-c("statistic", "df1", "df2", "pValue")
  testResults <- as.table(testResults)
  
  return (testResults)
  
}