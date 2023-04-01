#' Alexander-Govern Test
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
#' The formula used (Alexander & Govern, 1994, pp. 92-94):
#' \deqn{A = \sum_{j=1}^k z_j^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(A, df\right)}
#' 
#' With:
#' \deqn{z_j = c_j + \frac{c_j^3 + 3\times c_j}{b_j} - \frac{4\times c_j^7 + 33\times c_j^5 + 240\times c_j^3 + 855\times c_j}{10\times b_j^2 + 8\times b_j\times c_j^4 + 1000\times b_j}}
#' \deqn{c_j = \sqrt{a_j\times\ln\left(1 + \frac{t_j^2}{n_j - 1}\right)}}
#' \deqn{b_j = 48\times a_j^2}
#' \deqn{a_j = n_j - 1.5}
#' \deqn{t_j = \frac{\bar{x}_j - \bar{y}_w}{\sqrt{\frac{s_j^2}{n_j}}}}
#' \deqn{\bar{y}_w = \sum_{j=1}^k h_j\times \bar{x}_j}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{n} the total sample size
#' \item \eqn{k} the number of categories
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\chi^2\left(\dots, \dots\right)} the cumulative distribution function of the chi-square distribution.
#' }
#' 
#' @references 
#' Alexander, R. A., & Govern, D. M. (1994). A new and simpler approximation for ANOVA under variance heterogeneity. *Journal of Educational Statistics, 19*(2), 91–101. https://doi.org/10.2307/1165140
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
#' ts_alexander_govern(scores, groups)
#' 
#' @export
ts_alexander_govern <- function(scores, groups){
  
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("location", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("location", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("location", "var"))
  myRes <- merge(counts, means, by = 'location')
  myRes <- merge(myRes, vars, by = 'location')
  
  myRes$se <- sqrt(myRes$var/myRes$n)
  myRes$recSe <- 1/(myRes$se^2)
  myRes$w <- myRes$recSe/sum(myRes$recSe)
  
  muHat <- sum(myRes$w*myRes$mean)
  t <- (myRes$mean - muHat)/myRes$se
  v <- myRes$n - 1
  a <- v - 0.5
  b <- 48*a^2
  c <- sqrt(a*log(1 + t^2/v))
  z = c + (c^3 + 3*c)/b - (4*c^7 + 33*c^5 + 240*c^3 + 855*c)/(10*b^2 + 8*b*c^4 + 1000*b)
  A <- sum(z^2)
  k <- dim(myRes)[1]
  df <- k - 1
  pVal <- pchisq(A, df, lower.tail=FALSE)
  
  testResults <- data.frame(A, df, pVal)
  colnames(testResults)<-c("statistic", "df", "pValue")
  
  return (testResults)
  
}