#' Wilcox Test
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns 
#' A dataframe with:
#' \item{statistic}{the chi-square-statistic from the test}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Wilcox, 1988, pp. 110-111)
#' \deqn{H = \frac{\sum_{j=1}^k \left(W_j - \bar{W}\right)^2}{\hat{\theta}}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(H, df\right)}
#' With:
#' \deqn{W_j = b_j\times x_{n_j,j} + \frac{1 - b_j}{n_j}\times\sum_{i=1}^{n_j-1} x_{i,j}}
#' \deqn{\bar{W} = \frac{\sum_{j=1}^k W_j}{k}}
#' \deqn{b_j = \frac{1 + \sqrt{\frac{\left(n_j - 1\right)\times\left(n_j\times\hat{\theta}\right)}{s_j^2}}}{n_j}}
#' \deqn{\hat{\theta} = \text{max}\left\{\frac{s_1^2}{n_1}, \frac{s_2^2}{n_2}, \dots, \frac{s_k^2}{n_k}\right\}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\chi^2\left(\dots\right)} the cumulative density function of the chi-square distribution
#' }
#' 
#' The original article has an error in the formula for \eqn{b_j}. There are missing brackets. 
#' Using the population version in the article of \eqn{c_j} the formula used here was adapted.
#' 
#' @references 
#' Wilcox, R. R. (1988). A new alternative to the ANOVA F and new results on James’s second-order method. *British Journal of Mathematical and Statistical Psychology, 41*(1), 109–117. https://doi.org/10.1111/j.2044-8317.1988.tb00890.x
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
#' ts_wilcox(scores, groups)
#' 
#' @export 
ts_wilcox <- function(scores, groups){
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("group", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("group", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("group", "var"))
  myRes <- merge(counts, means, by = 'group')
  myRes <- merge(myRes, vars, by = 'group')
  
  d = max(myRes$var/myRes$n)
  bj = (1 + sqrt((myRes$n - 1)*(myRes$n*d - myRes$var)/myRes$var))/myRes$n
  
  k <- dim(myRes)[1]
  Wj = rep(0, k)
  for (j in 1:k) {
    Wj[j] = (1 - bj[j])/myRes$n[j] *sum(datFrame$scores[datFrame$groups==myRes$group[j]][1:(myRes$n[j]-1)])
    Wj[j] = Wj[j] + bj[j]*datFrame$scores[datFrame$groups==myRes$group[j]][myRes$n[j]]
  }
  
  Wm = sum(Wj)/k
  
  H = sum((Wj - Wm)**2)/d
  df = k - 1
  pValue = 1 - pchisq(H, df)
  
  testResults <- data.frame(H, df, pValue)
  colnames(testResults)<-c("statistic", "df", "pValue")
  
  return (testResults)
  
  
}