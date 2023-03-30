#' Welch One-Way ANOVA
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
#' The formula used is (Welch, 1951, pp. 334-335):
#' \deqn{F_w = \frac{\frac{1}{k - 1} \times \sum_{j=1}^k w_j\times\left(\bar{x}_j - \bar{y}_w\right)^2}{1 + 2\times\lambda\times\frac{k-2}{k^2 - 1}}}
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 = \frac{k^2 - 1}{3\times\lambda}}
#' \deqn{sig. = 1 - F\left(F_W, df_1, df_2\right)}
#' 
#' With:
#' \deqn{\lambda = \sum_{j=1}^k \frac{\left(1 - h_j\right)^2}{n_j - 1}}
#' \deqn{\bar{y}_w = \frac{\sum_{j=1}^k w_j\times\bar{x}_j}{\sum_{j=1}^k w_j} = \sum_{j=1}^k h_j\times\bar{x}_j}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{x_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{w_j} the weight for category j
#' \item \eqn{h_j} the adjusted weight for category j
#' \item \eqn{df_i} the i-th degrees of freedom
#' }
#' 
#' The formula can also be written as:
#' \deqn{F_W = \frac{\chi_{Cochran}^2}{k - 1 + 2\times\lambda\times\frac{k - 2}{k + 1}}}
#' Where \eqn{\chi_{Cochran}^2} is the test statistic of the Cochran one-way test
#' 
#' Cavus and Yazici (2020) make a difference between the Welch and the Welch-Aspin ANOVA. The only difference in the article is that with the Welch 2×(k-2) is used, 
#' while in the Welch-Aspin version 2×k-2. I think this is a mistake in their formula, since the article they refer to from Aspin is about two means.
#' 
#' Johansen F test (Johansen, 1980) will give the same results
#' @references 
#' Johansen, S. (1980). The Welch-James approximation to the distribution of the residual sum of squares in a weighted linear regression. *Biometrika, 67*(1), 85–92. https://doi.org/10.1093/biomet/67.1.85
#' 
#' Welch, B. L. (1947). The generalization of `Student’s’ problem when several different population variances are involved. *Biometrika, 34*(1/2), 28–35. https://doi.org/10.2307/2332510
#' 
#' Welch, B. L. (1951). On the comparison of several mean values: An alternative approach. *Biometrika, 38*(3/4), 330–336. https://doi.org/10.2307/2332579
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
#' ts_welch_owa(scores, groups)
#' 
#' @export
ts_welch_owa <- function(scores, groups){
  
  dfr = na.omit(data.frame(scores, groups))
  
  #overall mean and count
  xBar = mean(dfr$scores)
  n = nrow(dfr)
  
  #means and counts and variances per category
  xBars = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=mean)$x
  ns = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=length)$x
  vars = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=var)$x
  
  #number of categories
  k = length(ns)
  
  wj = ns/vars
  w = sum(wj)
  
  hj = wj/w
  yw = sum(hj*xBars)
  
  lambda = sum((1 - hj)**2/(ns - 1))
  
  Fvalue = (sum(wj*(xBars - yw)**2)/(k - 1)) / (1 + 2*(k - 2)/(k**2 - 1)*lambda)
  
  df1 = k - 1
  df2 = (k**2 - 1)/(3*lambda)
  
  pValue = 1 - pf(Fvalue, df1, df2)
  
  statistic = Fvalue
  results = data.frame(statistic, df1, df2, pValue)
  
  return(results)
  
}