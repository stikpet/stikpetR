#' Eta Squared
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns the effect size value
#' 
#' @description 
#' Eta squared is “the proportion of the variation in Y that is associated with membership of 
#' the different groups defined by X “ (Richardson, 2011, p. 136).
#' 
#' @details
#' The formula used is (Pearson, 1911, p. 254):
#' \deqn{\eta^2 = \frac{SS_b}{SS_t}}
#' With:
#' \deqn{SS_t = \sum_{j=1}^k\sum_{i=1}^{n_j}\left(x_{i,j}-\bar{x}\right)^2}
#' \deqn{SS_b = \sum_{j=1}^k\left(\bar{x}_j-\bar{x}\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = {\sum_{j=1}^k\sum_{i=1}^{n_j} x_{i,j}}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' }
#' 
#' There are variations on the formula that will give the same result, for example:
#' \deqn{\eta^2 = \frac{F\times\left(k - 1\right)}{F\times\left(k - 1\right) + n - k}}
#' or
#' \deqn{\eta^2 = \frac{F\times df_b}{F\times df_b + df_w}}
#' 
#' Eta-squared can be converted to Cohen f, using *es_convert(etasq, from="etasq", to="cohenf")*
#' 
#' Eta-squared can be converted to Epsilon square, using *es_convert(etasq, from="etasq", to="epsilonsq", ex1=n, ex2=k)*
#' 
#' **Alternatives**
#' 
#' *library(lsr)* 
#' 
#' etaSquared(aov(scores~groups))
#' 
#' *library(effectsize)*
#' 
#' anova_stats(aov(scores~groups))
#' 
#' eta_squared(aov(scores~groups))
#' 
#' @references 
#' Pearson, K. (1911). On a correction to be made to the correlation ratio η. *Biometrika, 8*(1/2), 254. https://doi.org/10.2307/2331454
#' 
#' Richardson, J. T. E. (2011). Eta squared and partial eta squared as measures of effect size in educational research. *Educational Research Review, 6*(2), 135–147. https://doi.org/10.1016/j.edurev.2010.12.001
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
#' es_eta_sq(scores, groups)
#' 
#' @export 
es_eta_sq <- function(scores, groups){
  
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("group", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("group", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("group", "var"))
  myRes <- merge(counts, means, by = 'group')
  myRes <- merge(myRes, vars, by = 'group')
  
  n = sum(myRes$n)
  xBar = sum(myRes$n*myRes$mean)/n
  SSb = sum(myRes$n*(myRes$mean - xBar)**2)
  SSt = var(datFrame$scores)*(n - 1)
  SSw = SSt - SSb
  
  es = SSb/SSt
  
  return(es)
}