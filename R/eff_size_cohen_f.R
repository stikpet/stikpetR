#' Cohen f
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns the effect size value
#' 
#' @details 
#' The formula used is (Cohen, 1988, p. 371):
#' \deqn{f = \sqrt{\frac{SS_b}{SS_w}}}
#' With:
#' \deqn{SS_b = \sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}
#' \deqn{SS_w = SS_t - SS_b}
#' \deqn{SS_t = \sum_{j=1}^k \sum_{i=1}^{n_j}\left(x_{i,j} - \bar{x}\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{j=1}^k n_j\times\bar{x}_j }{n} = \frac{\sum_{j=1}^k \sum_{i=1}^{n_j} x_{i,j}}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the number of scores in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{\bar{x}_j} the mean of the scores in category j
#' \item \eqn{SS_i} the sum of squares of i (sum of squared deviation of the mean)
#' \item \eqn{b} is between = factor = treatment = model
#' \item \eqn{w} is within = error (the variability within the groups)
#' }
#' 
#' Cohen shows as formula \eqn{f = \frac{\sigma_{\mu}}{\sigma}} with:
#' \deqn{\sigma_{\mu} = \sqrt{\frac{SS_b}{n}}}
#' \deqn{\sigma = \sqrt{\frac{SS_w}{n}}}
#' 
#' 
#' Often Cohen f is calculated using eta-squared. This can also be done 
#' using the conversion function: *es_convert(etasq, from="etasq", to="cohenf")*
#' 
#' **Conversions**
#' 
#' Cohen f can be converted to eta-squared using: *es_convert(f, from="cohenf", to="etasq")*
#' 
#' **Alternatives**
#' 
#' *library(effectsize)*
#' 
#' anova_stats(aov(scores~groups))
#' 
#' cohens_f(aov(scores~groups))
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
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
#' es_cohen_f(scores, groups)
#' 
#' @export
es_cohen_f <- function(scores, groups){
  dfr = na.omit(data.frame(scores, groups))
  
  #overall mean and count
  xBar = mean(dfr$scores)
  n = nrow(dfr)
  
  #means and counts per category
  xBars = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=mean)$x
  ns = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=length)$x
  
  #number of categories
  k = length(ns)
  
  SSb = sum(ns*(xBars - xBar)**2)
  SSt = var(dfr$scores)*(n - 1)
  SSw = SSt - SSb
  
  es = sqrt(SSb/SSw)
  
  return(es)
}
