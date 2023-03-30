#' Cohen d (for one-way ANOVA)
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns the effect size value
#' 
#' @details 
#' The formula used is (Cohen, 1988, p. 276):
#' \deqn{d = \frac{\bar{x}_{max} - \bar{x}_{min}}{\sigma}}
#' With:
#' \deqn{\sigma = \sqrt{\frac{SS_w}{n}}}
#' \deqn{SS_w = \sum_{j=1}^k \sum_{i=1}^{n_j} \left(x_{i,j}-\bar{x}_j\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n_j} the number of scores in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{\bar{x}_j} the mean of the scores in category j
#' \item \eqn{SS_w} the within sum of squares (sum of squared deviation of the mean)
#' }
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
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
#' es_cohen_d_ow(scores, groups)
#' 
#' @export
es_cohen_d_ow <- function(scores, groups){
  dfr = na.omit(data.frame(scores, groups))
  
  counts <- setNames(aggregate(dfr$scores~dfr$groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(dfr$scores~dfr$groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(dfr$scores~dfr$groups, FUN=var), c("category", "var"))
  res <- merge(counts, means, by = 'category')
  res <- merge(res, vars, by = 'category')
  
  SSw = sum(res$var*(res$n-1))
  
  k <- dim(res)[1]
  n = sum(res$n)
  
  
  d = (max(res$mean) - min(res$mean))/(sqrt(SSw/n))
  return(d)
  
}