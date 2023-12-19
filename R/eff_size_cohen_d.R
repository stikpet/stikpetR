#' Cohen d (for one-way ANOVA)
#' @description 
#' An effect size measure for a one-way ANOVA. It simply compares the largest possible difference between two categories means and divides this over the total variance.
#' 
#' Note that most often Cohen d is reported with pairwise tests, but that is actually Cohen d_z. That version is available using es_cohen_d_ps().
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' @returns the Cohen d value
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
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_cohen_d <- function(nomField, scaleField, categories=NULL){
  
  dfr = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$nomField %in% categories),]}
  
  counts <- setNames(aggregate(dfr$scaleField~dfr$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(dfr$scaleField~dfr$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(dfr$scaleField~dfr$nomField, FUN=var), c("category", "var"))
  res <- merge(counts, means, by = 'category')
  res <- merge(res, vars, by = 'category')
  
  SSw = sum(res$var*(res$n-1))
  
  k <- dim(res)[1]
  n = sum(res$n)
  
  
  d = (max(res$mean) - min(res$mean))/(sqrt(SSw/n))
  return(d)
  
}