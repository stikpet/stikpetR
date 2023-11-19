#' Root Mean Square Standardized Effect Size (RMSSE)
#' @description 
#' An effect size measure for a one-way ANOVA.
#' 
#' Similar as Hedges g, but for a one-way ANOVA. According to Wikipedia "this essentially presents the omnibus difference of the entire model adjusted by the root mean square" (2023).    
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' @returns the rmsse value
#' 
#' @details 
#' The formula used (Steiger & Fouladi, 1997, pp. 244-245):
#' \deqn{RMSSE = \sqrt{\frac{\delta}{\left(k - 1\right)\times n}} = \sqrt{\frac{\sum_{i=1}^k \alpha_i^2}{\left(k - 1\right)\times \sigma^2}}}
#' With:
#' \deqn{\delta = n\times\sum_{i=1}^k \left(\frac{\alpha_i}{\sigma}\right)^2}
#' \deqn{\alpha_i = \mu_i - \mu \approx \bar{x}_i - \bar{x}}
#' \deqn{\sigma \approx \sqrt{MS_w}}
#' \deqn{MS_w = \frac{SS_w}{df_w}}
#' \deqn{df_w = n - k}
#' \deqn{SS_w = SS_t - SS_b}
#' \deqn{SS_t = \sum_{j=1}^k \sum_{i=1}^{n_j}\left(x_{i,j} - \bar{x}\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{j=1}^k n_j\times\bar{x}_j }{n} = \frac{\sum_{j=1}^k \sum_{i=1}^{n_j} x_{i,j}}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the number of scores in category j
#' \item \eqn{\bar{x}_j} the mean of the scores in category j
#' \item \eqn{\bar{x} the overall mean}
#' \item \eqn{SS_w} the within sum of squares (sum of squared deviation of the mean)
#' \item \eqn{df_w} the within degrees of freedom
#' }
#' 
#' Note that the original article refers to \eqn{\sigma^2} as the error variance of the noncentral F-distribution.
#' This can be approximated with \eqn{MS_w} (Smith & Down, 2014, p. 2).
#' 
#' Zhang and Algina (2011) create a robust version of the RMSSE for one-way fixed effects anova.
#'
#' @references 
#' Smith, B., & Dowd, M. (2014). One-way analysis of variance (ANOVA). Dalhousie University. https://www.mathstat.dal.ca/~stat2080/Fall14/Lecturenotes/anova1.pdf
#' 
#' Steiger, J. H., & Fouladi, R. T. (1997). *Noncentrality interval estimation and the evaluation of statistical models*. In L. L. Harlow, S. A. Mulaik, & J. H. Steiger, What if there were no significance tests? (pp. 221–257). Lawrence Erlbaum Associates.
#' 
#' Wikipedia. (2023). Effect size. In Wikipedia. https://en.wikipedia.org/w/index.php?title=Effect_size&oldid=1175948622
#' 
#' Zhang, G., & Algina, J. (2011). A robust root mean square standardized effect size in one-way fixed-effects ANOVA. *Journal of Modern Applied Statistical Methods, 10*(1), 77–96. https://doi.org/10.22237/jmasm/1304222880
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_rmsse <- function(nomField, scaleField, categories=NULL){
  
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
  dfw = n - k
  
  MSw = SSw/dfw
  
  m = sum(res$mean*res$n)/n
  a = (res$mean - m)
  
  rmsse = sqrt(sum(a**2)/((k - 1)*MSw))
  
  return(rmsse)
  
}
