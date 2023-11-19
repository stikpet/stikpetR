#' Eta Squared
#' @description 
#' An effect size measure to indicate the the strength of the categories on the ordinal/scale field. A 0 would indicate no influence, and 1 a perfect relationship.
#' 
#' It is “the proportion of the variation in Y that is associated with membership of the different groups defined by X “ (Richardson, 2011, p. 136).
#' 
#' An alternative Epsilon Squared is an attempt to make eta-squared unbiased (applying a population correction ratio) (Kelley, 1935, p. 557). Although a popular belief is that omega-squared is preferred over epsilon-squared (Keselman, 1975), a later study actually showed that epsilon-squared might be preferred (Okada, 2013).
#' 
#' Tomczak and Tomczak (2014) recommend this this as one option to be used with a Kruskal-Wallis test, however I think they labelled epsilon-squared as eta-squared and the other way around.
#'
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param useRanks boolean, optional. Use ranks or use the scores as given in ordfield. Default is FALSE.
#' 
#' @returns
#' etaSq, float. The eta squared value
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
#' If ranks are used, the eta-squared can also be determined using (Tomczak & Tomczak, 2014, p. 24):
#' \deqn{\eta^2 = \frac{H}{n - 1}}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#'  \item \eqn{n}, the total sample size
#'  \item \eqn{k}, the number of categories
#'  \item \eqn{SS_b}, the between sum of squares (sum of squared deviation of the mean)
#'  \item \eqn{SS_t}, the total sum of squares (sum of squared deviation of the mean)
#'  \item \eqn{F}, the F-statistic
#'  \item \eqn{H}, H-statistic from Kruskal-Wallis H-test
#'  \item \eqn{df_i}, the degrees of freedom of i
#'  \item \eqn{x_{i,j}}, the i-th score in category j
#'  \item \eqn{n_j}, the number of scores in category j
#'  \item \eqn{\bar{x}_j}, the mean of the scores in category j
#'  \item \eqn{b}, is between = factor = treatment = model
#'  \item \eqn{w}, is within = error (the variability within the groups)
#'  }
#'  
#'  Eta-squared can be converted to Cohen f, using *es_convert(etasq, from="etasq", to="cohenf")*
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
#' Kelley, T. L. (1935). An unbiased correlation ratio measure. *Proceedings of the National Academy of Sciences, 21*(9), 554–559. doi:10.1073/pnas.21.9.554
#' 
#' Keselman, H. J. (1975). A Monte Carlo investigation of three estimates of treatment magnitude: Epsilon squared, eta squared, and omega squared. *Canadian Psychological Review / Psychologie Canadienne, 16*(1), 44–48. doi:10.1037/h0081789
#' 
#' Okada, K. (2013). Is omega squared less biased? A comparison of three major effect size indices in one-way anova. *Behaviormetrika, 40*(2), 129–147. doi:10.2333/bhmk.40.129
#' 
#' Pearson, K. (1911). On a correction to be made to the correlation ratio \eqn{\eta}. *Biometrika, 8*(1/2), 254. doi:10.2307/2331454
#' 
#' Richardson, J. T. E. (2011). Eta squared and partial eta squared as measures of effect size in educational research. *Educational Research Review, 6*(2), 135–147. doi:10.1016/j.edurev.2010.12.001
#' 
#' Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates revisited. An overview of some recommended measures of effect size. *Trends in Sport Sciences, 1*(21), 19–25.
#'
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_eta_sq <- function(catField, ordField, categories=NULL, levels=NULL, useRanks=FALSE){
  
  #replace levels if provided
  if (!is.null(levels)){
    myFieldOrd = factor(ordField, ordered = TRUE, levels = levels)
    ordField = as.numeric(myFieldOrd)
  }
  
  datFrame = na.omit(data.frame(catField, ordField))
  
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$catField %in% categories),]}
  
  if (useRanks){
    datFrame["ordField"] = rank(datFrame$ordField)}
  
  counts <- setNames(aggregate(datFrame$ordField~datFrame$catField, FUN=length), c("group", "n"))
  means <- setNames(aggregate(datFrame$ordField~datFrame$catField, FUN=mean), c("group", "mean"))
  vars <- setNames(aggregate(datFrame$ordField~datFrame$catField, FUN=var), c("group", "var"))
  myRes <- merge(counts, means, by = 'group')
  myRes <- merge(myRes, vars, by = 'group')
  
  n = sum(myRes$n)
  xBar = sum(myRes$n*myRes$mean)/n
  SSb = sum(myRes$n*(myRes$mean - xBar)**2)
  SSt = var(datFrame$ordField)*(n - 1)
  SSw = SSt - SSb
  
  es = SSb/SSt
  
  return(es)
}
