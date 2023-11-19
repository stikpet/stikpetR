#' Epsilon Squared
#' @description 
#' An effect size measure to indicate the the strength of the categories on the ordinal/scale field. A 0 would indicate no influence, and 1 a perfect relationship.
#' 
#' This is an attempt to make eta-squared unbiased (applying a population correction ratio) (Kelley, 1935, p. 557). Although a popular belief is that omega-squared is preferred over epsilon-squared (Keselman, 1975), a later study actually showed that epsilon-squared might be preferred (Okada, 2013).
#' 
#' Tomczak and Tomczak (2014) recommend this as one option to be used with a Kruskal-Wallis test, however I think they labelled epsilon-squared as eta-squared and the other way around.
#'
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param useRanks boolean, optional. Use ranks or use the scores as given in ordfield. Default is FALSE.
#' 
#' @returns
#' epsSq, float. The epsilon squared value
#' 
#' @details
#' The formula used (Kelley, 1935, p. 557):
#' \deqn{\epsilon^2 = \frac{n\times\eta^2 - k + \left(1 - \eta^2\right)}{n - k}}
#'
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\eta^2} eta squared
#' \item \eqn{k} the number of categories
#' \item \eqn{n} the sample size
#' }
#' 
#' There are quite some variations on this formula. 
#' 
#' For example Cureton (1966, p. 605):
#' \deqn{\epsilon^2 = 1 - \frac{n - 1}{n - k}\times\left(1 - \eta^2\right)}
#' Caroll and Nordholm (1975, p. 547):
#' \deqn{\epsilon^2 = \frac{F - 1}{F + \frac{n - k}{k - 1}}}
#' Albers and Lakens (2018, p. 194)
#' \deqn{\epsilon^2 =  \frac{F - 1}{F + \frac{df_w}{df_b}}}
#' Albers and Lakens (2018, p. 188)
#' \deqn{\epsilon^2 = \frac{SS_b - df_b\times MS_w}{SS_t}}
#'
#' **Conversions**
#'
#' To convert \eqn{\epsilon^2} to \eqn{\eta^2} use *es_conver(epsilonsq, from="epsilonsq", to="etasq", ex1=n, ex2=k)*
#'
#' To convert \eqn{\epsilon^2} to \eqn{\omega^2} use *es_convert(epsilonsq, from="etasq", to="omegasq", ex1=MS_w, ex2=SS_t)*
#'
#' **Alternatives**
#'
#' *library(effectsize)*
#'
#' anova_stats(aov(scores~groups))
#'
#' epsilon_squared(aov(scores~groups))
#'
#' @references
#' Albers, C., & Lakens, D. (2018). When power analyses based on pilot data are biased: Inaccurate effect size estimators and follow-up bias. *Journal of Experimental Social Psychology, 74*, 187–195. doi:10.1016/j.jesp.2017.09.004
#'
#' Carroll, R. M., & Nordholm, L. A. (1975). Sampling characteristics of Kelley’s \eqn{\epsilon} and Hays’ \eqn{\omega}. *Educational and Psychological Measurement, 35*(3), 541–554. doi:10.1177/001316447503500304
#'
#' Cureton, E. E. (1966). On correlation coefficients. *Psychometrika, 31*(4), 605–607. doi:10.1007/BF02289528
#'
#' Kelley, T. L. (1935). An unbiased correlation ratio measure. *Proceedings of the National Academy of Sciences, 21*(9), 554–559. doi:10.1073/pnas.21.9.554
#'
#' Keselman, H. J. (1975). A Monte Carlo investigation of three estimates of treatment magnitude: Epsilon squared, eta squared, and omega squared. *Canadian Psychological Review / Psychologie Canadienne, 16*(1), 44–48. doi:10.1037/h0081789
#'
#' Okada, K. (2013). Is omega squared less biased? A comparison of three major effect size indices in one-way anova. *Behaviormetrika, 40*(2), 129–147. doi:10.2333/bhmk.40.129
#' 
#' Pearson, K. (1911). On a correction to be made to the correlation ratio \eqn{\eta}. *Biometrika, 8*(1/2), 254. doi:10.2307/2331454
#' 
#' Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates revisited. An overview of some recommended measures of effect size. *Trends in Sport Sciences, 1*(21), 19–25.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_epsilon_sq <- function(catField, ordField, categories=NULL, levels=NULL, useRanks=FALSE){
  
  ct = tab_cross(ordField, catField, order1=levels, order2=categories, totals="include")
  nLvl = nrow(ct)-1 
  
  k = ncol(ct)-1
  n = ct[nLvl+1, k+1]
  
  e2 =  es_eta_sq(catField, ordField, categories, levels, useRanks)
  epsSq = (n * e2 - k + (1 - e2)) / (n - k)
  
  return(epsSq)
}
