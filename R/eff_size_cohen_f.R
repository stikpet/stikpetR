#' Cohen f
#' @description 
#' An effect size measure for regression analysis or an ANOVA test. It gives roughly the proportion of variance explained by the categorical variable.
#' 
#' The Cohen f is often used with ANOVA, while Cohen f-squared with regression.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' @param useRanks boolean, optional. Use of ranks or original scores. Default is FALSE
#' @returns the Cohen f value
#' 
#' @details 
#' The formula used (Cohen, 1988, p. 284):
#' \deqn{f = \sqrt{\frac{\eta^2}{1-\eta^2}}}
#' 
#' Where \eqn{\eta^2} is the value of eta-squared.
#' 
#' It can also be calculated using (Cohen, 1988, p. 371):
#' \deqn{f = \frac{\sigma_{\mu}}{\sigma}}
#' 
#' With:
#' \deqn{\sigma_{\mu} = \sqrt{\frac{SS_b}{n}}}
#' \deqn{\sigma = \sqrt{\frac{SS_w}{n}}}
#' 
#' Where \eqn{SS_i} is the sum of squared differences, see the Fisher one-way ANOVA for details on how to calculate these.
#' 
#' The \eqn{f^2} can be found in Cohen (1988, p. 410).
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
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_cohen_f <- function(nomField, scaleField, categories=NULL, useRanks=FALSE){
  eta2 = es_eta_sq(nomField, scaleField, categories, useRanks=useRanks)
  
  f2 = eta2 / (1 - eta2)
  
  return(f2**0.5)
}



