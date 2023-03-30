#' Omega Squared
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @param version which version to use for the calculation (see details)
#' @returns the effect size value
#' 
#' @details 
#' The formula used when version="hays1" is (Hays, 1973, p. 485):
#' \deqn{\hat{\omega}^2 = \frac{SS_b - df_b\times MS_w}{SS_t + MS_w}}
#' With:
#' \deqn{MS_w = \frac{SS_w}{df_w}}
#' \deqn{df_b = k - 1}
#' \deqn{df_w = n - k}
#' \deqn{SS_b = \sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}
#' \deqn{SS_w = SS_t - SS_b}
#' \deqn{SS_t = \sum_{j=1}^k \sum_{i=1}^{n_j}\left(x_{i,j} - \bar{x}\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{j=1}^k n_j\times\bar{x}_j }{n} = \frac{\sum_{j=1}^k \sum_{i=1}^{n_j} x_{i,j}}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' If version="hays2" the formula used is (Hays, 1973, p. 486):
#' \deqn{\hat{\omega} = \frac{\frac{\left(df_w - 2\right)\times F_F}{df_w}-1}{\frac{df_w + 1}{df_b} + \frac{\left(df_w - 2\right)\times F_F}{df_w}}}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the number of scores in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{\bar{x}_j} the mean of the scores in category j
#' \item \eqn{MS_i} the mean square of i
#' \item \eqn{SS_i} the sum of squares of i (sum of squared deviation of the mean)
#' \item \eqn{df_i} the degrees of freedom of i
#' \item \eqn{b} is between = factor = treatment = model
#' \item \eqn{w} is within = error (the variability within the groups)
#' \item \eqn{F_F} the test statistic of the Fisher/Classic one-way ANOVA
#' }
#' 
#' The formula appears in many different formats. Hays (1973, p. 486) shows:
#' \deqn{\hat{\omega}^2 = \frac{F_F - 1}{\frac{df_w + 1}{df_b}+ F_F}}
#' Which can also be found in Albers and Lakens (2018, p. 194).
#' 
#' Kirk (1996, p. 751) shows:
#' \deqn{\hat{\omega}^2 = \frac{df_b\times\left(F_F - 1\right)}{df_b\times\left(F_F - 1\right) + n}}
#' 
#' Instead of using the definitions for \eqn{df_b} and \eqn{df_w}, Caroll and Nordholm (1975, p. 547) show the formula as:
#' \deqn{\hat{\omega}^2 = \frac{F_F - 1}{\frac{n - k + 1}{k - 1}+ F_F}}
#' They also show Hays original formula (hays1) on p. 188.
#' 
#' Olejnik and Algina (2003, p. 435) use:
#' \deqn{\hat{\omega}^2 = \frac{SS_b - df_b\times MS_w}{SS_b + \left(n - df_b\right)\times MS_w}}
#' 
#' **Conversion**
#' 
#' To convert \eqn{\omega^2} to \eqn{\epsilon^2} use  *es_convert(omegasq, from="omegasq", to="epsilonsq", ex1=MS_w, ex2=SS_t)*
#' 
#' **Alternatives**
#' 
#' *library(effectsize)*
#' 
#' anova_stats(aov(scores~groups))
#' 
#' omega_squared(aov(scores~groups))
#' 
#' @references 
#' Albers, C., & Lakens, D. (2018). When power analyses based on pilot data are biased: Inaccurate effect size estimators and follow-up bias. *Journal of Experimental Social Psychology, 74*, 187–195. https://doi.org/10.1016/j.jesp.2017.09.004
#' 
#' Carroll, R. M., & Nordholm, L. A. (1975). Sampling characteristics of Kelley’s ε and Hays’ ω. *Educational and Psychological Measurement, 35*(3), 541–554. https://doi.org/10.1177/001316447503500304
#' 
#' Hays, W. L. (1973). *Statistics for the social sciences* (2nd ed.). Holt, Rinehart and Winston.
#' 
#' Kirk, R. E. (1996). Practical significance: A concept whose time has come. *Educational and Psychological Measurement, 56*(5), 746–759. https://doi.org/10.1177/0013164496056005002
#' 
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared statistics: Measures of effect size for some common research designs. *Psychological Methods, 8*(4), 434–447. https://doi.org/10.1037/1082-989X.8.4.434
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
#' es_omega_sq(scores, groups)
#' es_omega_sq(scores, groups, version="hays2")
#' 
#' 
#' @export 
es_omega_sq <- function(scores, groups, version=c("hays1", "hays2")){
  
  if (length(version)>1) {
    version="hays1"
  }
  
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
  
  k = length(myRes$n)
  
  dfw = n - k
  dfb = k - 1
  
  MSw = SSw/dfw
  MSb = SSb/dfb
  Fstat = MSb/MSw
  
  if (version=="hays1") {
    es = (SSb - dfb*MSw)/(SSt + MSw)
  }
  else if (version=="hays2") {
    es = ((dfw - 2)*Fstat/dfw - 1)/((dfw + 1)/dfb + (dfw - 2)*Fstat/dfw)
  }
  
  return(es)
}