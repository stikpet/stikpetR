#' Hedges g / Cohen ds (independent samples)
#' 
#' @param scores A vector with the scores data
#' @param groups A vector with the group data
#' @param appr c(NULL, 'hedges-exact', 'hedges-appr', 'durlak', 'xue') approximation to use, NULL will use biased version.
#' @return Hedges g value
#' 
#' @details
#' The formula used is (Hedges, 1981, p. 110):
#' \deqn{g = \frac{\bar{x}_1 - \bar{x}_2}{s_p}}
#' 
#' With:
#' \deqn{s_p = \frac{\sqrt{s_1^2 + s_2^2}}{n - 2}}
#' \deqn{SS_i = \sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' Hedges proposes the following exact bias correction (Hedges, 1981, p. 111):
#' \deqn{g_{c} = g \times\frac{\Gamma\left(m\right)}{\Gamma\left(m - \frac{1}{2}\right)\times\sqrt{m}}}
#' With:
#' \deqn{m = \frac{df}{2}}
#' \deqn{df = n_1 + n_2 - 2= n - 2}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{\Gamma\left(\dots\right)} the gamma function
#' }
#' 
#' The formula used for the approximation for this correction from Hedges (1981, p. 114) (appr="hedges"):
#' \deqn{g_c = g \times\left(1 - \frac{3}{4\times df - 1}\right)}
#' 
#' This approximation can also be found in Hedges and Olkin (1985, p. 81) and
#' Cohen (1988, p. 66)
#' 
#' The formula used for the approximation from Durlak (2009, p. 927) (appr="durlak"):
#' \deqn{g_c = g \times\frac{n - 3}{n - 2.25} \times\sqrt{\frac{n - 2}{n}}}
#' 
#' The formula used for the approximation from Xue (2020, p. 3) (appr="xue"):
#' \deqn{g_c = g \times \sqrt[12]{1 - \frac{9}{df} + \frac{69}{2\times df^2} - \frac{72}{df^3} + \frac{687}{8\times df^4} - \frac{441}{8\times df^5} + \frac{247}{16\times df^6}}}
#' 
#' 
#' Cohen (1988, pp. 66-67) d_s is the same as Hedges g (biased).
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Durlak, J. A. (2009). How to select, calculate, and interpret effect sizes. *Journal of Pediatric Psychology, 34*(9), 917–928. https://doi.org/10.1093/jpepsy/jsp004
#' 
#' Hedges, L. V. (1981). Distribution Theory for Glass’s Estimator of Effect Size and Related Estimators. *Journal of Educational Statistics, 6*(2), 107–128. https://doi.org/10.2307/1164588
#' 
#' Hedges, L. V., & Olkin, I. (1985). *Statistical methods for meta-analysis*. Academic Press.
#' 
#' Xue, X. (2020). Improved approximations of Hedges’ g*. https://doi.org/10.48550/arXiv.2003.06675
#' 
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' es_hedges_g_is(scores, groups)
#' 
#' @export
es_hedges_g_is <- function(scores, groups, corr = "none"){
  #make sure data is numeric
  scores = as.numeric(scores)
  
  #remove rows with missing values
  df = data.frame(scores, groups)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  X1 = df$score[df$group == df$group[1]]
  X2 = df$score[df$group != df$group[1]]
  
  var1 = var(X1)
  var2 = var(X2)
  
  mean1 = mean(X1)
  mean2 = mean(X2)
  
  n = length(df$group)
  n1 = sum(df$group==df$group[1])
  n2 = n - n1
  
  sd1 = sqrt(var1)
  sd2 = sqrt(var2)
  
  #Determine Sum of Squared (deviation from the mean) per category
  ss1 = sd1^2*(n1 - 1)
  ss2 = sd2^2*(n2 - 1)
  
  #Determine Hedges g (Cohen's d)
  n = n1 + n2
  g = (mean1 - mean2)/sqrt((ss1 + ss2)/(n - 2))
  
  c = 1
  
  if (corr=="hedges-exact") {
    if (n - 2 < 171) {
      c = gamma((n - 2)/2)/(sqrt((n - 2)/2)*gamma((n - 3)/2))
    }
    else {
      print("WARNING: exact method could not be computed due to large sample size, approximation used instead")
      c = 1 - 3/(4*(n - 2) - 9)
    }
  }
  else if(corr=="hedges-appr") {
    c = 1 - 3/(4*(n - 2) - 9)
  }
  else if(corr=="durlak") {
    c = (n - 3)/(n - 2.25)*sqrt((n - 2)/n)
  }
  else if(corr=="xue"){
    # Xue (2020, p. 3) approximation:
    df = n - 2
    c = (1 - 9/df + 69/(2*df^2) - 72/(df^3) + 687/(8*df^4) - 441/(8*df^5) + 247/(16*df^6))^(1/12)
    comment = "Xue approximation"
  }
  
  g = g*c
  
  return(g)
  
}