#' Hedges g (for Paired Samples)
#' 
#' @param var1 the scores on the first variable
#' @param var2 the scores on the second variable
#' @param dmu optional the difference according to null hypothesis (default is 0)
#' @param appr approximation to use (see details), default is "none"
#' @param within boolean to use a correction for correlated pairs
#' @returns 
#' A dataframe with:
#' \item{g}{the Hedges g value}
#' \item{version}{version that was used}
#' 
#' @details 
#' The formula used is (Hedges, 1981, p. 110):
#' \deqn{g = \frac{\bar{d}}{s_d}}
#' With:
#' \deqn{s_d = \sqrt{\frac{\sum_{i=1}^n\left(d_i - \bar{d}\right)^2}{n - 1}}}
#' \deqn{d_i = x_{i,1} - x_{i,2}}
#' \deqn{\bar{d} = \frac{\sum_{i=1}^n d_i}{n}}
#' 
#' **Symbols used:**
#' \itemize{
#' \item \eqn{n} the number of pairs (sample size)
#' \item \eqn{x_{i,1}} the i-th score of the first variable
#' \item \eqn{x_{i,2}} the i-th score of the second variable
#' }
#' 
#' This is the same as Cohen's d.
#' 
#' Hedges proposes the following exact bias correction (Hedges, 1981, p. 111):
#' \deqn{g_{c} = g \times\frac{\Gamma\left(m\right)}{\Gamma\left(m - \frac{1}{2}\right)\times\sqrt{m}}}
#' With:
#' \deqn{m = \frac{df}{2}}
#' \deqn{df = n - 1}
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
#' If *within=TRUE* the formula is changed to (Borenstein et al., 2009, p. 29):
#' \deqn{d_z = \frac{\bar{d}}{s_w}}
#' With:
#' \deqn{s_w = \frac{s_d}{\sqrt{2\times\left(1 - r_p\right)}}}
#' \deqn{r_p = \frac{\sum_{i=1}^n \left(x_{i,1} - \bar{x}_1\right) \times \left(x_{i,2} - \bar{x}_2\right)}{\left(n - 1\right)\times s_1\times s_2}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^n x_{i,j}}{n_i}}
#' 
#' 
#' **Alternatives**
#' 
#' *library(effsize)*
#' 
#' datF = na.omit(data.frame(var1, var2))
#' 
#' cohen.d(datF$var1, datF$var2, paired=TRUE, within=TRUE, hedges.correction=TRUE)
#' 
#' cohen.d(datF$var1, datF$var2, paired=TRUE, within=FALSE, hedges.correction=TRUE)
#' 
#' @examples 
#' var1 = c(8, 6, 20, 28, 60, 22, 26, 14, 30, 34, 36, 22, 10, NA, 96, 70, 62, 48, 38, 98, 82, 12, 70, 82, 90, 42)
#' var2 = c(0, 2, 2, 8, 12, 14, 14, 18, 18, 20, 22, 26, 32, 23, 32, 42, 44, NA, 48, 50, 52, 54, 54, 66, 68, 76)
#' es_hedges_g_ps(var1, var2, within=TRUE, appr="none")
#' es_hedges_g_ps(var1, var2, within=TRUE, appr="hedges")
#' es_hedges_g_ps(var1, var2, within=TRUE, appr="durlak")
#' es_hedges_g_ps(var1, var2, within=TRUE, appr="xue")
#' 
#' es_hedges_g_ps(var1, var2, within=FALSE, appr="none")
#' es_hedges_g_ps(var1, var2, within=FALSE, appr="hedges")
#' es_hedges_g_ps(var1, var2, within=FALSE, appr="durlak")
#' es_hedges_g_ps(var1, var2, within=FALSE, appr="xue")
#' 
#' 
#' @references 
#' Becker, B. J. (1988). Synthesizing standardized mean-change measures. *British Journal of Mathematical and Statistical Psychology, 41*(2), 257–278. https://doi.org/10.1111/j.2044-8317.1988.tb00901.x
#' 
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Effect sizes based on means. In *Introduction to Meta-Analysis* (pp. 21–32). John Wiley & Sons, Ltd. https://doi.org/10.1002/9780470743386
#' 
#' Durlak, J. A. (2009). How to select, calculate, and interpret effect sizes. *Journal of Pediatric Psychology, 34*(9), 917–928. https://doi.org/10.1093/jpepsy/jsp004
#' 
#' Hedges, L. V. (1981). Distribution theory for Glass’s estimator of effect size and related estimators. *Journal of Educational Statistics, 6*(2), 107–128. https://doi.org/10.2307/1164588
#' 
#' Xue, X. (2020). Improved approximations of Hedges’ g*. https://doi.org/10.48550/arXiv.2003.06675
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @export
es_hedges_g_ps <- function(var1, var2, dmu=0, appr=c("none", "hedges", "durlak", "xue"), within=TRUE){
  
  if (length(appr)>1) {
    appr="none"
  }
  
  datF = na.omit(data.frame(var1, var2))
  n = nrow(datF)
  
  d = datF$var1 - datF$var2
  dAvg = mean(d)
  
  s = sd(d)
  
  comment=""
  if (within) {
    s = s/sqrt(2*(1 - cor(datF$var1, datF$var2)))
    comment=", with correlation correction"
  }
  
  dz = (dAvg - dmu)/s
  df = n - 1
  m = df/2
  
  #check if exact chosen and can be calculated
  if (appr=="none" && m > 171){
    warning("Exact method could not be computed due to large sample size, Hedges approximation used instead")
    appr="hedges"    
  }
  
  
  if (appr=="none" && m <172){
    # Hedges g (exact)
    g = dz*gamma(m)/(gamma(m-0.5)*sqrt(m))
    comment = paste0("exact", comment)
  }
  else if(appr=="hedges"){
    # Hedges approximation
    g = dz*(1-3/(4*df-1))
    comment = paste0("Hedges approximation", comment)
  }
  else if(appr=="durlak"){
    # Durlak (2009, p. 927) approximation:
    g = dz*(n-3)/(n-2.25)*sqrt((n-2)/n)
    comment = paste0("Durlak approximation", comment)
  }
  else if(appr=="xue"){
    # Xue (2020, p. 3) approximation:
    g = dz*(1 - 9/df + 69/(2*df^2) - 72/(df^3) + 687/(8*df^4) - 441/(8*df^5) + 247/(16*df^6))^(1/12)
    comment = paste0("Xue approximation", comment)
  }

  #prepare results
  res <- data.frame(g, comment)
  colnames(res) <- c('g','version')
  
  return(res)
  
}