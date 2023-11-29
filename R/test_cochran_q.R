#' Cochran Q Test
#' @description 
#' A test for multiple binairy variables. The null hypothesis is that the proportion of successes is the same in all groups.
#' 
#' If the p-value (sig.) is below a certain threshold (usually .05) the assumption is rejected and at least one category has a significant different number of successes than at least one other group, in the population.
#' 
#' If the test is significant (below the threshold) a post-hoc Dunn test could be used, or pairwise McNemar-Bowker.

#' @param data dataframe with the binary scores
#' @param success indicator for what is considered a success (default is first value found)
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{the test statistic (chi-square value)}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Cochran, 1950, p. 259):
#' \deqn{Q = \frac{\left(k - 1\right)\times\sum_{j=1}^k\left(C_j - \bar{C}\right)^2}{k\times\sum_{i=1}^n R_i - \sum_{i=1}^n R_i^2}}
#' \deqn{sig. = 1 - \chi^2\left(Q, df\right)}
#' With:
#' \deqn{df = k - 1}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{C_j} the number of successes in category j
#' \item \eqn{k} the number of categories (factors)
#' \item \eqn{R_i} the number of succeses in case i 
#' \item \eqn{n} the number of cases
#' }
#' 
#' **Alternatives**
#' 
#' *library(nonpar)*
#' 
#' matr = cbind(var1, var2, var3, var4)
#' 
#' cochrans.q(matr)
#' 
#' *library(RVAideMemoire)*
#' 
#' myData.long<-reshape(dFr, varying=c("var1", "var2", "var3", "var4"), v.names="score", timevar="var", times=c("var1", "var2", "var3", "var4"),new.row.names = 1:1000, direction="long")
#' 
#' cochran.qtest(score~var |id, data=myData.long)
#' 
#' @references 
#' Cochran, W. G. (1950). The comparison of percentages in matched samples. *Biometrika, 37*(3/4), 256–266. https://doi.org/10.2307/2332378
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_cochran_q <- function(data, success = NULL){
  dFr = na.omit(data)  
  k = ncol(dFr)
  nf = nrow(dFr)
  
  if (is.null(success)){
    success=dFr[1,1]}
  
  C = rep(0, k)
  for (i in 1:k) {
    C[i] = sum(dFr[i]==success)
  }
  
  nSuc = sum(C)
  
  R = rep(0, nf)
  for (i in 1:nf) {
    R[i] = sum(dFr[i,]==success)
  }
  
  Q = k*(k - 1)*sum((C - mean(C))**2)/(k*sum(R) - sum(R**2))
  
  df = k - 1
  
  pValue = pchisq(Q, df, lower.tail = FALSE)
  
  statistic = Q
  results = data.frame(nf, statistic, df, pValue)
  colnames(results) = c("n", "statistic", "df", "p-value")
  
  return(results)
  
}