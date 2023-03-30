#' Cochran Q Test
#' 
#' @param data dataframe with the binary scores
#' @param success indicator for what is considered a success (default is 1)
#' @returns 
#' A dataframe with:
#' \item{statistic}{the test statistic}
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
#' @examples 
#' var1 = c(0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1)
#' var2 = c(0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
#' var3 = c(0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
#' var4 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1)
#' dFr = data.frame(var1, var2, var3, var4)
#' ts_cochran_q(dFr)
#' 
#' @references 
#' Cochran, W. G. (1950). The comparison of percentages in matched samples. *Biometrika, 37*(3/4), 256–266. https://doi.org/10.2307/2332378
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @export
ts_cochran_q <- function(data, success = 1){
  dFr = na.omit(data)  
  k = ncol(dFr)
  nf = nrow(dFr)
  
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
  results = data.frame(statistic, df, pValue)
  
  return(results)

}