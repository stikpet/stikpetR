#' Paired Samples Z Test
#' 
#' @param var1 the scores on the first variable
#' @param var2 the scores on the second variable
#' @param dmu difference according to null hypothesis (default is 0)
#' @param dsigma population standard deviation of the difference, if NULL sample results will be used
#' @returns 
#' A dataframe with:
#' \item{statistic}{the test statistic}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' 
#' The formula used is:
#' \deqn{z_p = \frac{\bar{d} - d_{H0}}{SE}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_p\right|\right)\right)}
#' With:
#' \deqn{\bar{d} = \mu_1 - \mu_2 \approx \bar{x}_1 - \bar{x}_2}
#' \deqn{SE = \sqrt{\frac{\sigma_d^2}{n}} \approx \sqrt{\frac{\sigma_s^2}{n}}}
#' \deqn{s_d^2 = \frac{\sum_{i=1}^n \left(d_i - \bar{d}_i\right)^2}{n - 1}}
#' \deqn{d_i = x_{i,1} - x_{i,2}}
#' \deqn{\bar{d} = \frac{\sum_{i=1}^n d_i}{n}}
#' 
#' **Symbols used:**
#' \itemize{
#' \item \eqn{n} the number of pairs (sample size)
#' \item \eqn{x_{i,1}} the i-th score of the first variable
#' \item \eqn{x_{i,2}} the i-th score of the second variable
#' \item \eqn{d_{H0}} the expected difference in the population
#' }
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' dfr = na.omit(data.frame(var1, var2))
#' 
#' ZTest(dfr$var1, dfr$var2, sd_pop=sqrt(var(dfr$var1-dfr$var2)), paired=TRUE)
#' 
#' @examples 
#' var1 = c(8, 6, 20, 28, 60, 22, 26, 14, 30, 34, 36, 22, 10, NA, 96, 70, 62, 48, 38, 98, 82, 12, 70, 82, 90, 42)
#' var2 = c(0, 2, 2, 8, 12, 14, 14, 18, 18, 20, 22, 26, 32, 23, 32, 42, 44, NA, 48, 50, 52, 54, 54, 66, 68, 76)
#' ts_z_ps(var1, var2)
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @export
ts_z_ps <- function(var1, var2, dmu=0, dsigma=NULL){
  
  datF = na.omit(data.frame(var1, var2))
  
  d = datF$var1 - datF$var2
  dAvg = mean(d)
  
  if (is.null(dsigma)) {
    dsigma = sd(d)
  }
  
  n = nrow(datF)
  se = dsigma/sqrt(n)
  
  z = (dAvg - dmu)/se
  
  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic=z
  results = data.frame(statistic, pValue)
  
  return(results)

}