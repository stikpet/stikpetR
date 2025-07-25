#' Z-test (Paired Samples)
#' @description 
#' This test is often used if there is a large sample size. For smaller sample sizes, a Student t-test is usually used.
#' 
#' The assumption about the population (null hypothesis) for this test is a pre-defined difference between two means, usually zero (i.e. the difference between the (arithmetic) means is zero, they are the same in the population). If the p-value (significance) is then below a pre-defined threhold (usually 0.05), the assumption is rejected.
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param dmu difference according to null hypothesis (default is 0)
#' @param dsigma population standard deviation of the difference, if NULL sample results will be used
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the number of scores}
#' \item{z}{the test statistic (z-value)}
#' \item{p-Value}{the significance (p-value)}
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
#' \item \eqn{\Phi\left(\dots\right)}, cumulative density function of the standard normal distribution.
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
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_z_ps <- function(field1, field2, dmu=0, dsigma=NULL){
  
  datF = na.omit(data.frame(field1, field2))
  
  d = datF$field1 - datF$field2
  dAvg = mean(d)
  
  if (is.null(dsigma)) {
    dsigma = sd(d)
  }
  
  n = nrow(datF)
  se = dsigma/sqrt(n)
  
  z = (dAvg - dmu)/se
  
  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic=z
  results = data.frame(n, statistic, pValue)
  colnames(results) = c("n", "z", "p-value")
  
  return(results)
  
}



