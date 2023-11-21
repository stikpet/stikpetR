#' Student t Test (Paired Samples)
#' @description 
#' The assumption about the population (null hypothesis) for this test is a pre-defined difference between two means, usually zero (i.e. the difference between the (arithmetic) means is zero, they are the same in the population). If the p-value (significance) is then below a pre-defined threhold (usually 0.05), the assumption is rejected.
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param dmu difference according to null hypothesis (default is 0)
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the number of scores}
#' \item{statistic}{the test statistic (t-value)}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used is:
#' \deqn{t_p = \frac{\bar{d} - d_{H0}}{SE}}
#' \deqn{sig. = 2\times\left(1 - \text{T}\left(\left|t_p\right|, df\right)\right)}
#' With:
#' \deqn{\bar{d} = \bar{x}_1 - \bar{x}_2}
#' \deqn{SE = \sqrt{\frac{\sigma_{s^2}}{n}}}
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
#' \item \eqn{\text{T}\left(\dots, \dots\right)} the cumulative distribution function of the Student t distribution
#' }
#' 
#' **Alternatives**
#' 
#' R's *stats* library
#' 
#' t.test(var1, var2, paired=TRUE)
#' 
#' t.test(var1, var2, paired=TRUE, mu=5)
#' 
#' @references 
#' Student. (1908). The probable error of a mean. *Biometrika, 6*(1), 1–25. doi:10.1093/biomet/6.1.1
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_student_t_ps <- function(field1, field2, dmu=0){
  
  datF = na.omit(data.frame(field1, field2))
  
  d = datF$field1 - datF$field2
  dAvg = mean(d)
  s = sd(d)
  n = nrow(datF)
  se = s/sqrt(n)
  
  t = (dAvg - dmu)/se
  
  df = n - 1
  
  pValue = 2*(1 - pt(abs(t), df, lower.tail = TRUE))
  
  results = data.frame(n, t, df, pValue)
  colnames(results) = c("n", "statistic", "df", "p-value")
  return(results)
  
}