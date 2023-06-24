#' One-Sample (Yuen or Yuen-Welch) Trimmed Mean Test
#' 
#' A variation on a one-sample Student t-test where the data is first trimmed, and the Winsorized variance is used.
#' 
#' The assumption about the population for this test is that the mean in the population is equal to the provide mu value. The test will show the probability of the found test statistic, or more extreme, if this assumption would be true. If this is below a specific threshold (usually 0.05) the assumption is rejected.
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized trimmed mean, otherwise the midrange will be used
#' @param trim optional proportion to trim from each side (so in total twice this will be trimmed)
#' @param se c("yuen", "wilcox") optional method to use to determine standard error (default is "yuen")
#' @return dataframe with the sample trimmed mean, hypothesized trimmed mean, the standard error, test statistic, degrees of freedom, p-value (sig.) and name of test used
#' 
#' @examples  
#' grade = c(4, 10, 2, 9, 5, 28, 8, 7, 9, 35, 40, 12, 8, 6, 16, 12, 14, 10, 18, 4, 11)
#' ts_trimmed_mean_os(grade, trim=0.1, mu=12)
#' ts_trimmed_mean_os(grade, trim=0.1, mu=12, se="wilcox")
#' 
#' @details 
#' The formula used is:
#' \deqn{\frac{\bar{x}_t - \mu_{H_0}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{\bar{x}_t = \frac{\sum_{i=g+1}^{n - g}y_i}{}}
#' \deqn{g = \lfloor n\times p_t\rfloor}
#' \deqn{m = n - 2\times g}
#' \deqn{SE = \sqrt{\frac{SSD_w}{m}\times\left(m - 1\right)}}
#' or:
#' \deqn{SE = \frac{\sqrt{SSD_w}}{\left(1 - 2\times p_t\right)\times\sqrt{n}}}
#' \deqn{SSD_w = g\times\left(y_{g+1} - \bar{x}_w\right)^2 + g\times\left(y_{n-g} - \bar{x}_w\right)^2 + \sum_{i=g+1}^{n - g} \left(y_i - \bar{x}_w\right)^2}
#' \deqn{\bar{x}_w = \frac{\bar{x}_t\times m + g\times\left(y_{g+1} + y_{n-g}\right)}{n}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_t} the trimmed mean of the scores
#' \item \eqn{x_w} The Winsorized mean
#' \item \eqn{SSD_w} the sum of squared deviations from the Winsorized mean 
#' \item \eqn{m} the number of scores in the trimmed data set from category i
#' \item \eqn{y_i} the i-th score after the scores are sorted from low to high
#' \item \eqn{p} the proportion of trimming on each side, we can define
#' }
#' 
#' The test is often also referred to as a Yuen test, or Yuen-Welch test.
#' 
#' The standard error can either be calculated using the first SE, which for example can be found in
#' Tukey and McLaughlin (1963, p. 342), and seems similar to the independent samples version of this 
#' test as proposed by Yuen (1974, p. 167)
#' 
#' The second version is used in the other libraries, and can be found in Wilcox (2012, p. 157), 
#' or Peró-Cebollero and Guàrdia-Olmos (2013, p. 409).
#' 
#' 
#' **Alternatives**
#' R stats library has a similar function: *t.test()*
#' 
#' The *DescTools* library has a similar function: *YuenTTest()*
#' 
#' The *PairedData* library has a similar functions: *yuen.t.test()* and *yuen1.test*
#' 
#' The *WRS2* library has a similar function: *yuen()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Peró-Cebollero, M., & Guàrdia-Olmos, J. (2013). The adequacy of different robust statistical tests in comparing two independent groups. *Psicológica*, 34, 407–424.
#' 
#' Tukey, J. W., & McLaughlin, D. H. (1963). Less vulnerable confidence and significance procedures for location based on a single sample: Trimming/Winsorization 1. *Sankhyā: The Indian Journal of Statistics, 25*(3), 331–352.
#' 
#' Wilcox, R. R. (2012). *Introduction to robust estimation and hypothesis testing* (3rd ed.). Academic Press.
#' 
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances. *Biometrika, 61*(1), 165–170. https://doi.org/10.1093/biomet/61.1.165
#'  
#' @export
ts_trimmed_mean_os <- function(data, mu=NULL, trim=0.1, se="yuen"){
  
  data = na.omit(data)
  
  #set mu to midrange if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }
  
  #untrimmed sample size (n):
  n = length(data)
  
  #number of scores not trimmed:
  nt = n - 2*round(n*trim)
  
  #trimmed mean
  mt = mean(data, trim=trim)
  
  #Winsorize the data
  scSort = sort(data)
  minReplace1 = scSort[round(n*trim)+1]
  maxReplace1 = scSort[n - round(n*trim)]
  data = replace(data, data<minReplace1, minReplace1)
  data = replace(data, data>maxReplace1, maxReplace1)
  
  #Winsorized variance
  var = var(data)
  
  if (se=="yuen"){
    se = sqrt(var)*sqrt(n - 1)/sqrt(nt*(nt - 1))
  }
  else if (se=="wilcox") {
    se = sqrt(var)/((1 - 2*trim)*sqrt(n))
  }
  
  trim.mean = mt
  t = (mt - 12)/se
  df = nt - 1
  pValue = 2*(1 - pt(abs(t), df))
  
  statistic=t
  testUsed = "one-sample trimmed mean test"
  results = data.frame(trim.mean, mu, se, statistic, df, pValue, testUsed)
  
  return(results)
}