#' One-Sample Z Test
#' 
#' @description 
#' This test is often used if there is a large sample size. For smaller sample sizes, a Student t-test
#' is usually used.
#' 
#' The assumption about the population (null hypothesis) for this test is a pre-defined mean, i.e. the (arithmetic) mean 
#' that is expected in the population. If the p-value (significance) is then below a pre-defined threhold 
#' (usually 0.05), the assumption is rejected.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/zq0IjOq-Jzo) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/zOneSample.html)
#' 
#' @param data A vector or dataframe with the data as numbers
#' @param mu optional hypothesized mean, otherwise the midrange will be used
#' @param sigma population standard deviation, if NULL sample results will be used
#' 
#' @returns
#' A dataframe with:
#' \item{mu}{the hypothesized mean}
#' \item{sample.mean}{the sample mean}
#' \item{statistic}{the test statistic}
#' \item{pValue}{the significance (p-value)}
#' \item{testUsed}{name of test used}
#'
#' @details
#' The formula used is:
#' \deqn{z = \frac{\bar{x} - \mu_{H_0}}{SE}}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#'
#' With:
#' \deqn{SE = \frac{\sigma}{\sqrt{n}}}
#' \deqn{\sigma \approx s = \sqrt{\frac{\sum_{i=1}^n\left(x_i - \bar{x}\right)^2}{n - 1}}}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^nx_i}{n}}
#'
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\Phi\left(\dots\right)} the cumulative distribution function of the standard normal distribution
#' \item \eqn{\bar{x}} the sample mean
#' \item \eqn{\mu_{H_0}} the hypothesized mean in the population
#' \item \eqn{SE} the standard error (i.e. the standard deviation of the sampling distribution)
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{s} the unbiased sample standard deviation
#' \item \eqn{x_i} the i-th score
#' }
#' 
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to create a binned frequency table or a visualisation:
#' \code{\link{tab_frequency_bins}}, to create a binned frequency table.
#' \code{\link{vi_boxplot_single}}, for a Box (and Whisker) Plot.
#' \code{\link{vi_histogram}}, for a Histogram.
#' \code{\link{vi_stem_and_leaf}}, for a Stem-and-Leaf Display.
#' 
#' After this you might want an effect size measure:
#' \code{\link{es_cohen_d_os}}, for for Cohen d'.
#' \code{\link{es_hedges_g_os}}, for Hedges g.
#' \code{\link{es_common_language_os}}, for the Common Language Effect Size.
#' 
#' Alternative Tests:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples 
#' #Example 1: dataframe
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = df2['Gen_Age']
#' ts_z_os(ex1)
#' ts_z_os(ex1, mu=22, sigma=12.1)
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' ts_z_os(ex2)
#' 
#' @export
ts_z_os <- function(data, mu=NULL, sigma=NULL){
  data = data.frame(data)
  data = na.omit(data)

  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }

  n = nrow(data)
  data = as.numeric(data[,1])
  m = mean(data)

  if (is.null(sigma)) {
    s =sd(data)
  }
  else{
    s = sigma
  }

  se = s/sqrt(n)
  z = (m - mu)/se

  pValue = 2*(1 - pnorm(abs(z)))

  sample.mean = m
  statistic = z
  testUsed = "one-sample Z"

  testResults <- data.frame(mu, sample.mean, statistic, pValue, testUsed)

  return(testResults)
}



