#' One-Sample Z t-Test
#' 
#' This test is often used if there is a large sample size. For smaller sample sizes, a Student t-test
#' is usually used.
#' 
#' The assumption about the population (null hypothesis) for this test is a pre-defined mean, i.e. the (arithmetic) mean 
#' that is expected in the population. If the p-value (significance) is then below a pre-defined threhold 
#' (usually 0.05), the assumption is rejected.
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized mean, otherwise the midrange will be used
#' @param sigma population standard deviation, if NULL sample results will be used
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
#' @examples
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' ts_z_os(data)
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @export
ts_z_os <- function(data, mu=NULL, sigma=NULL){

  data = na.omit(data)

  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }

  n = length(data)
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
