#' One-Sample Student t-Test
#' 
#' A test for a single (arithmetic) mean. 
#' 
#' The assumption about the population (null hypothesis) for this test is a pre-defined mean, i.e. the (arithmetic) mean that is expected in the population. If the p-value (significance) is then below a pre-defined threhold (usually 0.05), the assumption is rejected.
#'
#'
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized mean, otherwise the midrange will be used
#' @return dataframe with the hypothesized mean, sample mean, test statistic, degrees of freedom, p-value (sig.) and name of test used
#'
#' @examples
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' ts_student_t_os(data)
#'
#' @details
#' The formula used is:
#' \deqn{t = \frac{\bar{x} - \mu_{H_0}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#'
#' With:
#' \deqn{df = n - 1}
#' \deqn{SE = \frac{s}{\sqrt{n}}}
#' \deqn{s = \sqrt{\frac{\sum_{i=1}^n\left(x_i - \bar{x}\right)^2}{n - 1}}}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^nx_i}{n}}
#'
#' *Symbols used:*
#' \itemize{
#' \item \eqn{T\left(\dots, \dots\right)} the cumulative distribution function of the t-distribution
#' \item \eqn{\bar{x}} the sample mean
#' \item \eqn{\mu_{H_0}} the hypothesized mean in the population
#' \item \eqn{SE} the standard error (i.e. the standard deviation of the sampling distribution)
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{s} the unbiased sample standard deviation
#' \item \eqn{x_i} the i-th score
#' }
#'
#' The Student t test (Student, 1908) was described by Gosset under the pseudo name Student.
#'
#' **Alternative**
#' R stats library has a similar function: *t.test()*
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references
#' Student. (1908). The probable error of a mean. *Biometrika, 6*(1), 1–25. https://doi.org/10.1093/biomet/6.1.1
#'
#' @export
ts_student_t_os <- function(data, mu=NULL){
  data = data.frame(data)
  data = na.omit(data)

  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }

  n = nrow(data)
  m = mean(data[,1])
  s =sd(data[,1])
  se = s/sqrt(n)
  t = (m - mu)/se

  df = n - 1

  pValue = 2*(1 - pt(abs(t), df))

  sample.mean = m
  statistic = t
  testUsed = "one-sample Student t"

  testResults <- data.frame(mu, sample.mean, statistic, df, pValue, testUsed)
  colnames(testResults)<-c("mu", "sample mean", "statistic", "df", "p-value", "test used")

  return(testResults)
}
