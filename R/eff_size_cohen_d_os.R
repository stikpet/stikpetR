#' Cohen d' (for one-sample)
#'
#' @param data pandas series with the numeric scores
#' @param mu optional parameter to set the hypothesized mean. If not used the midrange is used
#' @return Cohen d'. mu is also printed if not provided.
#'
#' @examples
#' data <- c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, 90, 25, 40, 70, 65, 70, 98, 40, 65, 60)
#' es_cohen_d_os(data)
#' es_cohen_d_os(data, mu=70)
#'
#' @details
#' The formula used (Cohen, 1988, p. 46):
#' \deqn{d'=\frac{\bar{x}-\mu_{H_{0}}}{s}}
#' With:
#' \deqn{s = \sqrt{\frac{\sum_{i=1}^n \left(x_i - \bar{x}\right)^2}{n - 1}}}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^n x_i}{n}}
#'
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\bar{x}} the sample mean
#' \item \eqn{\mu_{H_0}} the hypothesized mean in the population
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{s} the unbiased sample standard deviation
#' \item \eqn{x_i} the i-th score
#' }
#'
#' Note to use a rule-of-thumb from Cohen d, first convert this to a regular Cohen d using
#' *es_convert(d', from="cohendos", to="cohend")*
#'
#' Then use *th_cohen_d(d)*
#'
#' Or convert it further to an Odds Ratio using:
#' *es_convert(d, from="cohend", to="or", ex1="chinn)* or *es_convert(d, from="cohend", to="or", ex1="borenstein)*
#'
#' Then use *th_odds_ratio(or)*
#'
#' **Alternative**
#'
#' The *lsr* library has a similar function: *cohensD()*
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'
#' @export
es_cohen_d_os <- function(data, mu=NULL){

  #set hypothesized mean to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
    print("Hypothesized mean used:")
    print(mu)}

  xBar = mean(data)
  s = sd(data)
  dos = (xBar - mu)/s

  return(dos)

}
