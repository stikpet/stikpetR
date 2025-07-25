#' Cohen d' (for one-sample)
#' 
#' This function will calculate Cohen d' (one-sample). An effect size measure that can be used with a test for a single mean (for example a one-sample Student t-test).
#' 
#' This function is shown in this [YouTube video](https://youtu.be/Ihn-f386m-U) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CohenD.html)
#' 
#' 
#' @param data pandas series with the numeric scores
#' @param mu optional parameter to set the hypothesized mean. If not used the midrange is used
#' 
#' 
#' @return Cohen d'. mu is also printed if not provided.
#' 
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
#' *es_convert(d', from="cohendos", to="cohend")*, then use *th_cohen_d(d)*
#'
#' Or convert it further to an Odds Ratio using, *es_convert(d, from="cohend", to="or", ex1="chinn")* or *es_convert(d, from="cohend", to="or", ex1="borenstein")*. Then use *th_odds_ratio(or)*
#' 
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to perform a test:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' \code{\link{ts_z_os}}, for One-Sample Z-Test.
#' 
#' After this you might want a rule-of-thumb for the effect size, first convert to regular Cohen d:
#' \code{\link{es_convert}}, to convert Cohen's d one-sample to Cohen d, use *fr = "cohendos"* and *to = "cohend"*.
#' \code{\link{th_cohen_d}}, for rules-of-thumb for Cohen d.
#' 
#' Alternative Effect Sizes:
#' \code{\link{es_hedges_g_os}}, for Hedges g.
#' \code{\link{es_common_language_os}}, for the Common Language Effect Size.
#' 
#' 
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples
#' #Example 1: Numeric dataframe
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = df2['Gen_Age']
#' es_cohen_d_os(ex1)
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' es_cohen_d_os(ex2)
#' 
#' @export
es_cohen_d_os <- function(data, mu=NULL){
  data = unlist(na.omit(data))

  #set hypothesized mean to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2}

  xBar = mean(data)
  s = sd(data)
  dos = (xBar - mu)/s

  return(dos)

}



