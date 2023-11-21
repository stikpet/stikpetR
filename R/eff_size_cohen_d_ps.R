#' Cohen d_z (for Paired Samples)
#' @description 
#' An effect size measure for paired samples. 
#'
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param within boolean to use a correction for correlated pairs
#' 
#' @return the Cohen d effect size
#'
#' @details
#' the formula used (Cohen, 1988, p. 48):
#' \deqn{d_z = \frac{\bar{d}}{s_d}}
#' 
#' With:
#' \deqn{s_d = \sqrt{\frac{\sum_{i=1}^n\left(d_i - \bar{d}\right)^2}{n - 1}}}
#' \deqn{d_i = x_{i,1} - x_{i,2}}
#' \deqn{\bar{d} = \frac{\sum_{i=1}^n d_i}{n}}
#'
#' **Symbols used:**
#' \itemize{
#' \item \eqn{n} the number of pairs (sample size)
#' \item \eqn{x_{i,1}} the i-th score of the first variable
#' \item \eqn{x_{i,2}} the i-th score of the second variable
#' }
#'
#' If *within=TRUE* the formula is changed to (Borenstein et al., 2009, p. 29):
#' \deqn{d_z = \frac{\bar{d}}{s_w}}
#' With:
#' \deqn{s_w = \frac{s_d}{\sqrt{2\times\left(1 - r_p\right)}}}
#' \deqn{r_p = \frac{\sum_{i=1}^n \left(x_{i,1} - \bar{x}_1\right) \times \left(x_{i,2} - \bar{x}_2\right)}{\left(n - 1\right)\times s_1\times s_2}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^n x_{i,j}}{n_i}}
#'
#'
#' **Alternatives**
#'
#' library(lsr)
#'
#' cohensD(var1, var2, method="paired")
#'
#' *library(effsize)*
#'
#' datF = na.omit(data.frame(var1, var2))
#'
#' cohen.d(datF$var1, datF$var2, paired=TRUE)
#'
#' cohen.d(datF$var1, datF$var2, paired=TRUE, within=FALSE)
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Effect sizes based on means. In *Introduction to Meta-Analysis* (pp. 21–32). John Wiley & Sons, Ltd. doi:10.1002/9780470743386
#'
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#'
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
es_cohen_d_ps <- function(field1, field2, within=TRUE){

  datF = na.omit(data.frame(field1, field2))

  d = datF$field1 - datF$field2
  dAvg = mean(d)

  s = sd(d)

  if (within) {

    n = nrow(datF)

    mx = mean(datF$field1)
    my = mean(datF$field2)

    sx = sd(datF$field1)
    sy = sd(datF$field2)

    sxy = sum((datF$field1 - mx)*(datF$field2 - my))

    r = sxy/((n - 1)*sx*sy)
    s = s/sqrt(2*(1 - r))
  }

  dz = dAvg/s

  return(dz)

}
