#' Common Language Effect Size (Paired-Samples)
#' @description 
#' the probability that a randomly selected score from the one population will be greater than a randomly sampled score from the other population.
#'
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param dmu difference according to null hypothesis (default is 0), only if method="mcgraw-wong"
#' @param method method to use for calculating CL (see details)
#' 
#' @returns
#' cl, float. the common language effect size measure value
#'
#' The formula used (McGraw & Wong, 1992, p. 363):
#' \deqn{CL = \Phi\left(z_{cl}\right)}
#' With:
#' \deqn{z_{cl} = \frac{\left|\bar{x}_1 - \bar{x}_2\right| - d_{H0}}{\sqrt{s_1^2 + s_2^2 - 2\times r_p \times s_1 \times s_2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^n x_{i,j}}{n_i}}
#' \deqn{r_p = \frac{\sum_{i=1}^n \left(x_{i,1} - \bar{x}_1\right) \times \left(x_{i,2} - \bar{x}_2\right)}{\left(n - 1\right)\times s_1\times s_2}}
#'
#' **Symbols used:**
#' \itemize{
#' \item \eqn{n} the total number of pairs
#' \item \eqn{x_{i,j}} the i-th score in the j-th variable
#' \item \eqn{r_p} the Pearson correlation coefficient
#' }
#'
#' This equation is used when method="mcgraw-wong"
#'
#' The formula used for the Dunlap method (Dunlap, 1994, p. 509):
#' \deqn{CL = \sin^{-1}\left(r\right) + \frac{1}{2}}
#'
#' This equation is used when method="dunlap".
#'
#' @references
#' Dunlap, W. P. (1994). Generalizing the common language effect size indicator to bivariate normal correlations. *Psychological Bulletin, 116*(3), 509-511. https://doi.org/10.1037/0033-2909.116.3.509
#'
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. *Psychological Bulletin, 111*(2), 361-365. https://doi.org/10.1037/0033-2909.111.2.361
#'
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
es_common_language_ps <- function(field1, field2, dmu=0, method=c("dunlap", "mcgraw-wong")){

  if (length(method)>1) {
    method="dunlap"
  }

  datF = na.omit(data.frame(field1, field2))

  n = nrow(datF)

  mx = mean(datF$field1)
  my = mean(datF$field2)

  sx = sd(datF$field1)
  sy = sd(datF$field2)

  sxy = sum((datF$field1 - mx)*(datF$field2 - my))

  r = sxy/((n - 1)*sx*sy)

  if (method=="dunlap") {
    cl = asin(r)/pi + 0.5
  }
  else if (method=="mcgraw-wong"){
    se = sqrt(sx**2 + sy**2 - 2*r*sx*sy)
    z = (abs(mx - my) - dmu)/se
    cl = pnorm(z)
  }


  return(cl)

}



