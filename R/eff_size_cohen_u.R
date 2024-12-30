#' Cohen U
#' @description 
#' Cohen (1988, p. 23) provided three measures that relate to Cohen's d.
#' \itemize{
#' \item  \eqn{U_1}, is the proportion of non-overlap between distributions
#' \item  \eqn{U_2}, is the proportion of overlap between distributions
#' \item  \eqn{U_3}, is the proportion of one group's scores below the mean of another group
#' }
#' 
#' \eqn{U_1} and \eqn{U_2} are probably the least used of these three.
#' 
#' By converting each back to Cohen's d, the rule-of-thumb from Cohen d could be used as classification.
#' A nice interactive visualisation of the relation between Cohen $U_3$ and the Common Language Effect size, can be found at https://rpsychologist.com/therapist-effects/. 
#' 
#' @param d the Cohen d value
#' @param version {"u3", "u2", "u1"}, Optional, the version of Cohen U to determine
#' 
#' @returns 
#' The Cohen U value
#' 
#' @details
#' The following formulas are used (Cohen, 1988, p. 23):
#' \deqn{U_3 = \Phi\left(d\right)}
#' \deqn{U_2 = \Phi\left(\frac{d}{2}\right)}
#' \deqn{U_1 = \Phi\left(\frac{2\times U_2 - 1}{U_2}\right)}
#' 
#' *Symbols used:*
#' \itemize{
#'  \item  \eqn{d}, Cohen's d value
#' \item \eqn{n_i} the number of scores in category i
#' \item  \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' @seealso 
#' \code{\link{es_convert}}, to convert an U to Cohen d use `fr="cohenu.", to="cohend"`.
#' 
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_cohen_u <- function(d, version='u3'){
  if (version=='u3'){u = pnorm(d)}
  else if (version=='u2' || version=='u1'){
    u = pnorm(d/2)
    if (version=='u1'){u = (2*u - 1)/u}
  }
  
  return (u)
}