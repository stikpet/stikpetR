#' Hodges-Lehmann Estimate (One-Sample)
#' @description 
#' The Hodges-Lehmann Estimate (Hodges & Lehmann, 1963) for a one-sample scenario, is the median of the Walsh averages. The Walsh averages (Walsh, 1949a, 1949b) are the average of each possible pair by taking one score and combining it with each of the other scores. Note that each is only counted once, so taking the second and fifth score is the same as taking the fifth and the second, so only one of these is used. It does also include self-pairs, e.g. the third score and third score.
#' 
#' It is in the one-sample case therefor a measure of central tendancy and sometimes referred to as the pseudo median.
#' 
#' @param scores list with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' 
#' @returns 
#' HL : float, the Hodges-Lehmann Estimate
#' 
#' @details
#' The formula used (Hodges & Lehmann, 1963, p. 599):
#' \deqn{HL = \text{median}\left(\frac{x_i + x_j}{2} | i \leq i \leq j \leq n\right)}
#' 
#' @references
#' Hodges, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests. *The Annals of Mathematical Statistics, 34*(2), 598–611. doi:10.1214/aoms/1177704172
#' 
#' Monahan, J. F. (1984). Algorithm 616: Fast computation of the Hodges-Lehmann location estimator. *ACM Transactions on Mathematical Software, 10*(3), 265–270. doi:10.1145/1271.319414
#' 
#' Walsh, J. E. (1949a). Applications of some significance tests for the median which are valid under very general conditions. Journal of the American Statistical Association, 44(247), 342–355. doi:10.1080/01621459.1949.10483311
#' 
#' Walsh, J. E. (1949b). Some significance tests for the median which are valid under very general conditions. *The Annals of Mathematical Statistics, 20*(1), 64–81. doi:10.1214/aoms/1177730091
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
me_hodges_lehmann_os <- function(scores, levels=NULL){
  if (!is.null(levels)){
    scores = factor(na.omit(scores), ordered = TRUE, levels = levels)
    scores = as.numeric(scores)
  }
  
  walsh = outer(scores, scores, "+") / 2
  HL = median(walsh[lower.tri(walsh, diag = TRUE)])  
  return (HL)
}