#' Hodges-Lehmann Estimate (One-Sample)
#' @description 
#' The Hodges-Lehmann Estimate (Hodges & Lehmann, 1963) for a one-sample scenario, is the median of the Walsh averages. The Walsh averages (Walsh, 1949a, 1949b) are the average of each possible pair by taking one score and combining it with each of the other scores. Note that each is only counted once, so taking the second and fifth score is the same as taking the fifth and the second, so only one of these is used. It does also include self-pairs, e.g. the third score and third score.
#' 
#' It is in the one-sample case therefor a measure of central tendancy and sometimes referred to as the pseudo median.
#' 
#' The measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Measures/HodgesLehmannOS.html)
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
#' @section Before, After and Alternatives:
#' Before this measure you might want an impression using a frequency table or a visualisation:
#' \code{\link{tab_frequency}}, for a frequency table
#' \code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
#' \code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.
#' 
#' After this you might want some other descriptive measures:
#' \code{\link{me_consensus}}, for the Consensus. 
#' \code{\link{me_median}}, for the Median.
#' \code{\link{me_quantiles}}, for Quantiles.
#' \code{\link{me_quartiles}}, for Quartiles / Hinges.
#' \code{\link{me_quartile_range}}, for Interquartile Range, Semi-Interquartile Range and Mid-Quartile Range.
#' 
#' or perform a test:
#' \code{\link{ts_sign_os}}, for One-Sample Sign Test.
#' \code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
#' \code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
#' 
#' @references
#' Hodges, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests. *The Annals of Mathematical Statistics, 34*(2), 598-611. doi:10.1214/aoms/1177704172
#' 
#' Monahan, J. F. (1984). Algorithm 616: Fast computation of the Hodges-Lehmann location estimator. *ACM Transactions on Mathematical Software, 10*(3), 265-270. doi:10.1145/1271.319414
#' 
#' Walsh, J. E. (1949a). Applications of some significance tests for the median which are valid under very general conditions. Journal of the American Statistical Association, 44(247), 342-355. doi:10.1080/01621459.1949.10483311
#' 
#' Walsh, J. E. (1949b). Some significance tests for the median which are valid under very general conditions. *The Annals of Mathematical Statistics, 20*(1), 64-81. doi:10.1214/aoms/1177730091
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' # Example 1: Dataframe
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' studentDf = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = studentDf[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' me_hodges_lehmann_os(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' me_hodges_lehmann_os(ex2)
#' 
#' #Example 3: Text data with
#' ex3 = c("a", "b", "f", "d", "e", "c")
#' order = c("a", "b", "c", "d", "e", "f")
#' me_hodges_lehmann_os(ex3, levels=order)
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



