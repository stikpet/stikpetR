#' Interquartile Range, Semi-Interquartile Range and Mid-Quartile Range
#' 
#' @description 
#' There are some measures of dispersion that instead of using the full range (i.e. maximum minus minimum), make use of the quartiles. The advantage of this, is that it is less influenced by extreme values.
#' 
#' The Interquartile Range (Galton, 1881, p. 245) is the range how big the difference is between the third and the first quartile. If Tukey's method for the quartiles is used (*method="tukey"*), referred to as hinges, this is then also known as H-spread (Tukey, 1977, p. 44)
#' 
#' Yule (1911, p. 147) used half the inter-quartile range and labelled this Semi-Interquartile Range which he preferred over the term Quartile Deviation..
#' 
#' There is also a measure of central tendency that uses the quartiles, the Mid-Quartile (Parzen, 1980, p. 19), which is the average of the first and second quartile. It is also sometimes referred to as the Mid-Quartile Range (see for example Luo et al. (2018, p. 2), who refer to Triola, but Triola doesn't add the 'range' (Triola, 2010, p. 120))
#' 
#' The function uses the *me_quartiles* function and any of the methods from that function can be used.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/mUl3LFzfsfg) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Measures/QuartileRanges.html)
#' 
#' @param data vector or dataframe with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param measure the specific measure to determine. Either `"iqr"` (default), `"siqr"`, `"qd"`, or `"mqr"`
#' @param method the method to use to determine the quartiles
#' 
#' @returns
#' A dataframe with:
#' \item{Q1}{the first (lower) quartile}
#' \item{Q3}{the third (upper/higher) quartile}
#' \item{range}{the range determined}
#' 
#' @details 
#' The formula used for the Interquartile Range is:
#' \deqn{IQR = Q_3 - Q_1}
#' This can be obtained by setting *range="iqr"*.
#' 
#' The IQR is mentioned in Galton (1881, p. 245) and the H-spread in Tukey (1977, p. 44).
#' 
#' The H-spread can be obtained by setting *range="iqr"* and *method="tukey"*.
#' 
#' The formula used for the Semi-Interquartile Range (Quartile Deviation) is (Yule, 1911, p. 147):
#' \deqn{SIQR = \frac{Q_3 - Q_1}{2}}
#' This can be obtained by setting *range="siqr"* or *range="qd"*.
#' 
#' The formula for the mid-quartile range used is:
#' \deqn{MQR = \frac{Q_3 + Q_1}{2}}
#' This can be obtained by setting *range="mqr"*.
#' 
#' This formula can be found in Parzen (1980, p. 19), but there are probably older references.
#' 
#' @section Before, After and Alternatives:
#' Before this measure you might want an impression using a frequency table or a visualisation:
#' \code{\link{tab_frequency}}, for a frequency table
#' \code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
#' \code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.
#' 
#' After this you might want some other descriptive measures:
#' \code{\link{me_consensus}}, for the Consensus. 
#' \code{\link{me_hodges_lehmann_os}}, for the Hodges-Lehmann Estimate (One-Sample).
#' \code{\link{me_median}}, for the Median.
#' \code{\link{me_quantiles}}, for Quantiles.
#' \code{\link{me_quartiles}}, for Quartiles / Hinges.
#' 
#' or perform a test:
#' \code{\link{ts_sign_os}}, for One-Sample Sign Test.
#' \code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
#' \code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
#' 
#' 
#' @references 
#' Galton, F. (1881). Report of the anthropometric committee. *Report of the British Association for the Advancement of Science, 51*, 225-272.
#' 
#' Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. *Statistical Methods in Medical Research, 27*(6), 1785-1805. doi:10.1177/0962280216669183
#' 
#' Parzen, E. (1980). *Data modeling using quantile and density-quantile functions*. Institute of Statistics, Texas A&M University.
#' 
#' Triola, M. F. (2010). *Elementary statistics* (11th ed). Addison-Wesley.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' Yule, G. U. (1911). *An introduction to the theory of statistics*. Charles Griffin.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Text dataframe
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' me_quartile_range(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' me_quartile_range(ex2)
#' 
#' @export
me_quartile_range <- function(data, levels=NULL,
                              measure=c("iqr", "siqr", "qd", "mqr"), 
                              method="cdf"){
  
  if (length(measure)>1){measure="iqr"}
  
  if (is.null(levels)){
    dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  Qs = me_quartiles(dataN, method=method)
  q1 = Qs$q1
  q3 = Qs$q3
  
  if (measure=="iqr") {
    r = q3 - q1
    if (method=="tukey" || method=="inclusive" || method=="tukey" || method=="vining" || method=="hinges") {rName = "Hspread"}
    else{rName = "IQR"}
  }
  else if (measure=="siqr" || measure=="qd"){
    r = (q3 - q1)/2
    rName = "SIQR"
  }
  else if (measure=="mqr"){
    r = (q3 + q1)/2
    rName = "MQR"
  }
  
  results = data.frame(q1, q3, r)
  colnames(results) = c("Q1", "Q3", rName)
  
  
  return (results)
  
}



