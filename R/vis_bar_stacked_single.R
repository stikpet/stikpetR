#' Single Stacked Bar-Chart
#' 
#' @description 
#' A regular bar-chart but with the bars on top of each other, instead of next to each other. 
#' This is called a compound bar chart, stacked bar chart (Wilkinson, 2005, p. 157) or 
#' component bar chart (Zedeck, 2014, p. 54). 
#' 
#' It can be defined as: "a bar chart showing multiple bars stacked at each x-axis category, 
#' each representing a value of the stacking variable" (Upton & Cook, 2014, p. 88).
#' 
#' This function is shown in this [YouTube video](https://youtu.be/j92bv5gFwpI) and the visualisation is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/bar-chart.html)
#' 
#' @param data the data from which to create the bar-chart
#' @param catCoding optional vector with the order for the bars
#' @param orientation optional to indicate horizontal or vertical chart Either `"h"` (default) or `"v"`
#' 
#' @return The chart.
#' 
#' @details 
#' This function basically uses barplot(...,beside = FALSE) from R's *graphics* library
#' 
#' @section Before, After and Alternatives:
#' Before the visualisation you might first want to get an impression using a frequency table:
#' \code{\link{tab_frequency}}, for a frequency table
#' 
#' After visualisation you might want some descriptive measures:
#' \code{\link{me_consensus}}, for the Consensus. 
#' \code{\link{me_hodges_lehmann_os}}, for the Hodges-Lehmann Estimate (One-Sample).
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
#' Alternatives for this visualisation could be:
#' \code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.
#' 
#' @references 
#' Upton, G. J. G., & Cook, I. (2014). *Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Wilkinson, L. (2005). *The grammar of graphics* (2nd ed). Springer.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @importFrom grDevices heat.colors
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Text dataframe
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' vi_bar_stacked_single(ex1, catCoding=order)
#' vi_bar_stacked_single(ex1, catCoding=order, orientation="v");
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' vi_bar_stacked_single(ex2);
#' 
#' @export
vi_bar_stacked_single <- function(data, catCoding=NULL, orientation=c("h", "v")){
  
  if (length(orientation)>1){orientation="h"}
  if (orientation=="h"){horizontal = TRUE}
  else {horizontal = FALSE}
  
  varname=deparse(substitute(data)) 
  if (!is.null(catCoding)){
    legendLabels = catCoding
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = catCoding)
    data = as.numeric(myFieldOrd)
  }
  #determine the counts (frequencies)
  freqs = table(data)
  
  #determine the number of categories (k)
  k = length(freqs)
  
  if (is.null(catCoding)){legendLabels=rownames(freqs)}
  
  freqs = freqs/sum(freqs)*100
  
  chart = barplot(as.matrix(freqs),
                  beside = FALSE,
                  legend.text = legendLabels,
                  horiz = horizontal,
                  col = heat.colors(k),
                  xlab = "percent",
                  ylab = varname)
  
  return(chart)
}



