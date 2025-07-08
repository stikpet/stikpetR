#' Dual-Axis Bar Chart
#' 
#' @description 
#' A dual axis bar-chart is a bar-chart with two vertical axis. In this function it will 
#' show both the count and cumulative proportion. 
#' 
#' This chart could be used with a single ordinal variable.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/J62B1zyTO3U) and the visualisation is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/bar-chart.html)
#' 
#' @param data the data from which to create a Pareto chart
#' @param varname a name for the data, if not provided the name of the data variable is used
#' 
#' @return a chart in the plot window
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
#' \code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_bar_dual_axis(ex1);
#' vi_bar_dual_axis(ex1, varname="marital status");
#' 
#' @export
vi_bar_dual_axis <- function(data, varname=NULL){
  #set variable name to data name if not provided
  if (is.null(varname)) {
    varname=deparse(substitute(data))
  }
  
  freqs<-table(data)
  k = length(freqs)
  cumFr <- cumsum(freqs)
  cumPerc <- cumFr /sum(freqs)
  
  op <- par(mar= c(5.1,4.1,4.1,4.1))
  barplot(freqs)
  par(new=TRUE)
  plot(cumPerc, type = 'b', 
       xlim=c(0.5, k +0.5), ylim=c(0,1), 
       col = "red", 
       axes = FALSE, 
       xlab = varname, ylab = "count")
  mtext("cumulative proportion", side = 4, line = 3)
  axis(side = 4)
  par(op)
  
}