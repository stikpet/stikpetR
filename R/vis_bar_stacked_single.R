#' Single Stacked Bar-Chart
#' 
#' @description 
#' A regular bar-chart but with the bars on top of each other, instead of next to each other. 
#' This is called a compound bar chart, stacked bar chart (Wilkinson, 2005, p. 157) or 
#' component bar chart (Zedeck, 2014, p. 54). 
#' 
#' It can be defined as: “a bar chart showing multiple bars stacked at each x-axis category, 
#' each representing a value of the stacking variable” (Upton & Cook, 2014, p. 88).
#' 
#' @param data the data from which to create the bar-chart
#' @param varname name of the variable, if not provided name of data is used.
#' @param height which values to show on axis, either percent or count
#' @return The chart.
#' 
#' @details 
#' This function basically uses barplot(...,beside = FALSE) from R's *graphics* library
#' 
#' @examples 
#' ordData <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' vi_bar_stacked_single(ordData)
#' vi_bar_stacked_single(ordData, height="count")
#' vi_bar_stacked_single(ordData, varname="scores", height="percent")
#' 
#' @seealso 
#' An alternative chart for a single ordinal variable could be a dual axis bar chart, see \code{\link{vi_bar_dual_axis}} 
#' 
#' @references 
#' Upton, G. J. G., & Cook, I. (2014). *Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Wilkinson, L. (2005). *The grammar of graphics* (2nd ed). Springer.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
vi_bar_stacked_single <- function(data, varname=NULL, height=c("percent", "count")){
  
  #set default if not provided
  if (length(height)>1) {
    height="percent"
  }
  
  #determine the counts (frequencies)
  freqs = table(data)
  
  #determine the number of categories (k)
  k = length(freqs)
  
  #set variable name to data name if not provided
  if (is.null(varname)) {
    varname=deparse(substitute(data))
  }
  
  if (height=="percent") {
    freqs = freqs/sum(freqs)*100
  }
  chart = barplot(as.matrix(freqs),
                  beside = FALSE,
                  legend.text = rownames(freqs),
                  horiz = TRUE,
                  col = heat.colors(k),
                  xlab = height,
                  ylab = varname)
  
  return(chart)
}