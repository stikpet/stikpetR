#' Dual-Axis Bar Chart
#' 
#' @param data the data from which to create a Pareto chart
#' @param varname a name for the data, if not provided the name of the data variable is used
#' @return a chart in the plot window
#' 
#' @description 
#' A dual axis bar-chart is a bar-chart with two vertical axis. In this function it will 
#' show both the count and cumulative proportion. 
#' 
#' This chart could be used with a single ordinal variable.
#' 
#' @examples 
#' ordData <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' vi_bar_stacked_single(ordData)
#' 
#' @seealso 
#' An alternative chart for a single ordinal variable could be a single stacked bar-chart, 
#' see \code{\link{vi_bar_stacked_single}} 
#' 
#' @section Alternatives:
#' It is also possible to do this with the library *ggplot2*. 
#' A video on how to use ggplot2 for this can be found [here](https://youtu.be/EhK9K9jwsco)
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
vi_bar_dual_axis <- function(data, varname=NULL){
  #set variable name to data name if not provided
  if (is.null(varname)) {
    varname=deparse(substitute(data))
  }
  
  freqs<-table(ordData)
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