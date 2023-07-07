#' Dual-Axis Bar Chart
#' 
#' @description 
#' A dual axis bar-chart is a bar-chart with two vertical axis. In this function it will 
#' show both the count and cumulative proportion. 
#' 
#' This chart could be used with a single ordinal variable.
#' 
#' @param data the data from which to create a Pareto chart
#' @param varname a name for the data, if not provided the name of the data variable is used
#' 
#' @return a chart in the plot window
#' 
#' @section Alternatives:
#' It is also possible to do this with the library *ggplot2*. 
#' A video on how to use ggplot2 for this can be found [here](https://youtu.be/EhK9K9jwsco)
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