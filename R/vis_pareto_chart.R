#' Pareto Chart
#' 
#' @description 
#' The Pareto Chart gets its name from the Pareto Principle, which is named after 
#' Vilfredo Pareto. This principle states that roughly 80% of consequencies come from 
#' 20% of causes (Pareto, 1896).
#' 
#' Unfortunately, there is no general agreed upon definition of a Pareto diagram. 
#' The most general description I’ve found was by Kemp and Kemp (2004) who mention it 
#' is a name for a bar chart if the order of the bars have no meaning 
#' (i.e. for a nominal variable), and they only mention that often the bars are then 
#' placed in decreasing order. According to some authors a Pareto diagram is any diagram 
#' with the bars in order of size (Joiner, 1995; WhatIs.com, n.d.), while others suggest 
#' that a line representing the cumulative relative frequencies should also be 
#' included (Weisstein, 2002). Upton and Cook (2014) also add that the bars should not 
#' have any gaps, but many other authors ignore this.
#' 
#' The following definition by the author is used: a bar chart where the bars are placed 
#' in descending order of frequency. Usually an ogive is added in the chart as well.
#' 
#' An ogive (oh-jive) is: "the graphs of cumulative frequencies" (Kenney, 1939).
#' 
#' A video on Pareto charts is available [here](https://youtu.be/kDp5zPfK-Po).
#' 
#' @param data the data from which to create a Pareto chart
#' @param varname a name for the data, if not provided the name of the data variable is used
#' 
#' @return a Pareto chart in the plot window
#' 
#' @references 
#' Joiner. (1995). Pareto charts: Plain & simple. Joiner Associates.
#' 
#' Kemp, S. M., & Kemp, S. (2004). *Business statistics demystified*. McGraw-Hill.
#' 
#' Kenney, J. F. (1939). *Mathematics of statistics; Part one*. Chapman & Hall.
#' 
#' Pareto, V. (1896). *Cours d’économie politique* (Vol. 1). Lausanne.
#' 
#' Upton, G. J. G., & Cook, I. (2014). *Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Weisstein, E. W. (2002). *CRC concise encyclopedia of mathematics* (2nd ed.). Chapman & Hall/CRC.
#' 
#' WhatIs.com. (n.d.). What is Pareto chart (Pareto distribution diagram)? - Definition from WhatIs.com. Retrieved April 20, 2014, from http://whatis.techtarget.com/definition/Pareto-chart-Pareto-distribution-diagram
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_pareto_chart(ex1);
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_pareto_chart(ex2);
#' 
#' @export
vi_pareto_chart <- function(data, varname=NULL){
  
  #set variable name to data name if not provided
  if (is.null(varname)) {
    varname=deparse(substitute(data))
  }
  
  freq = table(data)
  k = length(freq)
  freqSort = sort(freq, decreasing = TRUE)
  cumFr = cumsum(freqSort)
  cumPerc = cumFr /sum(freq)
  op = par(mar= c(5.1,4.1,4.1,4.1))
  barplot(freqSort)
  par(new=TRUE)
  plot(cumPerc, type = 'b', xlim=c(0.5, k+0.5), ylim=c(0,1), col = "red", axes = FALSE, xlab = varname, ylab = "count")
  mtext("cumulative percent", side = 4, line = 3)
  axis(side = 4)
  par(op)
  
}