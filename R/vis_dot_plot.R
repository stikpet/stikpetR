#' Dot Plot
#' 
#' @param data the data from which to create the dot plot
#' @param stackRatio ratio on how close the dots are to each other
#' @param dotSize indicator for how big the dots need to be
#' @return chart the dot plot
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_dotplot
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 aes
#' 
#' @description 
#' The Oxford Dictionary of Statistics defines a dot plot as "an alternative to a bar chart 
#' or line graph when there are very few data values. Each value is recorded as a dot, 
#' so that the frequencies for each value can easily be counted" (Upton & Cook, 2014, p. 129). 
#' 
#' This function uses ggplot2 geom_dotplot() to create a simple dot plot.
#' 
#' A [YouTube](https://youtu.be/qs1nh0CMiIY) video on pie charts.
#' 
#' @details 
#' In the definition a *bar chart* is mentioned. A bar chart can be defined as 
#' “a graph in which bars of varying height with spaces between them are used to 
#' display data for variables defined by qualities or categories” (Zedeck, 2014, p. 20). 
#' Together this indicates that a dot plot is used for categorical data.
#' 
#' However, Zedeck sees the dot plot as an alternative name for a scatterplot, 
#' which is for continuous data. A third version comes from the 
#' Cambridge Dictionary of Statistics: "A more effective display than a number of 
#' other methods, for example, pie charts and bar charts, for displaying quantitative 
#' data which are labelled" (Everitt, 2004, p. 123). They also show an example where 
#' we see a categorical variable on one axis, and a continuous variable on another.
#' 
#' This function was only for the original definition for categorical data.
#' 
#' @examples 
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_dot_plot(data, stackRatio=1.5, dotSize=2)
#' 
#' @references 
#' Everitt, B. (2004). *The Cambridge dictionary of statistics* (2nd ed). Cambridge University Press.
#' 
#' Upton, G. J. G., & Cook, I. (2014). *Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
vi_dot_plot <- function(data, stackRatio=1, dotSize=1){
  
  dFr = as.data.frame(data)
  
  chart= ggplot2::ggplot(dFr, aes(x = data)) + ggplot2::geom_dotplot(stackratio = stackRatio, dotsize=dotSize) + ggplot2::scale_y_continuous(NULL, breaks = NULL)
  
  return(chart)  
}