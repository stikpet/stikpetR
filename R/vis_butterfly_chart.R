# Suppress R CMD check note for ggplot variable bindings
utils::globalVariables(c("Var1", "Var2", "Freq"))

#' Butterfly Chart / Tornado Chart / Pyramid Chart
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_classic
#' 
#' @description
#' A special case of diverging bar charts when only comparing two categories. 
#' 
#' Depending on the ordering of the results different names exist. I've chosen to use 'butterfly' if no ordering is done, 'pyramid' if they are ordered from small to large, and 'tornado' when going from large to small.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/QeyqW5Vk69o) and the diagram is also discussed at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/PyramidChart.html
#' 
#' @param bin_field : dataframe field with categories for the rows
#' @param bin_field : dataframe field with categories for the columns
#' @param field_categories : optional list with selection of categories of field1
#' @param bin_categories : optional list with selection of categories of field2
#' @param variation : optional order of the bars. Either "butterfly" (default), "tornado", or "pyramid" 
#' @param roundHigh : optional to adjust number of tickmarks on horizontal axis
#' @param show : {'count', 'bin-percent', 'field-percent', 'overall-percent'}, optional show either counts or percentage
#' @param rotate : bool, optional rotate the bars so they appear horizontal. Default is True
#' @param xlbl : string, optional label for the field axis, if not set the name of the field is used.
#' @param colors : list, optional vector with two colors, one for each category to use
#' @param show_grid : bool, optional show a grid on the chart. Default is False
#' 
#' @returns
#' plot
#' 
#' @details
#' The term *butterfly chart* can for example be found in Hwang and Yoon (2021, p. 25).
#' 
#' The term *tornado diagrom* can be found in the guide from the Project Management Institute (2013, p. 338). The term *funnel chart* is also sometimes used (for example Jamsa (2020, p. 135)), but this is also a term sometimes used for a more analytical scatterplot used for some specific analysis.
#' 
#' The term *pyramid chart* can for example be found in Schwabish (2021, p. 185). It is very often used for comparing age distributions.
#' 
#' @references
#' Hwang, J., & Yoon, Y. (2021). Data analytics and visualization in quality analysis using Tableau. CRC Press.
#' 
#' Jamsa, K. (2020). Introduction to data mining and analytics: With machine learning in R and Python. Jones & Bartlett Learning.
#' 
#' Project Management Institute (Ed.). (2013). A guide to the project management body of knowledge (5th ed.). Project Management Institute, Inc.
#' 
#' Schwabish, J. (2021). Better data visualizations: A guide for scholars, researchers, and wonks. Columbia University Press.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' vi_butterfly_chart(df1[['mar1']], df1[['sex']], roundHigh=100)
#' vi_butterfly_chart(df1[['mar1']], df1[['sex']], variation="tornado", roundHigh=100)
#' vi_butterfly_chart(df1[['mar1']], df1[['sex']], variation="pyramid", roundHigh=100)
#' 
#' @export
vi_butterfly_chart <- function(field, 
                               bin_field,
                               field_categories=NULL, 
                               bin_categories=NULL, 
                               variation='butterfly', 
                               roundHigh=5, 
                               show='count', 
                               rotate=TRUE, 
                               xlbl=NULL,
                               colors = c('orange', 'blue'), 
                               show_grid=FALSE){
  
  
  percent <- switch(
    show,
    "count" = NULL,
    "bin-percent" = "column",
    "field-percent" = "row",
    "overall-percent" = "all",
    stop("Invalid value for 'show'")
  )
  
  ct<-tab_cross(field, bin_field, order1=field_categories, order2=bin_categories, percent=percent)  
  
  if (variation=="tornado") {
    ct <- addmargins(ct)
    ct <- ct[order(ct[,3]),]
    ct <- ct[0:(nrow(ct)-1), 0:2]
  }
  
  if (variation=="pyramid") {
    ct <- addmargins(ct)
    ct <- ct[order(ct[,3], decreasing=TRUE),]
    ct <- ct[2:(nrow(ct)), 0:2]
  }
  
  if(is.null(xlbl)){
    xlbl = deparse(substitute(field))}
  
  
  if (show=='count'){ylbl = 'frequency'}
  if (show=='bin-percent'){ylbl = paste0('percent of ', deparse(substitute(bin_field)))}
  if (show=='field-percent'){ylbl = paste0('percent of ', xlbl)}
  if (show=='overall-percent'){ylbl = 'percent of overall'}
  
  high <- roundHigh * (round(max(ct)/roundHigh, 0)+1)
  
  ctDf <- as.data.frame(ct)
  colnames(ctDf) = c("Var1", "Var2", "Freq")
  
  plot <- ggplot(ctDf, aes(x = Var1, fill = Var2, y = ifelse(test = Var2 == rownames(table(Var2))[1], yes = -Freq, no = Freq))) + 
    geom_bar(stat = "identity") + scale_fill_manual(values = colors) + 
    scale_y_continuous(limits = c(-high, high), breaks=seq(-high,high,roundHigh), labels=abs(seq(-high,high,roundHigh))) + 
    xlab(xlbl) + 
    ylab(ylbl) +
    guides(fill=guide_legend(title=xlbl)) 
  
  if (rotate){plot = plot + coord_flip()}
  if (!show_grid){plot = plot + theme_classic()}
  return(plot)
}



