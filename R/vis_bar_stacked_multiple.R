#' Multiple Stacked Bar-Chart
#' 
#' @description
#' To visualise an ordinal variable, it often makes sense to stack the results. Stacking the results creates a compound bar chart, or sometimes stacked bar chart (Wilkinson, 2005, p. 157) or component bar chart (Zedeck, 2014, p. 54). It can be defined as: “a bar chart showing multiple bars stacked at each x-axis category, each representing a value of the stacking variable” (Upton & Cook, 2014, p. 88).
#' 
#' Instead of one bar (see **vi_bar_stacked_single()**), we can create two or more (one for each group). This could then be considered a multiple compound bar-chart.
#' 
#' @param catField list or dataframe with the categories
#' @param ordField list or dataframe with the scores
#' @param levels optional list with the scores in order
#' @param ... optional, other parameters for use in barplot function
#' 
#' @returns 
#' multiple stacked bar-chart
#' 
#' @details
#' This function is more like a wrapper for the **barplot()** from R *graphics* library.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @references
#' Upton, G., & Cook, I. (2014). *Oxford: Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Wilkinson, L. (2005). *The grammar of graphics* (2nd ed). Springer.
#' 
#' @examples 
#' file1 = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(file1, sep=",", na.strings=c("", "NA"))
#' vi_bar_stacked_multiple(df1[['mar1']], df1[['accntsci']], ylab= "percent", col=1:5)
#' 
#' cats = c(1, 1, 2, 2, 2, 3, 3, 3, 3)
#' scor = c(1, 2, 1, 1, 2, 1, 1, 1, 2)
#' vi_bar_stacked_multiple(cats, scor, ylab= "percent", col=1:5)
#' 
#' 
#' @export
vi_bar_stacked_multiple <- function(catField, ordField, levels=NULL, ...){
  #replace the ordinal values if levels is provided
  if (!is.null(levels)){
    ordField = factor(ordField, ordered = TRUE, levels = levels)        
  }
  
  ct = table(ordField, catField)
  pct = prop.table(ct, 2)*100
  
  barplot(pct, legend = rownames(pct), ...)
  
}