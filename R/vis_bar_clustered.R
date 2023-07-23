#' Clustered / Multiple Bar Chart
#' 
#' @description
#' A bar-chart is defined as “a graph in which bars of varying height with spaces between them are used to display data for variables defined by qualities or categories” (Zedeck, 2014, p. 20).
#' 
#' The bars can be split into multiple bars based on another variable. This is then known as a multiple bar-chart (Kemp, 2004, p. 150) or clustered bar-chart (Brase, 2009, p. 50; Griffith, 2007, p. 168).
#' 
#' It can be defined as “a bar chart for comparing the frequencies of a categorical variable in two or more situations” (Upton & Cook, 2014, p. 283). 
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param order1 : optional list with order for categories of field1
#' @param order2 : optional list with order for categories of field2
#' @param percent : optional which percentages to show. Either "none" (default), "all", "row", "column"
#' 
#' @returns
#' clustered bar chart
#' 
#' @references
#' Brase, C. (2009). *Understandable statistics* (9th ed.). Houghton MIfflin.
#' 
#' Griffith, A. (2007). *SPSS for dummies*. Wiley.
#' 
#' Kemp, S. M., & Kemp, S. (2004). *Business statistics demystified*. McGraw-Hill.
#' 
#' Upton, G., & Cook, I. (2014). *Oxford: Dictionary of statistics* (3rd ed.). Oxford University Press.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example 1: Clustered Bar Chart in percentages
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' vi_bar_clustered(df1[['mar1']], df1[['sex']], percent="column")
#' 
#' #Example 2: Specified order
#' orderR = c("DIVORCED", "WIDOWED", "SEPARATED", "MARRIED", "NEVER MARRIED")
#' orderC = c("MALE", "FEMALE")
#' vi_bar_clustered(df1[['mar1']], df1[['sex']], order1=orderR, order2=orderC)
#' 
#' @export
vi_bar_clustered <- function(field1, 
                             field2, 
                             order1=NULL, 
                             order2=NULL, 
                             percent=c(NULL, "all", "row", "column")){
  
  tab = tab_cross(field2, field1, order1=order2, order2=order1, percent=percent, totals="exclude")
  
  barplot(tab, beside=TRUE,legend = rownames(tab))
  
}