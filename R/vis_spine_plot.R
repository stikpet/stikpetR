#' Spine Plot
#' 
#' @description 
#' A spine plot is similar to a multiple stacked bar-chart, but "the difference is that the bars fill the plot vertically so the shading gives us proportions instead of counts. Also, the width of each bar varies, reflecting the marginal proportion of observations in each workshop" (Muenchen, 2006, p. 286)
#' 
#' It is a chart you could use when with two nominal variables and do not have a clear independent and dependent variable. Otherwise a multiple/clustered bar-chart might be preferred.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection of categories of field1
#' @param categories2 : optional list with selection of categories of field2
#' 
#' @returns
#' spine plot
#' 
#' @details
#' The naming of this diagram is unfortunately not very clear. I use the term 'spine plot' as a special case of a Mosaic Plot. Mosaic Plots are often attributed to Hartigan and Kleiner (for example by Friendly (2002, p. 90)). Earlier versions are actually known, for example Walker (1874, p. PI XX). Hartigan and Kleiner (1981) start their paper with a Mosaic Plot for a cross table, but end it with showing Mosaic Plots for multiple dimension cross tables.
#' 
#' A Marimekko Chart is simply an alternative name for the Mosaic Plot, although according to Wikipedia "mosaic plots can be colored and shaded according to deviations from independence, whereas Marimekko charts are colored according to the category levels" (Wikipedia, 2022).
#' 
#' The term 'Spine Plot' itself is often attributed to Hummel, but I've been unable to hunt down his original article: Linked bar charts: Analysing categorical data graphically. Computational Statistics 11: 23-33.
#' 
#' @references
#' Carvalho, T. (2021, April 10). Marimekko Charts with Python's Matplotlib. Medium. https://towardsdatascience.com/marimekko-charts-with-pythons-matplotlib-6b9784ae73a1
#' 
#' Friendly, M. (2002). A brief history of the mosaic display. *Journal of Computational and Graphical Statistics, 11*(1), 89-107. https://doi.org/10.1198/106186002317375631
#' 
#' Hartigan, J. A., & Kleiner, B. (1981). Mosaics for contingency tables. In W. F. Eddy (Ed.), Proceedings of the 13th Symposium on the Interface (pp. 268-273). Springer. https://doi.org/10.1007/978-1-4613-9464-8_37
#' 
#' Muenchen, R. A. (2009). *R for SAS and SPSS Users*. Springer.
#' 
#' Walker, F. A. (1874). *Statistical atlas of the United States based on the results of the ninth census 1870*. Census Office.
#' 
#' Wikipedia. (2022). Mosaic plot. In Wikipedia. https://en.wikipedia.org/w/index.php?title=Mosaic_plot&oldid=1089465331
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' vi_bar_clustered(df1[['mar1']], df1[['sex']], percent="column")
#' 
#' @export
vi_spine_plot <- function(field1, 
                          field2,
                          categories1=NULL, 
                          categories2=NULL){
  
  field1Name <- deparse(substitute(field1))
  field2Name <- deparse(substitute(field2))
  ct = tab_cross(field2, field1, order1=categories2, order2=categories1)
  
  spineplot(ct, xlab=field1Name, ylab=field2Name)
}