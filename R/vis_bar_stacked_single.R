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
#' @param catCoding optional vector with the order for the bars
#' @param orientation optional to indicate horizontal or vertical chart Either `"h"` (default) or `"v"`
#' @return The chart.
#' 
#' @details 
#' This function basically uses barplot(...,beside = FALSE) from R's *graphics* library
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
vi_bar_stacked_single <- function(data, catCoding=NULL, orientation=c("h", "v")){
  
  if (length(orientation)>1){orientation="h"}
  if (orientation=="h"){horizontal = TRUE}
  else {horizontal = FALSE}
  
  varname=deparse(substitute(data)) 
  if (!is.null(catCoding)){
    legendLabels = catCoding
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = catCoding)
    data = as.numeric(myFieldOrd)
  }
  #determine the counts (frequencies)
  freqs = table(data)
  
  #determine the number of categories (k)
  k = length(freqs)
  
  if (is.null(catCoding)){legendLabels=rownames(freqs)}
  
  freqs = freqs/sum(freqs)*100
  
  chart = barplot(as.matrix(freqs),
                  beside = FALSE,
                  legend.text = legendLabels,
                  horiz = horizontal,
                  col = heat.colors(k),
                  xlab = "percent",
                  ylab = varname)
  
  return(chart)
}