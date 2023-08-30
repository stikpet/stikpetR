#' Split Histogram
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' 
#' @description 
#' Based on a categorical field the scores for each category are plotted in a separate histogram and each of the histograms is placed underneath each other.
#' 
#' See **vi_histogram()** for more details on histograms.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with categories to use
#' @param ... other parameters for use in geom_histogram function
#' 
#' @returns
#' The split histogram
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_histogram_split <- function(catField, scaleField, categories=NULL, ...){
  #get the original names of the variables
  argnames = sys.call()
  argnames = unlist(lapply(argnames[-1], as.character))
  
  #create a dataframe
  data = data.frame(catField, scaleField)
  data = na.omit(data)
  
  if (!is.null(categories)){
    data = data[data[, 1] %in% categories, ]
  }
  
  #create the split histogram
  colnames(data) = c("category", "score")
  ggplot(data, aes(x=score)) + 
    geom_histogram(...) + 
    facet_wrap(~ category, ncol=1) + labs(x = argnames[2])
  
}