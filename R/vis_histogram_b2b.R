# Suppress R CMD check note for ggplot variable bindings
utils::globalVariables(c("score"))

#' Back-to-Back Histogram
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 aes
#' 
#' @description 
#' This function creates a simple back-to-back histogram. This is sometimes also referred to as a Pyramid chart or a Dual-Sided histogram (Jelen, 2005).
#' 
#' The back-to-back histogram together with back-to-back stem-and-leaf and split box-plots are described as effective ways to compare two distributions by Lane and Sándor (2009).
#' 
#' The visualisation is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/histogram.html)
#' 
#' See **vi_histogram()** for more details on histograms.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with the two categories to use from bin_field. If not set the first two found will be used
#' @param bins optional list with upper bounds of bins to be used, or any of the pre-set options from hist() 
#' @param equal_bins optional boolean. use the same bins for each sample
#' @param density optional boolean. use frequency density instead of counts
#' @param ... other parameters for use in geom_histogram function
#' 
#' @returns
#' The back-to-back histogram
#' 
#' @references 
#' Jelen, B. (2005, December 24). Dual Sided Histogram in Excel. MrExcel. https://www.mrexcel.com/tech-tv/dual-sided-histogram-in-excel/
#' 
#' Lane, D. M., & Sándor, A. (2009). Designing better graphs by including distributional information and integrating words, numbers, and images. *Psychological Methods, 14*(3), 239–257. doi:10.1037/a0016620
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_histogram_b2b <- function(catField, scaleField, categories=NULL, bins=NULL, equal_bins=TRUE, density=FALSE, ...){
  # DATA PREPARATION
  # create dataframe and remove missing values
  df <- na.omit(data.frame(score=scaleField, category=catField))
  
  #the two categories
  if (!is.null(categories)){
    #use the provided categories
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    # use first two categories found
    cat1 = names(table(df[ ,2]))[1]
    cat2 = names(table(df[ ,2]))[2]
  }
  
  # DETERMINE BREAKS
  #split the scores across the categories    
  scoresCat1 = unname(unlist((subset(df, df[ ,2] == cat1)[1])))
  scoresCat2 = unname(unlist((subset(df, df[ ,2] == cat2)[1])))
  #combine this into one long list
  allScores = c(scoresCat1, scoresCat2)
  #if bins not set, use Sturges
  if (is.null(bins)){
    bins = "Sturges"}
  #use hist to get the breaks
  if (equal_bins){        
    h <- hist(allScores, breaks = bins, plot = FALSE)
    bins1 <- h$breaks
    bins2 <- bins1}
  else {
    h1 <- hist(scoresCat1, breaks = bins, plot = FALSE)
    bins1 <- h1$breaks
    h2 <- hist(scoresCat2, breaks = bins, plot = FALSE)
    bins2 <- h2$breaks
  }
  
  # CREATE THE PLOT
  if (density){
    hist1 = geom_histogram(
      data = df[df$category == cat1, ],
      aes(x = score, y = after_stat(density * sum(count)), fill = category),
      col = "black",
      breaks = bins1, 
      ...)
    
    hist2 = geom_histogram(
      data = df[df$category == cat2, ],
      aes(x = score, y = -after_stat(density * sum(count)), fill = category),
      col = "black",
      breaks = bins2, 
      ...)    
    
    ggplot() + hist1 + hist2 + scale_y_continuous(labels = abs) + labs(x = "score", y = "Density") + coord_flip()
  }
  else {
    hist1 = geom_histogram(
      data = df[df$category == cat1, ],
      aes(x = score, y = after_stat(count), fill = category),
      col = "black",
      breaks = bins1, 
      ...)
    
    hist2 = geom_histogram(
      data = df[df$category == cat2, ],
      aes(x = score, y = -after_stat(count), fill = category),
      col = "black",
      breaks = bins2, 
      ...)    
    ggplot() + hist1 + hist2 + 
      scale_y_continuous(labels = abs) + 
      labs(x = "score", y = "Frequency") + 
      coord_flip()
  }    
}