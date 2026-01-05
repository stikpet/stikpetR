#' Overlaid Histogram
#' 
#' @description 
#' This function creates a simple overlaid (overlapping) histogram. 
#' 
#' See **vi_histogram()** for more details on histograms.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with the two categories to use from bin_field. If not set the first two found will be used
#' @param bins optional list with upper bounds of bins to be used, or any of the pre-set options from hist() 
#' @param show c('count', 'relative'), show either counts or relative count on vertical scale
#' @param density c(FALSE, TRUE, 'auto'), show (relative) frequency or use (relative) frequency density.
#' @param title string, title on top of chart, default is 'histogram'
#' @param colors vector, two colors, one for each category, default is c(rgb(0, 0, 1, 1/4), rgb(1, 0, 0, 1/4))
#' 
#' @returns
#' The overlaid histograms
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_histogram_overlay <- function(catField, scaleField, categories=NULL, bins=NULL, show='count', density=FALSE, title='histogram', colors=c(rgb(0, 0, 1, 1/4), rgb(1, 0, 0, 1/4))){
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
  
  h <- hist(allScores, breaks = bins, plot = FALSE)
  bins <- h$breaks
  
  xLabel = deparse(substitute(scaleField))
  
  h1 <- hist(scoresCat1, breaks=bins, plot = FALSE)
  h2 <- hist(scoresCat2, breaks=bins, plot = FALSE)
  
  if (density=='auto'){
    if (h$equidist){density=FALSE}
    else {density=TRUE}}
  
  if (show=="count" && density==FALSE){
    yLabel = "Frequency"
    histFreq = TRUE
    # set densities to counts, in case unequal bin-widths
    h1$density = h1$counts
    h2$density = h2$counts
    ylimits = c(0, max(max(h1$counts), max(h2$counts)))
  }
  
  else if (show=="relative" && density==FALSE){                
    h1$counts <- h1$counts / sum(h1$counts) * 100
    h2$counts <- h2$counts / sum(h2$counts) * 100
    ylimits = c(0, max(max(h1$counts), max(h2$counts)))
    # set densities to counts, in case unequal bin-widths
    h1$density = h1$counts
    h2$density = h2$counts
    yLabel = "Percent of category"
    histFreq = TRUE        
  }
  
  else if (show=="count" && density==TRUE){
    h1$density <- h1$density * length(scoresCat1)
    h2$density <- h2$density * length(scoresCat2)
    ylimits = c(0, max(max(h1$density), max(h2$density)))
    yLabel = "Frequency density"
    histFreq = FALSE        
  } 
  
  else if (show=="relative" && density==TRUE){
    ylimits = c(0, max(max(h1$density), max(h2$density)))
    yLabel = "Relative Frequency density"
    histFreq = FALSE
  } 
  
  suppressWarnings(plot(h1, freq = histFreq, main=title, xlab = xLabel, ylab = yLabel, col = colors[1], ylim=ylimits))
  
  suppressWarnings(plot(h2, freq = histFreq, col = colors[2], add=TRUE))
  
  legend("topright",
         legend = c(cat1, cat2),
         fill = c(colors[1], colors[2]),
         border = c("black", "black"),
         bty = "n")
}