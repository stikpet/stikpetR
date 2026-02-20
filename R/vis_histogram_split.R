#' Split Histogram
#' 
#' @description 
#' Based on a categorical field the scores for each category are plotted in a separate histogram and each of the histograms is placed underneath each other.
#' 
#' See **vi_histogram()** for more details on histograms.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with the two categories to use from bin_field. If not set the first two found will be used
#' @param bins optional list with upper bounds of bins to be used, or any of the pre-set options from hist() 
#' @param show c('count', 'relative'), show either counts or relative count on vertical scale
#' @param density c('auto', FALSE, TRUE), show (relative) frequency or use (relative) frequency density.
#' @param shareY boolean, use same vertical axis range for each category
#' @param xlbl string, title for the horizontal axis. If set to NULL the name of the catField variable is used.
#' @param ... optional additional parameters to pass to plot()
#' 
#' @returns
#' The split histogram
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_histogram_split <- function(catField, scaleField, categories=NULL, bins=NULL, show='count', 
                               density='auto', shareY=TRUE, xlbl=NULL, ...){
  # DATA PREPARATION
  # create dataframe and remove missing values
  df <- na.omit(data.frame(score=scaleField, category=catField))
  #replace categories if provided
  if (!is.null(categories)){
    df = df[(df$category %in% categories),]}
  
  cats = unique(df[,2])
  k = length(cats)
  
  #if bins not set, use Sturges
  if (is.null(bins)){
    bins = "Sturges"}
  
  h <- hist(df[ ,1], breaks = bins, right=FALSE, plot = FALSE)
  bins <- h$breaks
  
  if (density=='auto'){
    if (h$equidist){density=FALSE}
    else {density=TRUE}}
  
  if (shareY){
    y_limit = 0
    for (i in 1:k){
      scores_i <- df[df$category == cats[i], 1]
      h1 <- hist(scores_i, breaks=bins, right=FALSE, plot = FALSE)
      if (show=="count" && density==TRUE){
        if (max(h1$density* length(scores_i)) > y_limit){ y_limit = max(h1$density* length(scores_i))}
      } 
      else if (show=="relative" && density==TRUE){
        if (max(h1$density) > y_limit){ y_limit = max(h1$density)}
      }
      else if (show=="relative" && density==FALSE){
        if (max(h1$counts / sum(h1$counts) * 100) > y_limit){ y_limit = max(h1$counts / sum(h1$counts) * 100)}
      }            
      else {
        if (max(h1$counts) > y_limit){ y_limit = max(h1$counts)}
      }
    }
  }
  
  if (is.null(xlbl)){
    xLabel = deparse(substitute(scaleField))}
  else {xLabel = xlbl}
  
  # split plotting area
  par(mfrow = c(k, 1))
  
  for (i in 1:k){
    if (i!=k){par(mar = c(2, 4, 2, 1))}
    else {par(mar = c(4, 4, 2, 1))}
    scores_i <- df[df$category == cats[i], 1]
    h1 <- hist(scores_i, breaks=bins, right=FALSE, plot = FALSE)
    show_freq = TRUE
    y_label = 'frequency'
    if (show=="count" && density==TRUE){
      h1$density <- h1$density * length(scores_i)
      show_freq = FALSE
      y_label = 'frequency density'}
    else if (show=="relative" && density==TRUE){
      show_freq = FALSE
      y_label = 'relative frequency density'}
    else if (show=="relative" && density==FALSE){
      h1$counts <- h1$counts / sum(h1$counts) * 100
      show_freq = TRUE
      y_label = 'percent'}
    
    
    plot(h1, freq=show_freq, main=cats[i], xlab=xLabel, ylab=y_label, ylim=c(0, y_limit), ...)
  }
}



