# Suppress R CMD check note for ggplot variable bindings
utils::globalVariables(c("Var1", "Var2", "Freq"))

#' Back-to-Back Histogram
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 after_stat
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' 
#' @description 
#' This function creates a simple back-to-back histogram. This is sometimes also referred to as a Pyramid chart or a Dual-Sided histogram (Jelen, 2005).
#' 
#' The back-to-back histogram together with back-to-back stem-and-leaf and split box-plots are described as effective ways to compare two distributions by Lane and Sándor (2009).
#' 
#' The function is shown in this [YouTube video](https://youtu.be/nJuel_XSCdo) and the visualisation is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/histogram.html)
#' 
#' See **vi_histogram()** for more details on histograms.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with the two categories to use from bin_field. If not set the first two found will be used
#' @param bins optional either a list with parameters to pass to tab_frequency_bins, a dataframe with upper and lower bounds, or simply NULL (default)
#' @param show c('count', 'relative'), show either counts or relative count on vertical scale
#' @param density c('auto', FALSE, TRUE), show (relative) frequency or use (relative) frequency density.
#' @param xlbl : string, optional label for the field axis, if not set the name of the field is used.
#' @param rotate : bool, optional rotate the bars so they appear horizontal. Default is True
#' @param colors vector, two colors, one for each category, default is c(rgb(0, 0, 1, 1/4), rgb(1, 0, 0, 1/4))
#' @param ... optional additional parameters to pass to geom_histogram
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
vi_histogram_b2b <- function(catField, scaleField, categories=NULL, bins=NULL, 
                             show='count', density='auto', xlbl=NULL, rotate=TRUE, 
                             colors=c('orange', 'blue'), ...){
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
  
  ##determine bins overall
  freq_table <- switch(
    class(bins)[1],   # use the first class if bins has multiple
    "NULL" = tab_frequency_bins(allScores),
    "list" = do.call(tab_frequency_bins, c(list(data = allScores), bins)),
    "data.frame" = tab_frequency_bins(allScores, bins = bins),
    stop("Invalid value for 'bins'")
  )
  
  bins = c(freq_table[, "lower bound"], max(freq_table["upper bound"]))
  
  if (is.null(xlbl)){
    xlbl = deparse(substitute(scaleField))
  }
  if (density=='auto'){
    if (all(diff(bins) == diff(bins)[1])){density=FALSE}
    else {density=TRUE}}
  
  if (show=="count" && density==FALSE){
    ylabel = 'frequency'
    hist1 = geom_histogram(data = df[df$category==cat1, ],
                           aes(x = score, y = after_stat(count), fill = category), 
                           breaks=bins, ...)
    hist2 = geom_histogram(data = df[df$category==cat2, ],
                           aes(x = score, y = -after_stat(count), fill = category), 
                           breaks=bins, ...)      
  }
  if (show=="count" && density==TRUE){
    ylabel = 'frequency density'
    hist1 = geom_histogram(data = df[df$category==cat1, ],
                           aes(x = score, y = after_stat(density*sum(count)), fill = category), 
                           breaks=bins, ...)
    hist2 = geom_histogram(data = df[df$category==cat2, ],
                           aes(x = score, y = -after_stat(density*sum(count)), fill = category), 
                           breaks=bins, ...)
  }
  if (show=="relative" && density==FALSE){
    ylabel = 'percent'
    hist1 = geom_histogram(data = df[df$category==cat1, ],
                           aes(x = score, y = after_stat(count/sum(count)*100), fill = category),
                           breaks=bins, ...)
    hist2 = geom_histogram(data = df[df$category==cat2, ],
                           aes(x = score, y = -after_stat(count/sum(count)*100), fill = category), 
                           breaks=bins, ...)
  }
  if (show=="relative" && density==TRUE){
    ylabel = 'probability density'
    hist1 = geom_histogram(data = df[df$category==cat1, ],
                           aes(x = score, y = after_stat(density), fill = category),
                           breaks=bins, ...)
    hist2 = geom_histogram(data = df[df$category==cat2, ],
                           aes(x = score, y = -after_stat(density), fill = category),
                           breaks=bins, ...)
  }
  plot = ggplot() + hist1 + hist2 + labs(x=xlbl, y=ylabel) + 
    scale_fill_manual(values = setNames(c(colors[1], colors[2]), c(cat1, cat2)))
  if (rotate){plot = plot + coord_flip()}
  
  plot
}