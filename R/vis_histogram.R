#' Histogram
#' 
#' @description
#' A histogram is a bit like a bar chart for a scale variable. You would create some bins, and then plot these as bars.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/v13UlQvbOvs) and the visualisation is described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/histogram.html)
#' 
#' @param data list or dataframe
#' @param show c('count', 'relative'), show either counts or relative count on vertical scale
#' @param density c('auto', FALSE, TRUE), show (relative) frequency or use (relative) frequency density.
#' @param xlbl string, optional label for horizontal axis (default is variable name of data)
#' @param ... other parameters for use in plot function
#' 
#' @returns
#' The histogram
#' 
#' @details
#' This function is just using some defaults for the **hist()** function from R's *graphics* library.
#' 
#' To set the bins, the *breaks* argument can be used. This could be a pre-set number based on a calculation, a specific rule (e.g. bins="sturges"), or a list with the cut-off points.
#' 
#' If your bins are of equal width, a true histogram than actually should show frequency densities (Pearson, 1895, p. 399). These are the frequencies divided by the bin-width. This can be done using *density=TRUE* parameter.
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to create a binned frequency table
#' \code{\link{tab_frequency_bins}}, to create a binned frequency table.
#' 
#' After this you might want some descriptive measures:
#' \code{\link{me_mode_bin}}, for Mode for Binned Data.
#' \code{\link{me_mean}}, for different types of mean.
#' \code{\link{me_variation}}, for different Measures of Quantitative Variation.
#' 
#' Or a perform a test:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' \code{\link{ts_z_os}}, for One-Sample Z Test.
#' 
#' Alternative Visualisations:
#' \code{\link{vi_boxplot_single}}, for a Box (and Whisker) Plot.
#' \code{\link{vi_stem_and_leaf}}, for a Stem-and-Leaf Display.
#' 
#' @references
#' Pearson, K. (1895). Contributions to the mathematical theory of evolution. II. Skew variation in homogeneous material. *Philosophical Transactions of the Royal Society of London. (A.)*, 186, 343-414. doi:10.1098/rsta.1895.0010
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_histogram <- function(data, show='count', density='auto', xlbl=NULL, ...){
  
  data = na.omit(data.frame(data))
  data = as.numeric(data[,1])
  
  if (is.null(xlbl)){
    xlbl = deparse(substitute(data))}
    
  h <- suppressWarnings(hist(data, plot = FALSE, ...))
  
  if (density=='auto'){
    if (h$equidist){density=FALSE}
    else {density=TRUE}}
  
  if (show=="count" && density==FALSE){
    yLabel = "frequency"
    histFreq = TRUE
    # set densities to counts, in case unequal bin-widths
    h$density = h$counts
  }
  
  else if (show=="relative" && density==FALSE){                
    h$counts <- h$counts / sum(h$counts) * 100
    # set densities to counts, in case unequal bin-widths
    h$density = h$counts
    yLabel = "percent"
    histFreq = TRUE        
  }
  
  else if (show=="count" && density==TRUE){
    h$density <- h$density * sum(h$counts)
    yLabel = "frequency density"
    histFreq = FALSE        
  } 
  
  else if (show=="relative" && density==TRUE){
    yLabel = "relative frequency density"
    histFreq = FALSE
  } 
  
  suppressWarnings(plot(h, freq=histFreq, xlab=xlbl, ylab=yLabel, ...))
  
}



