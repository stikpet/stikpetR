#' Binned Frequency Table 
#' 
#' @description
#' Bins data and creates a frequency table with frequency density.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/-1kJrdDkImI) and frequency tables are also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tables/FrequencyTable.html)
#' 
#' @param data list or dataframe
#' @param nbins optional, either the number of bins to create, or a specific method from the *tab_nbins()* function. Default is "sturges"
#' @param bins optional dataframe with lower and upper bounds
#' @param incl_lower optional boolean, to include the lower bound, otherwise the upper bound is included. Default is True
#' @param adjust optional value to add  or subtract to guarantee all scores will fit in a bin
#' 
#' @returns
#' dataframe with:
#' 
#' \item{lower bound}{lower bound of class}
#' \item{upper bound}{upper bound of class}
#' \item{frequency}{count of scores in bin}
#' \item{frequency density}{count divided by bin range}
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to determine the number of bins you use:
#' \code{\link{tab_nbins}}, to determine the number of bins.
#' 
#' After this you might want to visualise the result:
#' \code{\link{vi_boxplot_single}}, for a Box (and Whisker) Plot.
#' \code{\link{vi_histogram}}, for a Histogram.
#' \code{\link{vi_stem_and_leaf}}, for a Stem-and-Leaf Display.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Numeric Dataframe
#' ex1a = df2['Gen_Age']
#' tab_frequency_bins(ex1a)
#' 
#' ex1b = df2['Gen_Age']
#' myBins = data.frame(c(0, 20, 25, 30), c(20, 25, 30, 120))
#' tab_frequency_bins(ex1b, bins=myBins)
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' tab_frequency_bins(ex2, adjust=0.1)
#' 
#' @export
tab_frequency_bins <- function(data, nbins="sturges", bins=NULL, incl_lower=TRUE, adjust=1){
  
  data = data.frame(data)
  
  #remove missing values
  data = na.omit(data)
  
  if (is.null(bins)){
    if (is.integer(nbins)){
      k = nbins}
    else{
      k = tab_nbins(data, method=nbins)}
    
    #determine minimum and maximum
    mx = max(data)
    mn = min(data)
    
    #increase maximimum if to include the lower bound
    if (incl_lower){mx = mx + adjust}
    #decrease minimum if to include the upper bound
    else{mn = mn - adjust}
    
    #determine range and width
    r = mx - mn
    h = r/k
    
    #create the bins
    freq = data.frame(matrix(nrow=0, ncol=4))
    colnames(freq)<-c("lower bound", "upper bound", "frequency", "frequency density")
    for (i in 0:(k - 1)){
      lb = mn + i*h
      ub = lb + h
      if (incl_lower){
        f = sum(data<ub) - sum(data<lb)}
      else{
        f = sum(data<=ub) - sum(data<=lb)}
      fd = f / (ub - lb)
      freq[nrow(freq) + 1,] = c(lb, ub, f, fd)
    }
  }
  else{
    k = nrow(bins)
    freq = data.frame(matrix(nrow=0, ncol=4))
    colnames(freq)<-c("lower bound", "upper bound", "frequency", "frequency density")
    for (i in 0:(k - 1)){
      lb = bins[(i+1),1]
      ub = bins[(i+1),2]
      if (incl_lower){
        f = sum(data<ub) - sum(data<lb)}
      else{
        f = sum(data<=ub) - sum(data<=lb)}
      fd = f / (ub - lb)
      freq[nrow(freq) + 1,] = c(lb, ub, f, fd)
    }
  }
  
  
  return (freq)
}



