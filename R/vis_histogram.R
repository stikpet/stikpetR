#' Histogram
#' 
#' @description
#' A histogram is a bit like a bar chart for a scale variable. You would create some bins, and then plot these as bars.
#' 
#' @param data list or dataframe
#' @param xlbl optional label for the horizontal axis
#' @param ylbl optional label for the vertical axis
#' @param ... other parameters for use in hist function
#' 
#' @returns
#' The histogram
#' 
#' @details
#' This function is just using some defaults for the **hist()** function from R's *graphics* library.
#' 
#' To set the bins, the *breaks* argument can be used. This could be a pre-set number based on a calculation, a specific rule (e.g. bins="sturges"), or a list with the cut-off points.
#' 
#' If your bins are of equal width, a true histogram than actually should show frequency densities (Pearson, 1895, p. 399). These are the frequencies divided by the bin-width. This can be done using *freq=FALSE* parameter.
#' 
#' @references
#' Pearson, K. (1895). Contributions to the mathematical theory of evolution. II. Skew variation in homogeneous material. *Philosophical Transactions of the Royal Society of London. (A.)*, 186, 343–414. doi:10.1098/rsta.1895.0010
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: dataframe
#' ex1 = df2['Gen_Age']
#' vi_histogram(ex1);
#' vi_histogram(ex1, freq=FALSE);
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' vi_histogram(ex2);
#' 
#' @export
vi_histogram <- function(data, xlbl=NULL, ylbl=NULL, ...){
  
  data = data.frame(data)
  data = as.numeric(data[,1])
  hist(data, xlab = xlbl, ylab=ylbl, ...)
  
}