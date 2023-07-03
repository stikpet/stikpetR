#' Box (and Whisker) Plot
#' 
#' @description
#' A box plot is a little more complex visualisation than a histogram. It shows the five quartiles (e.g. minimum, 1st quartile, median, 3rd quartile, and maximum). It can also be adjusted to show so-called outliers.
#' 
#' @param data list or dataframe
#' @param varname optional name to display on vertical axis
#' 
#' @returns
#' boxplot
#' 
#' @details
#' This was actually a 'range chart' (Spear, 1952, p. 166) but somehow it is these days referred to as a box-and-whisker plot as named by Tukey (1977, p. 39)
#' 
#' The function uses the **boxplot()** function from the *graphics* library. If you want to modify more things you might want to use that function.
#' 
#' @references
#' 
#' Spear, M. E. (1952). *Charting statistics*. McGraw-Hill.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
vi_boxplot_single <- function(data, varname=NULL){
  if (is.null(varname)){varname=deparse(substitute(data))}
  
  data = data.frame(data)
  
  boxplot(data, ylab=varname, horizontal=TRUE)
  
}