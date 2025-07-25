#' Box (and Whisker) Plot
#' 
#' @description
#' A box plot is a little more complex visualisation than a histogram. It shows the five quartiles (e.g. minimum, 1st quartile, median, 3rd quartile, and maximum). It can also be adjusted to show so-called outliers.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/LAlnmEVuZXw) and the visualisation is described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/boxPlot.html)
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
#' \code{\link{vi_histogram}}, for a Histogram.
#' \code{\link{vi_stem_and_leaf}}, for a Stem-and-Leaf Display.
#' 
#' @references
#' Spear, M. E. (1952). *Charting statistics*. McGraw-Hill.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: dataframe
#' ex1 = df2['Gen_Age']
#' vi_boxplot_single(ex1);
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' vi_boxplot_single(ex2);
#' 
#' @export
vi_boxplot_single <- function(data, varname=NULL){
  if (is.null(varname)){varname=deparse(substitute(data))}
  
  data = data.frame(data)
  
  boxplot(data, ylab=varname, horizontal=TRUE)
  
}



