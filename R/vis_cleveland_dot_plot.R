#' Cleveland Dot Plot
#' 
#' @description 
#' A Cleveland dot plot (Cleveland & McGill, 1987) is a bar chart where instead of bars 
#' a dot is placed at the center of the top of the bar (and then the bars removed). 
#' It is a dot plot only showing the top dot.This requires less ink. 
#' 
#' The function simply uses the dotplot() function from the lattice library.
#' 
#' A video on (Cleveland) dot plots is available [here](https://youtu.be/qs1nh0CMiIY).
#' 
#' This function is shown in this [YouTube video](https://youtu.be/K1Jb7XCDcDg) and the visualisation is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/DotPlot.html)
#' 
#' @param data the data from which to create the plot
#' @param size the size of the dots (default is 2)
#' @return chart the Cleveland dot plot
#' 
#' @section Before, After and Alternatives:
#' Before the visualisation you might first want to get an impression using a frequency table:
#' \code{\link{tab_frequency}}
#' 
#' After visualisation you might want some descriptive measures:
#' \code{\link{me_mode}}, for the mode. 
#' \code{\link{me_qv}}, for Measures of Qualitative Variation.
#' 
#' or perform a test:
#' \code{\link{ts_pearson_gof}}, for Pearson Chi-Square Goodness-of-Fit Test. 
#' \code{\link{ts_freeman_tukey_gof}}, for Freeman-Tukey Test of Goodness-of-Fit. 
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_mod_log_likelihood_gof}}, for Mod-Log Likelihood Test of Goodness-of-Fit. 
#' \code{\link{ts_multinomial_gof}}, for Multinomial Goodness-of-Fit Test. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test. 
#' 
#' Alternatives for this visualisation could be:
#' \code{\link{vi_bar_simple}}, for Simple Bar Chart. 
#' \code{\link{vi_dot_plot}}, for Dot Plot.
#' \code{\link{vi_pareto_chart}}, for Pareto Chart.
#' \code{\link{vi_pie}}, for Pie Chart.
#' 
#' @references 
#' Cleveland, W. S., & McGill, R. (1984). Graphical perception: Theory, experimentation, and application to the development of graphical methods. *Journal of the American Statistical Association, 79*(387), 531–554. https://doi.org/10.2307/2288400
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @importFrom lattice dotplot
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_cleveland_dot_plot(ex1);
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_cleveland_dot_plot(ex2);
#' 
#' @export
vi_cleveland_dot_plot <- function(data, size=2){
  
  chart = lattice::dotplot(data, horizontal = FALSE, cex=size)
  
  return(chart)
}