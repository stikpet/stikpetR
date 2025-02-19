#' Simple Bar-Chart
#' 
#' @description 
#' A bar-chart is defined as “a graph in which bars of varying height with spaces between them are used to display data for variables defined by qualities or categories” (Zedeck, 2014, p. 20). 
#' 
#' A [YouTube](https://youtu.be/zT52FTyC6P8) video on pie charts.
#' 
#' @param data A vector or dataframe
#' @param varname Optional name for the variable
#' @param height Optional to indicate what the height should represent
#' 
#' @details 
#' 
#' The function uses the basic R's graphics library *barplot* function.
#' 
#' As a guideline for the size of the bar there is a rule of thumb known as the 'three quarter high rule' 
#' (Pitts, 1971). It means that the height of the vertical axis should be 3/4 of the length of the 
#' horizontal axis. So if the horizontal axis is 20 cm long, the vertical axis should be 
#' 3/4 * 20 = 15 cm high.
#' 
#' According to Singh (2009) vertical bars (instead of horizontal bars) are preferred since they are 
#' easier on the eye. However if you have long category names some names might become unreadable. 
#' A bar chart with the bars placed horizontally might then be preferred. 
#' 
#' One of the earliest found bar-charts from William Playfair (1786) has the bars placed horizontally. 
#' There is an earlier bar chart by Oresme (1486), but that is used more for a theoretical concept, 
#' than for descriptive statistics.
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
#' \code{\link{vi_cleveland_dot_plot}}, for Cleveland Dot Plot. 
#' \code{\link{vi_dot_plot}}, for Dot Plot.
#' \code{\link{vi_pareto_chart}}, for Pareto Chart.
#' \code{\link{vi_pie}}, for Pie Chart.
#' 
#' @references 
#' Oresme, N. (1486). *Tractatus de latitudinibus formarum*. (B. Pelacani da Parma, Ed.). Mathaeus Cerdonis.
#' 
#' Pitts, C. E. (1971). *Introduction to educational psychology: An operant conditioning approach*. Crowell.
#' 
#' Playfair, W. (1786). *The commercial and political atlas*. Debrett; Robinson; and Sewell.
#' 
#' Singh, G. (2009). *Map work and practical geography* (4th ed). Vikas Publishing House Pvt Ltd.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_bar_simple(ex1);
#' vi_bar_simple(ex1, varname="marital status", height="percent");
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_bar_simple(ex2);
#' 
#' @export
vi_bar_simple <- function(data, varname=NULL, height="count"){
  
  fr = table(data)
  
  if (height=="count") {
    barplot(fr, xlab=varname, ylab="Frequency");
  }
  else if (height=="percent"){
    n = sum(fr)
    barplot(fr/n*100, xlab=varname, ylab="Percent");
  }

}