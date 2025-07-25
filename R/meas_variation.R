#' Measures of Quantitative Variation
#' 
#' @description
#' Probably the most famous measure of dispersion is the standard deviation, but there are more. This function provides a variety of measures and allows the creation of your own version.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/hhqMByH1vIo) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Measures/QuantitativeVariation.html)
#' 
#' 
#' @param data list or dataframe
#' @param levels dictionary, optional coding to use
#' @param measure c("std", "var", "mad", "madmed", "medad", "stddm", "cv", "cd", "own"), optional the measure to determine. Default is "std"
#' @param ddof float, optional option to adjust the division in standard deviation or variance with. Default is 1.
#' @param center c("mean", "median", "mode") or float, optional if measure is "own" the value to use as center. Default is "mean"
#' @param azs c("square", "abs"), optional if measure is "own" the way to avoid a zero sum. Either by squaring or absolute value
#' 
#' 
#' @returns
#' A dataframe with:
#' * *value*, the value of the measure
#' * *measure*, description of the measure
#' 
#' 
#' @details
#' 
#' **Standard Deviation** (std)
#' 
#' The formula used is:
#' \deqn{s = \sqrt{\frac{\sum_{i=1}^n \left(x_i - \bar{x}\right)^2}{n - d}}}
#' 
#' Where \eqn{d} is the offset specified at *ddof*. By default this is 1, giving the sample standard deviation.
#' 
#' **Variance** (var)
#' 
#' The formula used is:
#' 
#' \deqn{s^2 = \frac{\sum_{i=1}^n \left(x_i - \bar{x}\right)^2}{n - d}}
#' 
#' Where \eqn{d} is the offset specified at *ddof*. By default this is 1, giving the sample standard deviation.
#' 
#' **Mean Absolute Deviation** (mad)
#' 
#' The formula used is:
#' \deqn{MAD = \frac{\sum_{i=1}^n \left| x_i - \bar{x}\right|}{n}}
#' 
#' **Mean Absolute Deviation from the Median** (madmed)
#' 
#' The formula used is:
#' 
#' \deqn{MAD = \frac{\sum_{i=1}^n \left| x_i - \tilde{x}\right|}{n}}
#' 
#' Where \eqn{\tilde{x}} is the median
#' 
#' **Median Absolute Deviation** (medad)
#' 
#' The formula used is:
#' \deqn{MAD = MED\left(\left| x_i - \tilde{x}\right|\right)}
#' 
#' **Decile Standard Deviation**
#' 
#' The formula used is (Siraj-Ud-Doulah, 2018, p. 310):
#' \deqn{s_{dm} = \sqrt{\frac{\sum_{i=1}^n \left(x_i - DM\right)^2}{n - d}}}
#' 
#' Where DM is the decile mean.
#' 
#' **Coefficient of Variation** (cv)
#' 
#' The formula used is (Pearson, 1896, p. 277):
#' \deqn{CV = \frac{s}{\bar{x}}}
#' 
#' **Coefficient of Diversity** (cd)
#' 
#' The formula used is (Siraj-Ud-Doulah, 2018, p. 310):
#' \deqn{CD = \frac{s_{dm}}{DM}}
#' 
#' **Own**
#' it's possible to create one's own method. Decide on a specific center. Default options are the mean, median and mode. Then on either to sum the squared deviations or the absolute differences.
#' 
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to create a binned frequency table or a visualisation:
#' \code{\link{tab_frequency_bins}}, to create a binned frequency table.
#' \code{\link{vi_boxplot_single}}, for a Box (and Whisker) Plot.
#' \code{\link{vi_histogram}}, for a Histogram.
#' \code{\link{vi_stem_and_leaf}}, for a Stem-and-Leaf Display.
#' 
#' After this you might want some other descriptive measures:
#' \code{\link{me_mode_bin}}, for Mode for Binned Data.
#' \code{\link{me_mean}}, for different types of mean.
#' 
#' Or a perform a test:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' \code{\link{ts_z_os}}, for One-Sample Z Test.
#' 
#' 
#' @references
#' Pearson, K. (1896). Contributions to the mathematical theory of evolution. III. Regression, Heredity, and Panmixia. *Philosophical Transactions of the Royal Society of London*. (A.), 1896, 253-318.
#' 
#' Siraj-Ud-Doulah, M. (2018). Alternative measures of standard deviation coefficient of variation and standard error. *International Journal of Statistics and Applications, 8*(6), 309-315. https://doi.org/10.5923/j.statistics.20180806.04
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' studentDf = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' # Example 1: Numeric dataframe
#' ex1 = studentDf[['Gen_Age']]
#' me_variation(ex1)
#' 
#' # Example 2: Mean Absolute Deviation of a Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5)
#' me_variation(ex2, measure='mad')
#' 
#' @export
me_variation <- function(data, 
                         levels=NULL, 
                         measure="std", 
                         ddof=1, 
                         center="mean", 
                         azs="square"){
  
  if (is.null(levels)){dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  dataN = na.omit(dataN)
  n = length(dataN)
  if (measure=="std"){
    if (ddof==1){lbl = "standard deviation (sample)"}
    else if (ddof==0){lbl = "standard deviation (population)"}
    else {lbl=paste0("standard deviation corrected with ", ddof)}
    m = mean(dataN)
    res = (sum((dataN - m)**2)/(n-ddof))**0.5}
  else if (measure=="var"){
    if (ddof==1){lbl = "variance (sample)"}
    else if (ddof==0){lbl = "variance (population)"}
    else {lbl=paste0("variance corrected with ", ddof)}
    m = mean(dataN)
    res = sum((dataN - m)**2)/(n-ddof)}
  else if (measure=="mad"){
    lbl = "mean absolute deviation"
    m = mean(dataN)
    res = sum(abs(dataN - m))/n}
  else if (measure=="madmed"){
    lbl = "mean absolute deviation around median"
    m = median(dataN)
    res = sum(abs(dataN - m))/n}
  else if (measure=="medad"){
    lbl = "median absolute deviation"
    m = median(dataN)
    res = median(abs(dataN - m))}
  else if (measure=="cv"){
    lbl = "coefficient of variation"
    mu = mean(dataN)
    s = me_variation(dataN, measure="std", ddof=ddof)[1,1]
    res = s/mu}
  else if (measure=="stddm"){
    lbl = "standard deviation with decile mean"
    mu = me_mean(dataN, version="decile")
    res = (sum((dataN - mu)**2)/(n-ddof))**0.5}
  else if (measure=="cd"){
    lbl = "coefficient of deviation"
    mu = me_mean(dataN, version="decile")
    s = (sum((dataN - mu)**2)/(n-ddof))**0.5
    res = s/mu}
  else{
    if (center=="mean"){
      lbl = "mean"
      mu = mean(dataN)}
    else if (center=="median"){
      lbl = "median"
      mu = median(dataN)}
    else if (center=="mode"){
      lbl = "mode"
      mu = me_mode(dataN)[1,1]}
    else {
      lbl = center
      mu = center}
    
    if (azs=="square"){
      lbl = paste("sum squared deviation around ", lbl)
      res = sum((dataN - mu)**2)}
    else if (azs=="abs"){
      lbl = paste("sum absolute deviation around ", lbl)
      res = sum(abs(dataN - mu))
    }
  }
  results = data.frame(res, lbl)
  colnames(results)<-c("value", "measure")    
  return (results)
}



