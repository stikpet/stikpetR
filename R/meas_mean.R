#' Mean
#' 
#' @description 
#' Different types of means can be determined using this function. 
#' 
#' The mean is a measure of central tendency, to indicate the center.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/ch2zRuLpw_A) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Measures/mean.html)
#' 
#' @param data, vector or dataframe with scores as numbers
#' @param levels : list, optional coding to use
#' @param version, optional mean to calculate. Either `"arithmetic"` (default), `"winsorized"`, `"trimmed"`, `"windsor"`, `"truncated"`, `"olympic"`, `"geometric"`, `"harmonic"`, `"midrange"`, or `"decile"`
#' @param trimProp, optional to indicate the total proportion to trim. Default at 0.1 i.e. 0.05 from each side.
#' @param trimFrac, optional parameter to indicate what to do if trimmed amount is non-integer. Either `"down"` (default), `"prop"`, `"linear"`
#' 
#' @returns 
#' res, the value of the mean
#' 
#' @details 
#' **Arithmetic Mean**
#' 
#' One of the three Pythagorean means, and the mean most people would assume if you ask them to calculate the mean.
#' 
#' It is the fulcrum of the distribution (Weinberg & Schumaker, 1962, p.19). One reference can for example be found in Aristotle (384-322 BC) (1850, p. 43). 
#' 
#' The formula:
#' \deqn{\bar{x} = \frac{\sum_{i=1}^n x_i}{n}}
#' 
#' **Harmonic Mean**
#' 
#' The second of the three Pythagorean means:
#' \deqn{H = \frac{n}{\sum_{i=1}^n \frac{1}{x_i}}}
#' 
#' **Geometric Mean**
#' 
#' The third of the three Pythagorean means:
#' \deqn{G = e^{\frac{1}{n}\times\sum_{i=1}^n \ln\left(x_i\right)}}
#' 
#' **Olympic Mean**
#' 
#' Simply ignore the maximum and minimum (only once) (Louis et al., 2023, p. 117):
#' \deqn{OM = \frac{\sum_{i=2}^{n-1} x_i}{n - 2}}
#' 
#' **Mid Range**
#' 
#' The average of the maximum and minimum (Lovitt & Holtzclaw, 1931, p. 91):
#' \deqn{MR = \frac{\min x + \max x}{2}}
#' 
#' **Trimmed**
#' 
#' With a trimmed (Windsor/Truncated) mean we trim a fixed amount of scores from each side (Tukey, 1962, p. 17). Let \eqn{p_t} be the proportion to trim, we then need to trim \eqn{n_t = \frac{p_t\times n}{2}} from each side.
#' 
#' If this \eqn{n_t} is an integer there isn't a problem, but if it isn't we have options. The  first option is to simply round down, i.e. \eqn{n_l = \lfloor n_t\rfloor}. The trimmed mean is then:
#' \deqn{\bar{x}_t = \frac{\sum_{i=n_t+1}^{n - n_l + 1} x_i}{n - 2\times n_l}}
#' 
#' This is used if *trimFrac = "down"* is set.
#' 
#' We could also use linear interpolation based on the number of scores to trim. We missed out on: \eqn{f = n_t - n_l} on each side. So the first and last value we do include should only count for \eqn{1 - f} each. The trimmed mean will then be:
#' \deqn{\bar{x}_t = \frac{\left(x_{n_t + 1} + x_{n - n_l + 1}\right)\times\left(1 - f\right) + \sum_{i=n_l+2}^{n - n_l} x_i}{n - 2\times n_t}}
#' 
#' This is used if *trimFrac = "prop"* is set.
#' 
#' Alternative, we could take the proportion itself and use linear interpolation on that. The found \eqn{n_l} will be \eqn{p_1 = \frac{n_l \times 2}{n}} of the the total sample size. While if we had rounded up, we had used \eqn{p_2 = \frac{\left(n_l + 1\right)\times 2}{n}} of the the total sample size. Using linear interpolation we then get:
#' \deqn{\bar{x}_t = \frac{p_t - p_1}{p_2 - p_1}\times\left(\bar{x}_{th}-\bar{x}_{tl}\right) + \bar{x}_{tl}}
#' Where \eqn{\bar{x}_{tl}} is the trimmed mean if \eqn{p_1} would be used as a trim proportion, and \eqn{\bar{x}_{th}} is the trimmed mean if \eqn{p_2} would be used.
#' 
#' This is used if *trimFrac = "linear"* is set.
#' 
#' **Winsorized Mean**
#' 
#' Similar as with a trimmed mean, but now the data is not removed, but replaced by the value equal to the nearest value that is still included (Winsor as cited in Dixon, 1960, p. 385).
#' \deqn{W = \frac{n_l \times \left(x_{n_l + 1} + x_{n - n_l}\right) + \sum_{n_l + 1}^{n - n_l} x_i}{n}}
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
#' \code{\link{me_variation}}, for different Measures of Quantitative Variation.
#' 
#' Or a perform a test:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' \code{\link{ts_z_os}}, for One-Sample Z Test.
#' 
#' 
#' @references 
#' Aristotle. (1850). *The nicomachean ethics of Aristotle* (R. W. Browne, Trans.). Henry G. Bohn.
#' 
#' Dixon, W. J. (1960). Simplified estimation from censored normal samples. *The Annals of Mathematical Statistics, 31*(2), 385–391. https://doi.org/10.1214/aoms/1177705900
#' 
#' Louis, P., Núñez, M., & Xefteris, D. (2023). Trimming extreme reports in preference aggregation. *Games and Economic Behavior, 137*, 116–151. https://doi.org/10.1016/j.geb.2022.11.003
#' 
#' Lovitt, W. V., & Holtzclaw, H. F. (1931). *Statistics*. Prentice Hall.
#' 
#' Tukey, J. W. (1962). The future of data analysis. *The Annals of Mathematical Statistics, 33*(1), 1–67. https://doi.org/10.1214/aoms/1177704711
#' 
#' Weinberg, G. H., & Schumaker, J. A. (1962). *Statistics An intuitive approach*. Wadsworth Publishing.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Numeric dataframe
#' ex1 = df2['Gen_Age']
#' me_mean(ex1)
#' 
#' #Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' me_mean(ex2)
#' 
#' #Example 3: Ordinal Pandas Series
#' ex3 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' me_mean(ex3, levels=order)
#' 
#' @export
me_mean <- function(data, 
                    levels=NULL, 
                    version="arithmetic", 
                    trimProp=0.1, 
                    trimFrac="down"){
  
  if (is.null(levels)){dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  data = unlist(na.omit(dataN))
  
  if (version=="arithmetic"){
    res = mean(data)}
  else if (version %in% c("winsorized", "trimmed", "windsor", "truncated")){
    data = sort(data)
    n = length(data)
    nt1 = n*trimProp/2
    nl = floor(nt1)
    
    if (version=="winsorized"){
      data[1:nl] = data[nl+1]
      data[(n-nl):n] = data[n-nl]
      res = mean(data)}
    else{
      if (trimFrac=="down"){
        res = mean(data[(nl+1):(n-nl)])}
      else if (trimFrac=="prop"){
        fr = nt1 - nl
        res = (data[nl+1]*(1 - fr) + data[n-nl]*(1 - fr) + sum((data[(nl+2):(n-nl-1)])))/(n - nt1*2)}
      else if (trimFrac=="linear"){
        p1 = nl*2/n
        p2 = (nl + 1)*2/n
        m1 = mean(data[(nl+1):(n-nl)])
        m2 = mean(data[(nl+2):(n-nl-1)])
        res = (trimProp - p1)/(p2 - p1)*(m2 - m1)+m1}
    }
  }
  else if (version=="olympic"){
    res = (sum(data) - max(data) - min(data))/(length(data)-2)}
  else if (version=="geometric"){
    n = length(data)
    res = exp(sum(log(data))/n)}
  else if (version=="harmonic"){
    n = length(data)
    res = n/sum(1/data)}
  else if (version=="midrange"){
    res = (max(data) + min(data))/2}
  else if (version=="decile"){
    qs = me_quantiles(data, k=10)
    dm = 0
    for (i in 2:10){
      dm = dm + qs[i]
    }
    res = dm/9
    
  }
  
  return (res)
}