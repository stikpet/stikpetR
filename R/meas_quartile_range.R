#' Interquartile Range, Semi-Interquartile Range and Mid-Quartile Range
#' 
#' @param data dataframe with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param range the specific range to determine
#' @param method the method to use to determine the quartiles
#' @returns
#' A dataframe with:
#' \item{Q1}{the first (lower) quartile}
#' \item{Q3}{the third (upper/higher) quartile}
#' \item{range}{the range determined}
#' 
#' @description 
#' Three different ranges that can be used with quartiles. The Interquartile Range (IQR) is 
#' simply the third minus the first quartile.
#' 
#' The Semi-Interquartile range (a.k.a. Quartile Deviation) divides the IQR by 2.
#' 
#' The Mid-Quartile Range adds the two quartiles and then divides by 2.
#' 
#' A special case is the H-spread, which is if Hinges are used the IQR. This can be obtained by 
#' setting *method="tukey"*.
#' 
#' The function uses the *me_quartiles* function and any of the methods from that function can be used.
#' 
#' @details 
#' The formula used for the Interquartile Range is:
#' \deqn{IQR = Q_3 - Q_1}
#' This can be obtained by setting *range="iqr"*.
#' 
#' The IQR is mentioned in Galton (1881, p. 245) and the H-spread in Tukey (1977, p. 44).
#' 
#' The H-spread can be obtained by setting *range="iqr"* and *method="tukey"*.
#' 
#' The formula used for the Semi-Interquartile Range (Quartile Deviation) is (Yule, 1911, p. 147):
#' \deqn{SIQR = \frac{Q_3 - Q_1}{2}}
#' This can be obtained by setting *range="siqr"* or *range="qd"*.
#' 
#' The formula for the mid-quartile range used is:
#' \deqn{MQR = \frac{Q_3 + Q_1}{2}}
#' This can be obtained by setting *range="mqr"*.
#' 
#' This formula can be found in Parzen (1980, p. 19), but there are probably older references.
#' 
#' @examples 
#' data = runif(n=10, min=1, max=50)
#' me_quartile_range(data)
#' me_quartile_range(data, method="tukey")
#' me_quartile_range(data, range="siqr")
#' me_quartile_range(data, range="mqr")
#' 
#' @seealso 
#' For more details on the quartiles calculation see \code{\link{me_quartiles}}.
#' 
#' @references 
#' Galton, F. (1881). Report of the anthropometric committee. Report of the British Association for the Advancement of Science, 51, 225–272.
#' 
#' Parzen, E. (1980). *Data modeling using quantile and density-quantile functions*. Institute of Statistics, Texas A&M University.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' Yule, G. U. (1911). *An introduction to the theory of statistics*. Charles Griffin.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
me_quartile_range <- function(data, levels=NULL,
                              range=c("iqr", "siqr", "qd", "mqr"), 
                              method="cdf"){
  
  if (length(range)>1){range="iqr"}
  
  if (is.null(levels)){
    dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  Qs = me_quartiles(dataN, method=method)
  q1 = Qs$q1
  q3 = Qs$q3
  
  if (range=="iqr") {
    r = q3 - q1
    if (method=="tukey" || method=="inclusive" || method=="tukey" || method=="vining" || method=="hinges") {rName = "Hspread"}
    else{rName = "IQR"}
  }
  else if (range=="siqr" || range=="qd"){
    r = (q3 - q1)/2
    rName = "SIQR"
  }
  else if (range=="mqr"){
    r = (q3 + q1)/2
    rName = "MQR"
  }
  
  results = data.frame(q1, q3, r)
  colnames(results) = c("Q1", "Q3", rName)
  
  
  return (results)
  
}