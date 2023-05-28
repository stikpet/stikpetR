#' Quantiles
#' 
#' @param data the scores from which to determine the quantiles
#' @param k the k-quantiles to determine (default is 4)
#' @param method specific method to use to determine quantiles (default is 'cdf')
#' @return a named vector with the different quantiles
#' 
#' @description 
#' Quantiles split the data into k sections, each containing n/k scores. They can be seen as a 
#' generalisation of various 'tiles'. For example 4-quantiles is the same as the quartiles, 
#' 5-quantiles the same as quintiles, 100-quantiles the same as percentiles, etc.
#' 
#' Quite a few different methods exist to determine these. See the details for more information.
#' 
#' @details 
#' To determine the quantile, first determine the index of the quantile (\eqn{iQ_x}). 
#' Unfortunately there are quite some different methods for this. Below a few:
#' 
#' |method|\eqn{iQ_x}|
#' |----------------|----------------|
#' |SAS-1 (HF-4)|\eqn{n\times p_i}|
#' |SAS-2 (HF-3)|\eqn{\lfloor n\times p_i \rceil}|
#' |SAS-3 (HF-1)|\eqn{\lceil n\times p_i \rceil}|
#' |Minitab (HF-6 / Snedecor)|\eqn{\left(n + 1\right)\times p_i}|
#' |CDF (SAS-5 / HF-2)|\eqn{\begin{cases} n\times p_i+\frac{1}{2} & \text{ if } n\times p_i= \left [ n\times p_i \right ]\\ \lceil n\times p_i\rceil & \text{ if } n\times p_i \neq \left [ n\times p_i \right ] \end{cases}} |
#' |Hogg-Ledolter-1|\eqn{\begin{cases} n\times p_i+\frac{1}{2} & \text{ if } n\times p_i +\frac{1}{2}= \left [ n\times p_i +\frac{1}{2} \right ]\\ \lfloor n\times p_i+\frac{1}{2}\rfloor +\frac{1}{2}& \text{ if } n\times p_i +\frac{1}{2}\neq \left [ n\times p_i +\frac{1}{2}\right ]  \end{cases}} |
#' |Hogg-Ledolter-2 (HF-5)|\eqn{n\times p_i + \frac{1}{2}} |
#' |Excel (HF-7)|\eqn{\left(n - 1\right)\times p_i + 1} |
#' |HF-8| \eqn{\left(n + \frac{1}{3}\right)\times p_i + \frac{1}{3}} |
#' |HF-9| \eqn{\left(n + \frac{1}{4}\right)\times p_i + \frac{3}{8}} |
#' 
#' Then note that the first and last should be the minimum and maximum of the data, i.e.
#' \eqn{iQ_0 = 1} and \eqn{iQ_k = n}
#' 
#' If \eqn{iQ_x} is an integer, then:
#' \deqn{Q_x = x_{iQ_x}}
#' 
#' Otherwise, use linear interpolation:
#' \deqn{Q_x = \left(iQ_x - iLQ_x\right)\times\left(x_{iHQ_x} - x_{iLQ_x}\right) + x_{iLQ_x}}
#' With:
#' \deqn{iHQ_x = \lceil iQ_x \rceil = iLQ_x + 1}
#' \deqn{iLQ_x = \lfloor iQ_x \rfloor}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_i} the i-the score of x after the scores have been sorted
#' \item \eqn{n} the sample size (i.e. the number of scores)
#' \item \eqn{\left[\dots\right]} round to the nearest integer, if 0.5 round up
#' \item \eqn{\lceil\dots\rceil} round up to the nearest integer
#' \item \eqn{\lfloor\dots\rfloor} round down to the nearest integer
#' \item \eqn{\lfloor\dots\rceil} round to the nearest even integer
#' }
#' 
#' For the quartiles, simply set 'k=4'. The p values will then be 0/4, 1/4, 2/4, 3/4, and 4/4.
#' 
#' **SAS method 1**
#' 
#' Simply uses \eqn{iQ_x = n\times p}.It is the same as Hyndman and Fan method 4.
#' 
#' **SAS method 3** 
#' 
#' simply rounds each of the results up to the nearest integer. It is the same as Hyndman and Fan method 1.
#' 
#' **SAS method 2** 
#' 
#' is very similar as method 3, except if \eqn{n\times p} ends with 0.5. It then suggests to round to 
#' the nearest even integer. So for example 3.5 gets rounded to 4, but 2.5 gets rounded to 2.
#' It is the same as Hyndman and Fan method 3.
#' 
#' **SAS method 4** 
#' 
#' first adds 1 to the sample size, so uses \eqn{\left(n + 1\right)\times p}. 
#' 
#' The SAS method 4 is also used by **Minitab** and can also be found in Snedecor (1940, p. 43).
#' It is the same as Hyndman and Fan method 6.
#' 
#' **SAS method 5** 
#' 
#' the same as SAS method 3, except \eqn{n\times p} is an integer. In that case it simply 
#' takes the average of the index and the next value i.e. \eqn{n\times p + 1}.
#' It is the same as Hyndman and Fan method 2, and Langford (2006) refers to this as the 'CDF' method.
#' 
#' The methods used by SAS are found in the procedures guide (SAS, 2006, p. 626).
#' 
#' **Hogg and Ledolter**
#' 
#' Hogg and Ledolter (1992, p. 21) use that for the median \eqn{\frac{n}{2} + \frac{1}{2}} is used.
#' In their first version they suggest to simply average the two values the index falls between 
#' if it is not an integer.
#' This means we can't use the linear interpolation directly. We (in case it is not an integer), first 
#' simply round down and then add 0.5.  
#' It could be that this method was also described in Hazen (1914), but I couldn't find the exact 
#' page number in that document.
#' 
#' The second version of Hogg and Ledolter is to simply use the linear interpolation. 
#' This is the same as Hyndman and Fan method 5.
#' 
#' **Excel**
#' 
#' It appears that MS Excel has a unique method. It can also be found in Freund and Perles (1987, p. 201).
#' This is the same as Hyndman and Fan method 7.
#' 
#' **Hyndman-Fan**
#' 
#' Hyndman and Fan (1996) discuss 9 different methods for quantiles, which can be used to determine quartiles.
#' The R stats library's quantile() function uses their numbering. 
#' HF-1 = SAS-3, HF-2 = CDF / SAS-5, HF-3 = SAS-2, HF-4 = SAS-1, HF-5 = Hogg-Ledolter v2, 
#' HF-6 = Minitab / SAS-4, and HF-7 = Excel.
#' 
#' This leaves their 8th and 9th definitions which are based on approximating the incomplete 
#' beta function. 
#' 
#' @section Alternatives:
#' 
#' The R *stats* library has a very similar function *quantile()*. They use the numbering from 
#' Hyndman and Fan.
#' 
#' @examples 
#' 
#' data = seq(10)
#' me_quantiles(data, 4, method="sas1")
#' me_quantiles(data, 4, method="sas2")
#' me_quantiles(data, 4, method="sas3")
#' me_quantiles(data, 4, method="minitab")
#' me_quantiles(data, 4, method="cdf")
#' me_quantiles(data, 4, method="hl1")
#' me_quantiles(data, 4, method="hl2")
#' me_quantiles(data, 4, method="excel")
#' me_quantiles(data, 4, method="hf8")
#' me_quantiles(data, 4, method="hf9")
#' 
#' @seealso 
#' A special case of quantiles are the quartiles, see \code{\link{me_quartiles}} for more info on those.
#' 
#' @references 
#' Freund, J. E., & Perles, B. M. (1987). A new look at quartiles of ungrouped data. *The American Statistician, 41*(3), 200–203. https://doi.org/10.1080/00031305.1987.10475479
#' 
#' Hazen, A. (1914). Storage to be provided in impounding municipal water supply. *Transactions of the American Society of Civil Engineers, 77*(1), 1539–1640. https://doi.org/10.1061/taceat.0002563
#' 
#' Hogg, R. V., & Ledolter, J. (1992). *Applied statistics for engineers and physical scientists* (2nd int.). Macmillan.
#' 
#' Hyndman, R. J., & Fan, Y. (1996). *Sample quantiles in statistical packages. The American Statistician, 50*(4), 361–365. https://doi.org/10.2307/2684934
#' 
#' Langford, E. (2006). Quartiles in elementary statistics. *Journal of Statistics Education, 14*(3), 1–17. https://doi.org/10.1080/10691898.2006.11910589
#' 
#' SAS. (1990). *SAS procedures guide: Version 6* (3rd ed.). SAS Institute.
#' 
#' Snedecor, G. W. (1940). *Statistical methods applied to experiments in agriculture and biology* (3rd ed.). The Iowa State College Press.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
me_quantiles <- function(data, k=4, 
                         method=c("cdf", "hl1", "hl2", "minitab", "excel", "sas1", 
                                  "sas2", "sas3", "hf8", "hf9")){
  
  sData = sort(na.omit(data))
  
  p = c(0, seq(k)/k)
  n = length(sData)
  
  if (method=="sas1") {
    h = n*p
  }
  else if (method=="sas2"){
    h = rep(0, k+1)
    for (i in 1:k+1) {
      if (n*p[i]*2 %% 1 ==0) {
        if (floor(n*p[i]) %% 2 == 0) {
          h[i] = floor(n*p[i])
        }
        else{
          h[i] = ceiling(n*p[i])
        }
      }
      else{
        h[i] = round(n*p[i])
      }
    }
  }
  else if (method=="sas3"){
    h = ceiling(n*p)
  }
  else if (method=="minitab"){
    h = (n + 1)*p
  }
  else if (method=="cdf"){
    h = rep(0, k+1)
    for (i in 1:k+1) {
      if (n*p[i]==round(n*p[i])) {
        h[i] = n*p[i] + 1/2
      }
      else{
        h[i] = ceiling(n*p[i])
      }
    }
  }
  else if (method=="hl1"){
    h = rep(0, k+1)
    for (i in 1:k+1) {
      if (n*p[i] + 1/2 == round(n*p[i] + 1/2)) {
        h[i] = n*p[i] + 1/2
      }
      else{
        h[i] = floor(n*p[i] + 1/2) + 1/2
      }
    }
  }
  else if (method=="hl2"){
    h = n*p + 1/2
  }
  else if (method=="excel"){
    h = (n - 1)*p + 1
  }
  else if (method=="hf8"){
    h = (n + 1/3)*p + 1/3
  }
  else if (method=="hf9"){
    h = (n + 1/4)*p + 3/8
  }
  
  #Set 0% and 100% to minimum and maximum value
  h[1] = 1
  h[k + 1] = n
  
  q = rep(0, k+1)
  #use linear interpolation if needed.
  for (i in 1:(k+1)) {
    if (h[i] %% 1 != 0) {
      hh = ceiling(h[i])
      lh = floor(h[i])
      q[i] = (h[i] - lh)*(sData[hh] - sData[lh]) + sData[lh]
    }
    else{
      q[i] = sData[h[i]]
    }
  }
  
  names(q) = paste0(p*100, "%")
  return (q)
}