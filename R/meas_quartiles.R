#' Quartiles / Hinges
#' 
#' @param data ordinal input data
#' @param method optional which method to use to calculate quartiles. Default is "cdf", see details for more info.
#' @returns
#' A dataframe with:
#' \item{Q1}{the first (lower) quartile}
#' \item{Q3}{the third (upper/higher) quartile}
#' 
#' @description 
#' The quartiles are at quarters of the data. The median is at 50%, and the quartiles at 25% and 75%. Note that there are 
#' five quartiles, the minimum value is the 0-quartile, at 25% the first (or lower) quartile, at 50% the median a.k.a. the 
#' second quartile, at 75% the third (or upper) quartile, and the maximum as the fourth quartile.
#' 
#' Tukey (1977) also introduced the term Hinges and sorted the values in a W shape, where the bottom parts of the W are then 
#' the hinges.
#' 
#' There are quite a few different methods to determine the quartiles. This function has 12 different ones. See the details 
#' for a description.
#' 
#' @details
#' First the index of \eqn{Q_1} and \eqn{Q_3} needs to be determined. This depends on the method used 
#' and the sample size. The table below shows the calculation for the indexes (\eqn{iQ_1} and \eqn{iQ_3}). 
#' the 'm' is determined by:
#' \deqn{m = n \text{ mod } 4}
#' 
#' |method|m=0|m=1|m=2|m=3|
#' |---|----------|----------|----------|----------|
#' |Inclusive / Tukey / Vining| \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} |
#' |Exclusive / Joarder-Firozzaman | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} |  \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} |
#' |CDF / SAS 5 | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\lceil \frac{1}{4}\times n \rceil}, \eqn{\lceil \frac{3}{4}\times n \rceil} | \eqn{\lceil \frac{1}{4}\times n \rceil}, \eqn{\lceil \frac{3}{4}\times n \rceil} | \eqn{\lceil \frac{1}{4}\times n \rceil}, \eqn{\lceil \frac{3}{4}\times n \rceil} |
#' |Mendenhall-Sincich | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} | \eqn{\left[\frac{n + 1}{4}\right]}, \eqn{\lfloor \frac{3n + 3}{4}\rfloor} | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} |
#' |Lohninger | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} | \eqn{\left[\frac{n+1}{4}\right]}, \eqn{\left[\frac{3n+3}{4}\right]} |
#' |Hogg-Ledolter v1 | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 3}{4}} | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 1}{4}} |
#' |Hogg-Ledolter v2 | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 3}{4}} | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} | \eqn{\frac{n + 2}{4}}, \eqn{\frac{3n + 2}{4}} |
#' |Minitab / SAS 4 / Snedecor| \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} | \eqn{\frac{n + 1}{4}}, \eqn{\frac{3n + 3}{4}} |
#' |Excel| \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} | \eqn{\frac{n + 3}{4}}, \eqn{\frac{3n + 1}{4}} |
#' |SAS 1| \eqn{\frac{n}{4}}, \eqn{\frac{3n}{4}} | \eqn{\frac{n}{4}}, \eqn{\frac{3n}{4}} | \eqn{\frac{n}{4}}, \eqn{\frac{3n}{4}} | \eqn{\frac{n}{4}}, \eqn{\frac{3n}{4}} |
#' |SAS 2| \eqn{\left[\frac{n}{4}\right]}, \eqn{\left[\frac{3n}{4}\right]} | \eqn{\left[\frac{n}{4}\right]}, \eqn{\left[\frac{3n}{4}\right]} | \eqn{\lfloor\frac{n}{4}\rceil}, \eqn{\lfloor\frac{3n}{4}\rceil} | \eqn{\left[\frac{n}{4}\right]}, \eqn{\left[\frac{3n}{4}\right]} |
#' |SAS 3| \eqn{\lceil\frac{n}{4}\rceil}, \eqn{\lceil\frac{3n}{4}\rceil} | \eqn{\lceil\frac{n}{4}\rceil}, \eqn{\lceil\frac{3n}{4}\rceil} | \eqn{\lceil\frac{n}{4}\rceil}, \eqn{\lceil\frac{3n}{4}\rceil} | \eqn{\lceil\frac{n}{4}\rceil}, \eqn{\lceil\frac{3n}{4}\rceil} |
#' |Hyndman-Fan v8| \eqn{\frac{3n+5}{12}}, \eqn{\frac{9n+7}{12}} | \eqn{\frac{3n+5}{12}}, \eqn{\frac{9n+7}{12}} | \eqn{\frac{3n+5}{12}}, \eqn{\frac{9n+7}{12}} | \eqn{\frac{3n+5}{12}}, \eqn{\frac{9n+7}{12}} |
#' |Hyndman-Fan v9| \eqn{\frac{4n+7}{16}}, \eqn{\frac{12n+9}{16}} | \eqn{\frac{4n+7}{16}}, \eqn{\frac{12n+9}{16}} | \eqn{\frac{4n+7}{16}}, \eqn{\frac{12n+9}{16}} | \eqn{\frac{4n+7}{16}}, \eqn{\frac{12n+9}{16}} |
#' 
#' Once the index is known, the quartiles can be determined using the method shown in the table below using:
#' \deqn{fr = iHQ_x - iQ_x}
#' 
#' |method|fr=0|fr=0.25|fr=0.5|fr=0.75|
#' |---|----------|----------|----------|----------|
#' |\eqn{Q_x}| \eqn{x_{iLQ_x}} | \eqn{\frac{x_{iHQ_x + 3\times x_{iLQ_x}}}{4}} | \eqn{\frac{x_{iHQ_x + x_{iLQ_x}}}{2}} | \eqn{\frac{3\times x_{iHQ_x + x_{iLQ_x}}}{4}} |
#' 
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
#' The names of the different methods were adapted from Langford (2006). McAlister (1897, p. 374) uses 
#' the terms higher and lower quartile, while Galton (1881, p. 245) uses upper and lower.
#' 
#' The final formula for \eqn{Q_x} is based on linear interpolation. 
#' \deqn{Q_x = \frac{iQ_x - iLQ_x}{iHQ_x - iLQ_x}\times\left(x_{iHQ_x} - x_{iLQ_x}\right) + x_{iLQ_x}}
#' Since \eqn{iHQ_x - iLQ_x = 1} this can be further simplified to:
#' \deqn{Q_x = \left(iQ_x - iLQ_x\right)\times\left(x_{iHQ_x} - x_{iLQ_x}\right) + x_{iLQ_x}}
#' 
#' If the data is not numeric the interpolation will not be used and the two categories the quartile fell between 
#' is shown.
#' 
#' **Inclusive method vs. Exclusive**
#' 
#' For both methods divides the data into two, using the median as the cut-off. 
#' The quartiles are then the medians of each of the two parts. If the sample size is an even number.
#' The median will fall between two indexes, i.e. between \eqn{\frac{n}{2}} and \eqn{\frac{n + 2}{2}}.
#' The first half will therefor have the indexes 1 to \eqn{\frac{n}{2}}. The median of this will be:
#' \deqn{\frac{\frac{n}{2} + 1}{2} = \frac{\frac{n+2}{2}}{2}=\frac{n+2}{2}}.
#' The second half will start at \eqn{\frac{n + 2}{2}} and end at \eqn{n}. The number of scores in the 
#' upper half is then:
#' \deqn{n - \frac{n + 2}{2} + 1 = \frac{n - 2}{2} + 1 = \frac{n}{2}}
#' The same as in the lower half, as it should be. The median of the upper half will then be:
#' \deqn{\frac{n + 2}{2} + \frac{\frac{n}{2} + 1}{2} - 1 = \frac{n + 2 + \frac{n}{2} + 1 - 2}{2} = \frac{\frac{3n + 2}{2}}{2} = \frac{3n + 2}{4}}
#' 
#' However, this creates a problem if the median is a specific value (when n is odd). We could 
#' either exclude the median itself then in each of the two parts, or include it. Hence, the names 
#' 'inclusive' and 'exclusive.
#' 
#' If we include it, the number of values in each half will be:
#' \deqn{\frac{n + 1}{2}}
#' The median of the first half is then:
#' \deqn{\frac{\frac{n + 1}{2} + 1}{2} = \frac{\frac{n + 3}{2}}{2}= \frac{n + 3}{4}}
#' And for the second half:
#' \deqn{\frac{n + 1}{2} + \frac{\frac{n + 1}{2} + 1}{2} - 1 = \frac{n + 1 + \frac{n + 1}{2} + 1 - 2}{2} = \frac{\frac{3n + 1}{2}}{2} = \frac{3n + 1}{4}}
#' 
#' If we exclude it, the number of values in each half will be:
#' \deqn{\frac{n - 1}{2}}
#' The median of the first half is then:
#' \deqn{\frac{\frac{n - 1}{2} + 1}{2} = \frac{\frac{n + 1}{2}}{2}= \frac{n + 1}{4}}
#' And for the second half:
#' \deqn{\frac{n + 3}{2} + \frac{\frac{n + 1}{2} + 1}{2} - 1 = \frac{n + 3 + \frac{n + 1}{2} + 1 - 2}{2} = \frac{\frac{3n + 3}{2}}{2} = \frac{3n + 3}{4}}
#' 
#' The inclusive method can be found in Tukey (1977), but also Siegel and Morgan (1996, p. 77) and Vining (1998, p. 44).
#' While the exclusive method could be found in Moore and McCabe (1989, p. 33) or Joarder and Firozzaman (2001, p. 88).
#' 
#' **SAS versions**
#' 
#' Another approach to the issue is simply using \eqn{n\times p} where \eqn{p} is the percentile, which 
#' for the first quartile is then 0.25, and for the third 0.75.
#' 
#' Simply using \eqn{iQ_x = n\times p} is what **SAS method 1 does**. 
#' \deqn{iQ_1 = \frac{1}{4}\times n = \frac{n}{4}}
#' \deqn{iQ_3 = \frac{3}{4}\times n = \frac{3n}{4}}
#' 
#' **SAS method 3** simply rounds each of the results up to the nearest integer:
#' \deqn{iQ_1 = \lceil\frac{1}{4}\times n\rceil = \lceil\frac{n}{4}\rceil}
#' \deqn{iQ_3 = \lceil\frac{3}{4}\times n\rceil = \lceil\frac{3n}{4}\rceil}
#' 
#' **SAS method 2** is very similar, except if \eqn{n\times p} ends with 0.5. It then suggests to round to 
#' the nearest even integer. So for example 3.5 gets rounded to 4, but 2.5 gets rounded to 2.
#' Note that we only get \eqn{n\times p} ending with .5 if n is a multiple of 4 + 2, i.e. \eqn{n \text{mod} 4 = 2}
#' 
#' **SAS method 4** first adds 1 to the sample size, so uses \eqn{\left(n + 1\right)\times p}. This leads for 
#' the two quartiles to:
#' \deqn{iQ_1 = \frac{1}{4}\times\left(n + 1\right) = \frac{n + 1}{4}}
#' \deqn{iQ_3 = \frac{3}{4}\times\left(n + 1\right) = \frac{3n + 3}{4}}
#' 
#' The SAS method 4 is also used by **Minitab** and can also be found in Snedecor (1940, p. 43).
#' 
#' **SAS method 5** is the same as SAS method 3, except \eqn{n\times p} is an integer. This will 
#' only occur if \eqn{n} is a multiple of 4. In that case it simply takes the average of the index and the 
#' next value i.e. \eqn{n\times p + 1}. For the first quartile this leads to:
#' \deqn{iQ_1 = \frac{\frac{n}{4} + \frac{n}{4} + 1}{2} = \frac{\frac{2n + 4}{4}}{2}= \frac{\frac{n + 2}{2}}{2} = \frac{n + 2}{4}}
#' And for the third:
#' \deqn{iQ_3 = \frac{\frac{3n}{4} + \frac{3n}{4} + 1}{2} = \frac{\frac{6n + 4}{4}}{2}= \frac{\frac{3n + 2}{2}}{2} = \frac{3n + 2}{4}}
#' 
#' The methods used by SAS are found in the procedures guide (SAS, 2006, p. 626).
#' 
#' **Mendenhall and Sincich**
#' 
#' Similar as SAS method 4, however Mendenhall and Sincich (1992, p. 35) round the values. 
#' For the lower (first) quartile use:
#' \deqn{iQ_1 = \left[\frac{1}{4}\times \left(n + 1\right)\right] = \left[\frac{n + 1}{4}\right]}
#' However, if the third quartile index ends with 0.5 round down. The index will become ending with 
#' a 0.5 if \eqn{n \text{ mod } 4 = 1}, since then \eqn{n + 1} will have a modulo 2, which divided by 4 
#' gives the .5. So the third quartile is then:
#' \deqn{iQ_3 = \lfloor \frac{3}{4}\times\left(n + 1\right)\rfloor = \lfloor \frac{3\times n + 3}{4}\rfloor}
#' Otherwise:
#' \deqn{iQ_3 = \left[\frac{3}{4}\times \left(n + 1\right)\right] = \left[\frac{3\times n + 3}{4}\right]}
#' 
#' **Lohninger**
#' 
#' Lohninger (n.d.) does the same as Mendenhall and Sincich, but simply always uses 'regular' roundings.
#' 
#' \deqn{\left[iQ_1 = \frac{1}{4}\times \left(n + 1\right)\right] = \left[\frac{n + 1}{4}\right]}
#' and for the third quartile
#' \deqn{\left[iQ_3 = \frac{3}{4}\times \left(n + 1\right)\right] = \left[\frac{3\times n + 3}{4}\right]}
#' 
#' **Hogg and Ledolter**
#' 
#' Hogg and Ledolter (1992, p. 21) use that for the median \eqn{\frac{n}{2} + \frac{1}{2}} is used.
#' So use for the first quartile the same approach:
#' \deqn{iQ_1 = \frac{n}{4} + \frac{1}{2} = \frac{n + 2}{4}}
#' For the third:
#' \deqn{iQ_3 = \frac{3n}{4} + \frac{1}{2} = \frac{3n + 2}{4}}
#' If however this falls between two indexes, average the two they fall between.
#' This means we can't use the linear interpolation directly, since the two equations could end with 
#' .25 and .75 as well, and then still need to just determine the average of the two.
#' If it ends in .25 we simply add 1/4, and if it ends in .75 we subtract it.
#' So if \eqn{n \text{ mod } 4 = 1}:
#' \deqn{iQ_1 = \frac{n + 2}{4} + \frac{1}{4} = \frac{n + 3}{4}}
#' \deqn{iQ_3 = \frac{3n + 2}{4} + \frac{1}{4} = \frac{3n + 3}{4}}
#' And if it ends in .75, i.e. \eqn{n \text{ mod } 4 = 3}:
#' \deqn{iQ_1 = \frac{n + 2}{4} - \frac{1}{4} = \frac{n + 1}{4}}
#' \deqn{iQ_3 = \frac{3n + 2}{4} - \frac{1}{4} = \frac{3n + 1}{4}}
#' 
#' It could be that this method was also described in Hazen (1914), but I couldn't find the exact 
#' page number in that document.
#' 
#' **Excel**
#' 
#' It appears that MS Excel has a unique method. It can also be found in Freund and Perles (1987, p. 201):
#' \deqn{iQ_1 = \frac{1}{4}\times \left(n - 1\right) + 1 = \frac{n + 3}{4}}
#' \deqn{iQ_3 = \frac{3}{4}\times \left(n - 1\right) + 1 = \frac{3\times n + 1}{4}}
#' 
#' **Hyndman-Fan**
#' 
#' Hyndman and Fan (1996) discuss 9 different methods for quantiles, which can be used to determine quartiles.
#' The R stats library's quantile() function uses their numbering. 
#' HF-1 = SAS-3, HF-2 = CDF / SAS-5, HF-3 = SAS-2, HF-4 = SAS-1, HF-5 = Hogg-Ledolter v2, 
#' HF-6 = Minitab / SAS-4, and HF-7 = Excel.
#' 
#' This leaves their 8th and 9th definitions. For their 8th (method="hf8"):
#' \deqn{iQ_1 = \left(n + \frac{1}{3}\right)\times\frac{1}{4}+\frac{1}{3} = \frac{3n+5}{12}}
#' \deqn{iQ_3 = \left(n + \frac{1}{3}\right)\times\frac{3}{4}+\frac{1}{3} = \frac{9n+7}{12}}
#' 
#' For their 9th method (method="hf9"):
#' \deqn{iQ_1 = \left(n + \frac{1}{4}\right)\times\frac{1}{4}+\frac{3}{8} = \frac{4n+7}{16}}
#' \deqn{iQ_1 = \left(n + \frac{1}{4}\right)\times\frac{3}{4}+\frac{3}{8} = \frac{12n+9}{16}}
#' 
#' @examples 
#' ex1 = c(1, 2, 3, 4, 5, 6, 7)
#' me_quartiles(ex1, method="inclusive")
#' me_quartiles(ex1, method="exclusive")
#' me_quartiles(ex1, method="tukey")
#' me_quartiles(ex1, method="cdf")
#' me_quartiles(ex1, method="ms")
#' me_quartiles(ex1, method="lohninger")
#' me_quartiles(ex1, method="vining")
#' me_quartiles(ex1, method="jf")
#' me_quartiles(ex1, method="hl1")
#' me_quartiles(ex1, method="hl2")
#' me_quartiles(ex1, method="minitab")
#' me_quartiles(ex1, method="excel")
#' me_quartiles(ex1, method="sas1")
#' me_quartiles(ex1, method="sas2")
#' me_quartiles(ex1, method="sas3")
#' 
#' @seealso 
#' Quartiles are a special case of quantiles, see \code{\link{me_quantiles}} for more info on those.
#' 
#' The quartiles are used for different ranges (interquartile, semi-interquartile and mid-quartile). 
#' See \code{\link{me_quartile_range}} for more info on those.
#' 
#' @references 
#' Freund, J. E., & Perles, B. M. (1987). A new look at quartiles of ungrouped data. *The American Statistician, 41*(3), 200–203. https://doi.org/10.1080/00031305.1987.10475479
#' 
#' Galton, F. (1881). Report of the anthropometric committee. *Report of the British Association for the Advancement of Science, 51*, 225–272.
#' 
#' Hazen, A. (1914). Storage to be provided in impounding municipal water supply. *Transactions of the American Society of Civil Engineers, 77*(1), 1539–1640. https://doi.org/10.1061/taceat.0002563
#' 
#' Hogg, R. V., & Ledolter, J. (1992). *Applied statistics for engineers and physical scientists* (2nd int.). Macmillan.
#' 
#' Hyndman, R. J., & Fan, Y. (1996). *Sample quantiles in statistical packages. The American Statistician, 50*(4), 361–365. https://doi.org/10.2307/2684934
#' 
#' Joarder, A. H., & Firozzaman, M. (2001). *Quartiles for discrete data. Teaching Statistics, 23*(3), 86–89. https://doi.org/10.1111/1467-9639.00063
#' 
#' Langford, E. (2006). Quartiles in elementary statistics. *Journal of Statistics Education, 14*(3), 1–17. https://doi.org/10.1080/10691898.2006.11910589
#' 
#' Lohninger, H. (n.d.). Quartile. Fundamentals of Statistics. Retrieved April 7, 2023, from http://www.statistics4u.com/fundstat_eng/cc_quartile.html
#' 
#' McAlister, D. (1879). The law of the geometric mean. *Proceedings of the Royal Society of London, 29*(196–199), 367–376. https://doi.org/10.1098/rspl.1879.0061
#' 
#' Mendenhall, W., & Sincich, T. (1992). S*tatistics for engineering and the sciences* (3rd ed.). Dellen Publishing Company.
#' 
#' Moore, D. S., & McCabe, G. P. (1989). *Introduction to the practice of statistics*. W.H. Freeman.
#' 
#' Parzen, E. (1979). Nonparametric statistical data modeling. *Journal of the American Statistical Association, 74*(365), 105–121. https://doi.org/10.1080/01621459.1979.10481621
#' 
#' SAS. (1990). *SAS procedures guide: Version 6* (3rd ed.). SAS Institute.
#' 
#' Siegel, A. F., & Morgan, C. J. (1996). *Statistics and data analysis: An introduction* (2nd ed.). J. Wiley.
#' 
#' Snedecor, G. W. (1940). *Statistical methods applied to experiments in agriculture and biology* (3rd ed.). The Iowa State College Press.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' Vining, G. G. (1998). *Statistical methods for engineers*. Duxbury Press.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
me_quartiles <- function(data, 
                         method=c("cdf", "inclusive", "exclusive", "tukey", "ms", "lohninger", 
                                  "vining", "jf", "hl1", "hl2", "minitab", "excel", "sas1", 
                                  "sas2", "sas3", "hf8", "hf9")){
  
  if (length(method)>1) {
    #default is recommended method by Langford (2006)
    method="cdf"
  }
  
  sData = sort(data)
  n = length(sData)
  m = n %% 4
  
  
  #NOTE: R uses bankings rounding (i.e. rounds .5 to the nearest even number)
  #use as.integer(x + 0.5) to us 'regular' rounding
  
  if (method=="inclusive" || method=="tukey" || method=="vining") {
    #See also Siegel and Morgan (1996, p. 77), or Tukey's hinges (Tukey, 1977, p. 33), or Vining (1998, p. 44) .
    if (m==0 || m ==2){
      iQ1 = (n + 2)/4
      iQ3 = (3*n + 2)/4
    }
    else{
      iQ1 = (n + 3)/4
      iQ3 = (3*n + 1)/4
    }
  }
  
  else if (method=="exclusive" || method=="jf") {
    #also referred to as Moore and McCabe (1989, p. 33), or Joarder and Firozzaman (2001, p. 88) 
    if (m==0 || m ==2){
      iQ1 = (n + 2)/4
      iQ3 = (3*n + 2)/4
    }
    else{
      iQ1 = (n + 1)/4
      iQ3 = (3*n + 3)/4
    }
  }
  
  else if (method=="cdf" || method=="sas5" || method=="hf2"){
    if (m==0){
      iQ1 = (n + 2)/4
      iQ3 = (3*n + 2)/4
    }
    else{
      iQ1 = ceiling(n/4)
      iQ3 = ceiling(3*n/4)
    }
  }
  
  else if (method=="ms"){
    #Mendenhall and Sincich (1992, p. 35)
    iQ1 = as.integer((n + 1)/4+0.5)
    if (m == 1){
      iQ3 = floor(3*(n + 1)/4)
    }
    else{
      iQ3 = as.integer(3*(n + 1)/4 + 0.5)
    }
  }
  
  else if (method=="lohninger"){
    # Lohninger (n.d.)
    iQ1 = as.integer((n + 1)/4+0.5)
    iQ3 = as.integer(3*(n + 1)/4+0.5)
  }
  
  else if (method=="hl1"){
    #Hogg and Ledolter (1992, p. 21)
    if (m==0 || m==2){
      iQ1 = (n+2)/4
      iQ3 = (3*n+2)/4
    }
    else if (m==1){
      iQ1 = (n+3)/4
      iQ3 = (3*n+3)/4
    }
    else {
      iQ1 = (n+1)/4
      iQ3 = (3*n+1)/4
    }
  }
  else if (method=="hl2" || method=="hf5"){
    #Hogg and Ledolter (1992, p. 21)
    iQ1 = (n+2)/4
    iQ3 = (3*n+2)/4
  }
  
  else if (method=="minitab" || method =="sas4" || method=="hf6"){
    #See also Snedecor (1940, p. 43)
    iQ1 = (n+1)/4
    iQ3 = (n+1)*3/4
  }
  
  else if (method=="excel" || method=="hf7"){
    #See also Freund and Perles (1987, p. 201)
    iQ1 = 1 + (n-1)/4
    iQ3 = 1 + (n-1)*3/4
  }
  
  else if (method=="sas1" || method=="hf4"){
    #See also Parzen (1979, p. ???)
    iQ1 = 1/4*n
    iQ3 = 3/4*n
  }
  
  else if (method=="sas2" || method=="hf3"){
    iQ1 = round(n*1/4)
    iQ3 = round(n*3/4)
  }
  
  else if (method=="sas3" || method=="hf1"){
    #See also 
    iQ1 = ceiling(1/4*n)
    iQ3 = ceiling(3/4*n)
  }
  
  else if (method=="hf8"){
    #See also Hyndman and Fan (1996, p. 363)
    iQ1 = (n + 1/3)*1/4 + 1/3
    iQ3 = (n + 1/3)*3/4 + 1/3
  }
  
  else if (method=="hf9"){
    #See also Hyndman and Fan (1996, p. 364)
    iQ1 = (n + 1/4)*1/4 + 3/8
    iQ3 = (n + 1/4)*3/4 + 3/8
  }
  
  #make sure integer, if indeed integer
  q1Int = FALSE
  if (iQ1 == round(iQ1)) {
    iQ1 = round(iQ1)
    q1Int = TRUE}
  q3Int = FALSE
  if (iQ3 == round(iQ3)) {
    iQ3 = round(iQ3)
    q3Int = TRUE}
  
  #find the quartiles
  numbers = is.numeric(sData[1]) 
  
  if (q1Int){
    Q1 = sData[iQ1]
  }
  else {
    iQ1low = floor(iQ1)
    iQ1high = ceiling(iQ1)
    
    if (sData[iQ1low]==sData[iQ1high]) {
      Q1 = sData[iQ1low]
    }
    else{
      if (numbers) {
        r = sData[iQ1high] - sData[iQ1low]
        Q1 = sData[iQ1low] + (iQ1 - iQ1low)/(iQ1high - iQ1low)*r
      }
      else{
        Q1 = paste0("between ", sData[iQ1low], " and ", sData[iQ1high])
      }
    }
  }
  
  if (q3Int){
    Q3 = sData[iQ3]
  }
  else {
    iQ3low = floor(iQ3)
    iQ3high = ceiling(iQ3)
    
    if (sData[iQ3low]==sData[iQ3high]) {
      Q3 = sData[iQ3low]
    }
    else{
      if (numbers) {
        r = sData[iQ3high] - sData[iQ3low]
        Q3 = sData[iQ3low] + (iQ3 - iQ3low)/(iQ3high - iQ3low)*r
      }
      else{
        Q3 = paste0("between ", sData[iQ3low], " and ", sData[iQ3high])
      }
    }
  }
  
  results = data.frame(Q1, Q3)
  
  return(results)
}

