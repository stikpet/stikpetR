#' Frequency Table
#' 
#' @param data A vector with the data
#' @returns 
#' Dataframe with the folowing columns:
#' \item{data}{the categories}
#' \item{Freq}{the absolute count}
#' \item{Percent}{the percentage based on the total including missing values}
#' \item{ValidPercent}{the percentage based on the total excluding missing values, only if missing values are present}
#' \item{CumulativePercent}{the cumulative percentages}
#' 
#' @description 
#' A frequency table is defined as "a table showing (1) all of the values for a variable in a dataset, 
#' and (2) the frequency of each of those responses. Some frequency tables also show a cumulative frequency and 
#' proportions of responses" (Warne, 2017, p. 512). 
#' 
#' A frequency table can help to get impression of your survey data of a binary, nominal, or ordinal variable. 
#' It could also help with a scale variable, provided there are not too many options. If, for example, you have asked for 
#' age, a list going from 1 to 90 with different ages and frequencies, will probably not be so helpful.
#' 
#' If you have many options in the scale variable, the data is often binned (e.g. 0 < 10, 10 < 20, etc.), which creates 
#' then an ordinal variable, of which a frequency table can then be helpful. See binning for more information on this.
#' A frequency table can show different types of frequencies. Various options are discussed in the details.
#' 
#' A YouTube video with explanation on this test is available [here](https://youtu.be/DPmwWxYYCp4)
#'  
#' @details 
#' The column **Frequency** shows how many respondents answered each option. We can tell that 100 people in this survey 
#' chose the option 'very scientific'. This is also known as the **absolute frequency** and defined as “the number of 
#' occurrences of a particular phenomenon” (Zedeck, “Frequency”, 2014, p. 144).
#' 
#' The **Percent** column shows the percentages, based on the grand total, so including the missing values. 
#' Percentages can be defined as “a way of expressing ratios in terms of whole numbers. A ratio or fraction is converted to 
#' a percentage by multiplying by 100 and appending a "percentage sign" %” (Weisstein, 2002, p. 2200).
#' 
#' The **Valid Percent** shows the percentage, based on the valid total, so excluding the missing values. 
#' Most often the ‘Percent’ shown in reports are actually Valid Percent, but the word ‘Valid’ is then simply left out.
#' 
#' Percentages show the number of cases that could be expected if there would be 100 cases in total, hence per-cent which means 
#' 'per 100'. If your sample size is very small, be careful about using percentages. If it is less than 100, it means that you 
#' are 'blowing up' your differences, while percentages are more commonly used to 'scale down'.
#' 
#' The term **relative frequency** is also sometimes used. This is the frequency divided by the total number of cases. 
#' Note that this should then always produce a decimal value between 0 and 1 (inclusive). 
#' Multiply this by 100 and you get the percentage, multiply it by 1000 and you get permille (‰), multiply it by 360 and 
#' you get the degrees of a circle, etc.
#' 
#' In general the formula for a percentage is:
#' \(PR_i = \frac{F_i}{n}\times 100 \)
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{PR_i} the percentage of category i
#' \item \eqn{F_i} the (absolute) frequency of category i
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies (either including or excluding the missing values)
#' }
#' 
#' The **cumulative frequency** (not shown in table) can be defined as: “the total (absolute) frequency up to the upper 
#' boundary of that class” (Kenney, 1939, p. 16). This would only be useful if there is an order to the categories, 
#' so we can say that for example 299 respondents found accounting pretty scientific or even more. Which is why these 
#' cumulative frequencies will not have a meaningful interpretation for a nominal variable (e.g. 28 students study 
#' business or less?).
#' 
#' The **Cumulative Percent** is the running total of the Valid Percent, it is the addition of all previous and the current 
#' category’s valid percentages.
#' 
#' The cumulative frequency can be calculated using:
#' \deqn{CF_i = \sum_{j=1}^i F_j}
#' Or using recursion:
#' \deqn{CF_i = F_i + CF_{i-1}}
#' For the cumulative percent the same formulas as for cumulative frequency can be used, but replacing \eqn{F_i} with \eqn{PR_i}. 
#' It can also be determined using the cumulative frequency:
#' \deqn{CPR_i = \frac{CF_i}{n}}
#' 
#' When the categories are ranges of values (bins), the frequency density could become helpful. 
#' It can be defined as: “the number of occurrences of an event divided by the bin size…” (Zedeck, 2014, pp. 144–145).
#' See the binned tables for more information about this.
#' 
#' @examples  
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", NA, "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' tab_frequency(data)
#' 
#' @references 
#' Kenney, J. F. (1939). *Mathematics of statistics; Part one*. Chapman & Hall.
#' 
#' Warne, R. T. (2017). *Statistics for the social sciences: A general linear model approach*. Cambridge University Press.
#' 
#' Weisstein, E. W. (2002). *CRC concise encyclopedia of mathematics* (2nd ed.). Chapman & Hall/CRC.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
tab_frequency <- function(data){
  #counts without missing values
  freq = table(data)
  #counts with missing values
  freq2 = table(data, exclude=NULL)
  
  #Missing values present?
  noMissing <- sum(freq) == sum(freq2)

  #Proportions without missing values
  prop = prop.table(freq)  
  #Proportions with missing values
  prop2 = prop.table(freq2)
  
  #Cumulative Frequency
  cumf  = cumsum(freq)
  #Cumulative Relative Frequency
  cumrf = cumsum(prop)
  
  
  f = as.data.frame(freq2,stringsAsFactors=FALSE)
  rf = as.data.frame(prop2,stringsAsFactors=FALSE)
  vp = as.data.frame(prop,stringsAsFactors=FALSE)
  crf = as.data.frame(cumrf,stringsAsFactors=FALSE)
  myTable= f
  myTable$Percent = rf$Freq*100
  
  if (noMissing) {
    myTable$CumulativePercent = crf$cumrf*100
    myTable = rbind(myTable,c("TOTAL",sum(myTable$Freq),sum(myTable$Percent),NA))
  }
  else{
    vp2 = rbind(vp,c(NA,NA))
    crf2 = rbind(crf,c(NA))
    myTable$ValidPercent = vp2$Freq*100
    myTable$CumulativePercent = crf2$cumrf*100
    myTable = rbind(myTable,c("TOTAL",sum(myTable$Freq),sum(myTable$Percent),sum(myTable$ValidPercent, na.rm = TRUE),NA))
  }
  
  return(myTable)
}