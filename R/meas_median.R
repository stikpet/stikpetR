#' Median
#' 
#' @description 
#' Function to determine the median of a set of data. The median can be defined as "the middle value in a distribution, below and above which lie values with equal total frequencies or probabilities" (Porkess, 1991, p. 134). This means that 50% of the respondents scored equal or higher to the median, and also 50% of the respondents scored lower or equal.
#' 
#' @param data vector with the data
#' @param levels optional list to indicate what values represent
#' @param tieBreaker optional which to return if median falls between two values. Either `"between"` (default), `"low"`, or `"high"`
#' 
#' @return the median
#' 
#' @details 
#' The formula that is used, assuming the data has been sorted, is:
#' \deqn{\tilde{x} = \begin{cases} x_{MI} & \text{ if } MI= \left \lfloor MI \right \rfloor\\ \frac{x_{MI-0.5} + x_{MI+0.5}}{2} & \text{ if } MI\neq \left \lfloor MI \right \rfloor \end{cases}}
#' With:
#' \deqn{MI = \frac{n + 1}{2}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} the sample size
#' \item \eqn{x_i} the i-th score of X, assuming X has been sorted.
#' \item \eqn{MI} the index of the median
#' \item \eqn{\tilde{x}} the median
#' }
#' 
#' If the number of scores is an odd number, and the median falls between two categories. With the 
#' *tieBreaker* it can then be set to return the lower value (low), upper (high), or average (between).
#' 
#' Some old references to the median are Pacioli (1523) in Italian, Cournot (1843, p. 120) in French, 
#' and Galton (1881, p. 246) in English.
#' 
#' @section Alternatives:
#' 
#' R's *stats* library has a function *median*, and also *DescTools*, *missMethods*, and probably 
#' more. Note that many of these do not make it possible to return the "between" score if the data 
#' are ordered factors.
#' 
#' @references 
#' Cournot, A. A. (1843). *Exposition de la théorie des chances et des probabilités*. L. Hachette.
#' 
#' Galton, F. (1881). Report of the anthropometric committee. *Report of the British Association for the Advancement of Science, 51*, 225–272.
#' 
#' Pacioli, L. (1523). *Summa de arithmetica geometria proportioni: Et proportionalita*. Paganino de Paganini.
#' 
#' Porkess, R. (1991). *The HarperCollins dictionary of statistics*. HarperPerennial.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Text dataframe
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' me_median(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' me_median(ex2)
#' 
#' #Example 3: Text data with between median
#' ex3 = c("a", "b", "f", "d", "e", "c")
#' order = c("a", "b", "c", "d", "e", "f")
#' me_median(ex3, levels=order)
#' me_median(ex3, levels=order, tieBreaker="low")
#' me_median(ex3, levels=order, tieBreaker="high")
#' 
#' #Example 4: Numeric data with between median
#' ex4 = c(1, 2, 3, 4, 5, 6)
#' me_median(ex4)
#' me_median(ex4, tieBreaker="low")
#' me_median(ex4, tieBreaker="high")
#' 
#' @export
me_median <- function(data, levels=NULL, tieBreaker=c("between", "low", "high")){
  
  if (length(tieBreaker)>1) {tieBreaker = "between"}
  
  if (!is.null(levels)){
    dataN = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(dataN)
  }
  else{dataN = data}
  
  sData = sort(na.omit(dataN))
  n = length(dataN)
  medIndex = (n + 1)/2
  
  if (medIndex == round(medIndex)) {med = sData[medIndex]}
  else{
    medLow = sData[medIndex - 0.5]
    medHigh = sData[medIndex + 0.5]
    
    if (tieBreaker=="between") {
      medN = (medLow + medHigh)/2
    }
    else if (tieBreaker=="low"){
      medN = medLow
    }
    else if (tieBreaker=="high"){
      medN = medHigh
    }
  }
  
  if (is.null(levels)){
    med = medN
  }
  else{
      if (medN == round(medN)) {
        med = levels[medN]
      }
    else{
      med = paste0("between ", levels[medLow], " and ", levels[medHigh])
    }
  }
  return(med)
}