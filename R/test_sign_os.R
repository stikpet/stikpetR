#' one-sample sign test
#' 
#' @description 
#' This function will perform one-sample sign test.
#' 
#' @param data A vector or dataframe
#' @param levels optional vector with levels in order
#' @param mu optional hypothesized median, otherwise the midrange will be used
#' 
#' @returns 
#' Dataframe with:
#' \item{mu}{the mean tested}
#' \item{p-value}{he significance (p-value)}
#' \item{test}{a description of the test used}
#' 
#' @details 
#' The test statistic is calculated using (Stewart, 1941, p. 236):
#' \deqn{p = 2\times B\left(n, \text{min}\left(n_+, n_-\right), \frac{1}{2}\right)}
#' 
#' *Symbols used:*
#' 	\itemize{
#' 	\item \eqn{B\left(\dots\right)} is the binomial cumulative distribution function
#' 	\item \eqn{n} is the number of cases
#' 	\item \eqn{n_+} is the number of cases above the hypothesized median
#' 	\item \eqn{n_-} is the number of cases below the hypothesized median
#' 	\item \eqn{min} is the minimum value of the two values
#' 	}
#' 
#' The test is described in Stewart (1941), although there are earlier uses. 
#' 
#' The paired version for example was already described by Arbuthnott (1710)
#' 
#' @section Alternatives:
#' 
#' The library *BSDA* has a **SIGN.test()** function
#' 
#' The library *DescTools* has a **SignTest()** function
#' 
#' The library *nonpar* has a **signtest()** function
#' 
#' @references 
#' Arbuthnott, J. (1710). An argument for divine providence, taken from the constant regularity observ’d in the births of both sexes. *Philosophical Transactions of the Royal Society of London, 27*(328), 186–190. https://doi.org/10.1098/rstl.1710.0011
#' 
#' Stewart, W. M. (1941). A note on the power of the sign test. *The Annals of Mathematical Statistics, 12*(2), 236–239. https://doi.org/10.1214/aoms/1177731755
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: Text dataframe
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' ts_sign_os(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' ts_sign_os(ex2) 
#' 
#' @export
ts_sign_os <- function(data, levels=NULL, mu = NULL){
  
  if (!is.null(levels)){
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    data = as.numeric(myFieldOrd)
  }
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }
  
  #Determine count of cases below hypothesized median
  group1 = data[data<mu]
  group2 = data[data>mu]
  n1 = length(group1)
  n2 = length(group2)
  
  #Select the lowest of the two
  myMin = min(n1,n2)
  
  #Determine total number of cases (unequal to hyp. median)
  n = n1+n2
  
  #Determine the significance using binomial test
  pValue = 2*pbinom(myMin, n,0.5)
  
  testUsed = "one-sample sign test"
  testResults <- data.frame(mu, pValue, testUsed)
  colnames(testResults)<-c("mu", "p-value", "test")
  
  return(testResults)
  
}

