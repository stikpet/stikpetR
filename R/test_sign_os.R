#' one-sample sign test
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized median, otherwise the midrange will be used
#' @return dataframe with mu, the significance (p-value) and test used
#' 
#' @examples 
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' ts_sign_os(data)
#' ts_sign_os(data, mu = 2)
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
#' The paired version for example was already described by Arbuthnott (1710)
#' 
#' **Alternatives**
#' 
#' *library(BSDA)* 
#' 
#' SIGN.test()
#' 
#' *library(DescTools)* 
#' 
#' SignTest()
#' 
#' *library(nonpar)*
#' 
#' signtest(data, m=2)
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Arbuthnott, J. (1710). An argument for divine providence, taken from the constant regularity observ’d in the births of both sexes. *Philosophical Transactions of the Royal Society of London, 27*(328), 186–190. https://doi.org/10.1098/rstl.1710.0011
#' 
#' Stewart, W. M. (1941). A note on the power of the sign test. *The Annals of Mathematical Statistics, 12*(2), 236–239. https://doi.org/10.1214/aoms/1177731755
#'  
#' @export
ts_sign_os <- function(data, mu = NULL){
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
  
  return(testResults)
  
}

