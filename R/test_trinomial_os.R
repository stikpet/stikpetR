#' One-Sample Trinomial Test
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized median, otherwise the midrange will be used
#' @return dataframe with mu, the number of scores above mu, below and tied, the significance (p-value) and test used
#' 
#' @examples 
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' ts_trinomial_os(data)
#' ts_trinomial_os(data, mu = 2)
#' 
#' @details 
#' The p-value is calculated using (Bian et al., 2009, p. 6):
#' \deqn{p = 2\times \sum_{i=n_d}^n \sum_{j=0}^{\lfloor \frac{n - i}{2} \rfloor} \text{tri}\left(\left(j, j+i, n - i\right), \left(p_{pos}, p_{neg}, p_0\right) \right)}
#' 
#' With:
#' \deqn{p_0 = \frac{n_0}{n}}
#' \deqn{p_{pos} = p_{neg} = \frac{1 - p_0}{n}}
#' \deqn{\left|n_{pos} - n_{neg}\right|}
#' 
#' *Symbols used:*
#' 	\itemize{
#' 	\item \eqn{n_0} the number of scores equal to the hypothesized median
#' 	\item \eqn{n_{pos}} the number of scores above the hypothesized median
#' 	\item \eqn{n_{neg}} the number of scores below the hypothesized median
#' 	\item \eqn{p_0} the probability of the a score in the sample being equal to the hypothesized median
#' 	\item \eqn{p_{pos}} the population proportion of a score being above the hypothesized median
#' 	\item \eqn{p_{neg}} the population proportion of a score being below the hypothesized median
#' 	\item \eqn{\text{tri}\left(…,… \right)} the trinomial probability mass function
#' 	}
#' 
#' The paired version of the test is described in Bian et al. (1941), while Zaiontz (n.d.) mentions it can 
#' also be used for one-sample situations.
#' 
#' **Alternatives**
#' 
#' I'm not aware of any other library that has a similar function.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Bian, G., McAleer, M., & Wong, W.-K. (2009). A trinomial test for paired data when there are many ties. *SSRN Electronic Journal*. https://doi.org/10.2139/ssrn.1410589
#' 
#' Zaiontz, C. (n.d.). Trinomial test. Real Statistics Using Excel. Retrieved March 2, 2023, from https://real-statistics.com/non-parametric-tests/trinomial-test/
#'  
#' @export
ts_trinomial_os <- function(data, mu=NULL){
  testUsed = "one-sample trinomial test"
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }

  pos = sum(data>mu)
  neg = sum(data<mu)
  ties = sum(data==mu)
  
  n = pos + neg + ties
  
  p0 = ties/n
  p1 = (1 - p0)/2
  k = abs(pos - neg)
  
  sig=0
  for (z in k:n) {
    for (i in 0:floor((n - z)/2)) {
      sig = sig + dmultinom(c(i, i+z, n - i - (i+z)), prob = c(p1, p1, p0))
    }
  }
  
  pValue = sig*2  
  
  results = data.frame(mu, pos, neg, ties, pValue, testUsed)
  
  return(results)
  
}
