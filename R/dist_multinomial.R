#' Multinomial Probability Mass Function
#' @description
#' This is a function for the multinomial probability. It returns the probability of a distribution as given in F for a sample size of sum of F, where the probability for each category is given as in P. It is a generalization of the binomial distribution.
#' 
#' The distribution is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Distributions/Multinomial.html)
#' 
#' 
#' @param F list with the observed counts
#' @param P list with the probabilities for each category
#' @param method optional the calculation method to use. Either "loggamma" (default), "factorial", "gamma", "mprob".
#' 
#' 
#' @returns
#' A float with the requested probability
#' 
#' 
#' @details 
#' If *method=factorial* the following formula is used:    
#' \deqn{mpmf\left(F, P\right) = \frac{n!}{\prod_{i=1}^{k} \left(F_i!\right)} \times \prod_{i=1}^{k} P_i^{F_i}}
#' This formula was most likely already used by for example Edgeworth (1905), but can for example also be found in Berry and Mielke (1995, p. 769)
#' 
#' If *method=gamma*:
#' \deqn{mpmf\left(F, P\right) = \frac{\Gamma\left(1+n\right)}{\prod_{i=1}^n \Gamma\left(1+F_i\right)} \times \prod_{i=1}^{k} P_i^{F_i}} 
#' If *method=loggamma*:
#' \deqn{mpmf\left(F, P\right) = e^{\ln\left(mpmf\left(F, P\right)\right)}}
#' \deqn{\ln\left(mpmf\left(F, P\right)\right) = \ln\left(\Gamma\left(n + 1\right)\right) + \sum_{i=1}^k F_i\times\ln\left(P_i\right) - \ln\left(\Gamma\left(F_i + 1\right)\right)}
#' This formula can for example be found in Arnold (2018).
#' 
#' If *method=mprob* the algorithm from García-Pérez (1999) is used:
#' \enumerate{
#' \item Determine \eqn{F^*}, the counts in descending order, and move the elements in \eqn{P} accordingly creating \eqn{P^*}.
#' \item Set \eqn{pmf = 1}, \eqn{t=P_1^*}, \eqn{i=2}, \eqn{x=0}, and \eqn{m=F_1^*}
#' \item Set \eqn{l = F_i^*}. For \eqn{r=1} to \eqn{l} do:
#' \itemize{
#' \item update \eqn{x = x + 1}
#' \item if \eqn{x > F_1^*} then set \eqn{t = 1} (else nothing)
#' \item update \eqn{pmf = pmf\times t\times P_i^*\times \frac{r + m}{r}}
#' }
#' \item If \eqn{i=k}, then go to step 5, otherwise update \eqn{i = i + 1}, \eqn{m = m + F_i^*} and go to step 3
#' \item If \eqn{x < F_1^*} then for \eqn{r=x + 1} to \eqn{F_1^*} update \eqn{pmf = pmf\times P_1^*}
#' }
#' 
#' This distribution is used in a Multinomial Goodness-of-Fit Test. The stikpetR library has a function \code{\link{ts_multinomial_gof}} for this, but it uses the dmultinomial function from R.
#' 
#' 
#' @references 
#' Arnold, J. (2018, December 3). Maximum Likelihood for the multinomial distribution (bag of words) \[Blog\]. Jakuba. https://blog.jakuba.net/maximum-likelihood-for-multinomial-distribution/
#' 
#' Berry, K. J., & Mielke, P. W. (1995). Exact cumulative probabilities for the multinomial distribution. *Educational and Psychological Measurement, 55*(5), 769–772. doi:10.1177/0013164495055005008
#' 
#' Edgeworth, F. Y. (1905). The law of error. *Transactions of the Cambridge Philosophical Society, 20*, 36–66.
#' 
#' García-Pérez, M. A. (1999). MPROB: Computation of multinomial probabilities. *Behavior Research Methods, Instruments, & Computers, 31*(4), 701–705. doi:10.3758/BF03200749
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' freq = c(3, 6, 2, 9)
#' prob = c(0.2, 0.3, 0.1, 0.4)
#' di_mpmf(freq, prob)
#' 
#' @export
di_mpmf <- function(F, P, method='loggamma'){
  
  k = length(F)
  n = sum(F)
  
  if (method=='factorial'){
    mpmf = factorial(n)/prod(factorial(F)) * prod(P**F)    
  }
  else if (method=='gamma'){
    mpmf = gamma(1 + n)/prod(gamma(F+1)) * prod(P**F)
    }
  
  else if (method=='loggamma'){
    mpmf = exp(lgamma(n+1) + sum(F*log(P) - lgamma(F+1)))
    }
  
  else if (method=='mprob'){
    # Get the order of vec1 in descending order
    sorted_indices <- order(F, decreasing = TRUE)
    
    # Reorder vec1 and vec2 based on the sorted_indices
    F_s <- F[sorted_indices]
    P_s <- P[sorted_indices]
    
    mpmf = 1
    t = P_s[1]
    i = 2
    x = 0
    m = F_s[1]    
    
    #loop from 2 to k
    for (i in 2:k){
      #set l
      l = F_s[i]
      #for r from 1 to l
      for (r in 1:l){
        #update x
        x = x + 1
        
        #check x greater than F_1^*
        if (x > F_s[1]){
          t = 1}
        
        #update pmf
        mpmf = mpmf*t*P_s[i]*(r+m)/r
      }
      
      #update m    
      m = m + F_s[i]
    }  
    
    if (x < F_s[1]){
      for (r in x+1:F_s[1]){
        mpmf = mpmf*P_s[1]
      }
    }
  }
  
  return (mpmf)
}

#' Multinomial Cumulative Distribution Function
#' @description
#' This is a function for the cumulative multinomial probability. It returns the probability of a distribution as given in F for a sample size of sum of F, where the probability for each category is given as in P, or a distribution even more rare. It is a generalization of the binomial distribution.
#' 
#' The distribution is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Distributions/Multinomial.html)
#' 
#' 
#' @param F list with the observed counts
#' @param P list with the probabilities for each category
#' @param method optional the calculation method to use. Either "loggamma" (default), "factorial", "gamma", "mprob".
#' 
#' 
#' @returns
#' A float with the requested probability
#' 
#' 
#' @details 
#' The function first determines all possible arrangements over k categories that sum to n, using the **find_combinations()** function. It then uses the **di_pmf()** function to determine the probability for each of these, and sums those that are less or equal to the sample version.
#' 
#' This distribution is used in a Multinomial Goodness-of-Fit Test. The stikpetR library has a function \code{\link{ts_multinomial_gof}} for this, but it uses the dmultinomial function from R.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples 
#' freq = c(3, 6, 2, 9)
#' prob = c(0.2, 0.3, 0.1, 0.4)
#' di_mcdf(freq, prob)
#' 
#'  
#' @export
di_mcdf <- function(F, P, method='loggamma'){
  n = sum(F)
  k = length(F)
  mpmf = di_mpmf(F, P, method)
  
  combinations = he_find_combinations(n, k)
  mcdf = 0
  for (i in combinations){
    pmf = di_mpmf(i, P, method)
    if (pmf <= mpmf){
      mcdf = mcdf + pmf}
  }
  return (mcdf)
}

#' Find Combinations
#' @description
#' Helper function for the multinomial cumulative distribution. Will return all possible combinations to distribute n items over k categories.
#' 
#' @param n int with the sample size
#' @param k int with the number of categories
#' @returns
#' A float with the requested probability
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
he_find_combinations <- function(n, k) {
  # Helper function
  helper <- function(n, k, current_arrangement, result_env) {
    if (k == 1) {
      result_env$result <- append(result_env$result, list(c(current_arrangement, n)))
      return()
    }
    for (i in 0:n) {
      helper(n - i, k - 1, c(current_arrangement, i), result_env)
    }
  }
  
  result_env <- new.env()
  result_env$result <- list()
  helper(n, k, numeric(0), result_env)
  return(result_env$result)
}