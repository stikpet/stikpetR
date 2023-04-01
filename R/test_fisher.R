#' Fisher Exact test
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @return the two-tailed p-value/sig.
#' 
#' @examples 
#' var1 <- c("female", "female","female","female","female","female","female","female",
#'           "female","female","female", "male", "male", "male", "male", "male", "male", 
#'           "male", "male", "male", "male", "male", "male", "male", "male", "male", "male",
#'           "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#'           "male", "male", "male", "male", "male")
#' var2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other",
#'           "nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl",
#'           "other", "other", "other", "other", "other", "other", "other", "other", "other", 
#'           "other", "other", "other", "other", "other", "other")
#' ts_fisher(var1, var2)
#' 
#' @details 
#' The p-value is calculated as follows.
#' First determine the probability of the sample data cross table using the hypergeometric distribution.
#' \deqn{p_{s} = h\left(a, C_1, C_2, R_1\right)}
#' 
#' Second determine the minimum and maximum value the top-left cell could have.
#' \deqn{a_{min} = max\left(0, C_1 + R_1 - n\right)}
#' \deqn{a_{max} = min\left(C_1, R_1\right)}
#' 
#' Third determine the probability for each possible value of the top-left cell, and add those equal or less to \eqn{p_s}
#' \deqn{sig. = \sum_{i=a_{min}}^{a_{max}} \begin{cases} p_i & \text{ if } p_i \leq p_s \\ 0 & \text{ if } p_i > p_s \end{cases}}
#' With:
#' \deqn{p_i = h\left(i, C_1, C_2, R_1\right)}
#'  
#' *Symbols used:*
#' \itemize{
#' \item \eqn{p_s} the probability of the sample cross table 
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{C_i} the sum of the i-th column in the cross table 
#' \item \eqn{R_i} the sum of the i-th row in the cross table 
#' \item \eqn{h\left(\dots\right)} the probability mass function of the hypergeometric distribution 
#' \item \eqn{n} the total sample size (sum of all counts in the cross table) 
#' }
#' 
#' The test is described by Fisher (1950, p. 96).
#' 
#' For larger than 2x2 tables the Fisher-Freeman-Halton test could be used.
#' 
#' **Alternative**
#' 
#' R's *stats* library has a similar function: *fisher.test()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Fisher, R. A. (1950). *Statistical methods for research workers* (11th rev.). Oliver and Boyd.
#'  
#' @export
ts_fisher <- function(var1, var2){
  
  data = data.frame(var1, var2)
  
  #remove missing values
  data = na.omit(data)
  
  #Create a cross table first
  ct = table(data)
  
  Rs = rowSums(ct)
  Cs = colSums(ct)
  n = sum(Rs)
  
  amin = max(0, Cs[1] + Rs[1] - n)
  amax = min(Rs[1], Cs[1])
  
  pSample = dhyper(ct[1,1], Cs[1], Cs[2], Rs[1])
  
  pValue = 0
  for (i in amin:amax) {
    pFori = dhyper(i, Cs[1], Cs[2], Rs[1])
    
    if (pFori <= pSample) {
      pValue = pValue + pFori
    }
    
  }
  
  return(pValue)
}