#' One-Sample Binomial Test
#' 
#' @param data A vector with the data
#' @param codes Optional vector with the two codes to use
#' @param p0 Optional hypothesized proportion for the first category (default is 0.5)
#' @param twoSidedMethod Optional method to be used for 2-sided significance (see details)
#' @returns 
#' Dataframe with:
#' \item{pValue}{two-sided p-value}
#' \item{testUsed}{a description of the test used}
#' 
#' @description 
#' 
#' Performs a one-sample (exact) binomial test. 
#' This test can be useful with a single binary variable as input. The null hypothesis is usually that
#' the proportions of the two categories in the population are equal (i.e. 0.5 for each). 
#' If the p-value of the test is below the pre-defined alpha level (usually 5% = 0.05) the null hypothesis
#' is rejected and the two categories differ in proportion significantly.
#' 
#' The input for the function doesn't have to be a binary variable. A nominal variable can also be used
#'  and the two categories to compare indicated. 
#' 
#' A significance in general is the probability of a result as in the sample, or more extreme, if the null
#' hypothesis is true. For a two-tailed binomial test the 'or more extreme' causes a bit of a complication.
#' There are different methods to approach this problem. See the details for more information.
#' 
#' A [YouTube](https://youtu.be/9OGCi1Q7tBQ) video on the binomial test.
#' 
#' @details 
#' A one sided p-value is calculated first:
#' \deqn{sig_{one-tail} = \text{Bin}\left(n, n_{min}, p_0^*\right)}
#' With:
#' \deqn{n_min = min\left\{n_s, n_f\right\}}
#' \deqn{p_0^* = \begin{cases}p_0 & \text{ if } n_{min}=n_s \\1 - p_0 & \text{ if } n_{min}= n_f\end{cases}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} is the number of cases
#' \item \eqn{n_s} is the number of successes
#' \item \eqn{n_f} is the number of failures
#' \item \eqn{p_0} is the probability of a success according to the null hypothesis
#' \item \eqn{p_0^*} is the probability adjusted in case failures is used
#' \item \eqn{\text{Bin}\left(\dots\right)} the binomial cumulative distribution function
#' }
#' 
#' For the two sided significance three options can be used.
#' 
#' *Option 1: Equal Distance Method (twoSidedMethod="eqdist")*
#' \deqn{sig_{two-tail} = B\left(n, n_{min}, p_0^* \right) + 1 - B\left(n, \left \lfloor 2 \times n_0 \right \rfloor - n_{min} - 1, p_0^*\right)}
#' With:
#' \deqn{n_0 = \left\lfloor n\times p_0\right\rfloor}
#' 
#' This method looks at the number of cases. In a sample of \eqn{n} people, 
#' we’d then expect \eqn{n_0 = \left\lfloor n\times p_0\right\rfloor} successes (we round the result down to the nearest integer).
#' We only had \eqn{n_{min}}, so a difference of \eqn{n_0-n_{min}}. 
#' The ‘equal distance method’ now means to look for the chance of having \eqn{k} or less, 
#' and \eqn{n_0+n_0-n_{min}=2\times n_0-n{_min}} or more. 
#' Each of these two probabilities can be found using a binomial distribution. 
#' Adding these two together than gives the two-sided significance. 
#' 
#' *Option 2: Small p-method (twoSidedMethod="smallp")*
#' \deqn{sig_{two-tail} = B\left(n, n_{min}, p_0^*\right) + \sum_{i=n_{min}+1}^n \begin{cases} 0 & \text{ if } b\left(n, i, p_0^*\right)> b\left(n, n_{min}, p_0^*\right) \\ b\left(n, i, p_0^*\right)& \text{ if } x \leq  b\left(n, i, p_0^*\right)> b\left(n, n_{min}, p_0^*\right) \end{cases}}
#' With:
#' \eqn{b\left(\dots\right)} as the binomial probability mass function.
#' 
#' This method looks at the probabilities itself. 
#' \eqn{b\left(n, n_{min},p_0^*\right)} is the probability of having exactly \eqn{n_{min}} 
#' out of a group of n, with a chance \eqn{p_0^*} each time. 
#' The method of small p-values now considers ‘or more extreme’ any number between 0 and n (the sample size) 
#' that has a probability less or equal to this. 
#' This means we need to go over each option, determine the probability and check if it is lower or equal. 
#' So, the probability of 0 successes, the probability of 1 success, etc. 
#' The sum for all of those will be the two-sided significance. 
#' We can reduce the work a little since any value below \eqn{n_{min}} , 
#' will also have a lower probability, so we only need to sum over the ones above it and add the one-sided significance to the sum of those.
#' 
#' *Option 3: Double single (twoSidedMethod="double")*
#' \deqn{sig_{two-tail} = 2\times sig_{one-tail}}
#' 
#' Fairly straight forward. Just double the one-sided significance.
#' 
#' @section Alternatives:
#' 
#' R *stats* library: binom.test()
#' 
#' @examples
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male")
#' ts_binomial_os(data, c("Female", "Male"), p0 = 0.5, twoSidedMethod="eqdist")
#' 
#' @seealso 
#' Effect size measures that could go with the binomial test are Cohen g (\code{\link{es_cohen_g}}), 
#' Cohen h (\code{\link{es_cohen_h_os}}), or the alternative ratio (\code{\link{es_alt_ratio}})
#' 
#' Other tests for a binary variable are the Score test (\code{\link{ts_score_os}}), and Wald test (\code{\link{ts_wald_os}})
#' 
#' [Companion Website](https://peterstatistics.com/CrashCourse/2-SingleVar/Binary/Binary-2a-Testing.html) on test for a single binary variable.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_binomial_os <- function(data, 
                           codes=NULL, 
                           p0 = 0.5, 
                           twoSidedMethod=c("eqdist", "double", "smallp")){
  
  if (length(twoSidedMethod)>1) {twoSidedMethod = "eqdist"}
  
  testUsed = "one-sample binomial"
  
  data = na.omit(data)
  
  #if no codes provided use first found
  if (is.null(codes)) {
    n1 = unname(table(data)[1])
    n2 = sum(table(data)) - n1
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
  }
  
  #Determine total sample size
  n<-n1 + n2
  
  minCount = n1
  ExpProp = p0
  if (n2 < n1){
    minCount = n2
    ExpProp = 1 - p0
  }
  
  #one sided test
  sig1 = pbinom(minCount,n,ExpProp)  
  
  #two sided
  if (twoSidedMethod=="double"){
    sigR = sig1
    testUsed = paste(testUsed, ", with double one-sided method", sep='')}
  else if(twoSidedMethod=="eqdist"){
    #Equal distance
    ExpCount = n * ExpProp
    Dist = ExpCount - minCount
    RightCount = ExpCount + Dist
    sigR = 1 - pbinom(RightCount - 1,n,ExpProp)
    testUsed = paste(testUsed, ", with equal-distance method", sep='')}
  else {
    #Method of small p
    binSmall = dbinom(minCount, n, ExpProp)
    sigR = 0
    for (i in (minCount + 1):n){
      binDist = dbinom(i, n, ExpProp)
      if (binDist <= binSmall){
        sigR = sigR + binDist}
    }
    testUsed = paste(testUsed, ", with small p method", sep='')
  }
  
  pValue = sig1 + sigR
  
  if (pValue>1) {
    pValue=1
  }
  
  testResults <- data.frame(pValue, testUsed)
  colnames(testResults)<-c("p-value (2-sided)", "test")
  
  return (testResults)
  
}