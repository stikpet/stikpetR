#' Freeman-Tukey-Read Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test, 
#' but also an exact multinomial, G-test, Freeman-Tukey, Neyman, Mod-Log Likelihood and Cressie-Read test are possible.
#' 
#' This is actually a family (class) of tests, similar as the Cressie-Read. Weights can be chosen. the default 
#' will give the same results as the default for Cressie-Read with lambda = 0.5. Setting the weights to 4/5, 8/5, 
#' 16/15, 8/15 gives the same results as Cressie-Read with lambda = 3/2. The Pearson chi-square test is the 
#' same when setting weights to 1, 2, 1 and setting the weight simply to 4 gives the original Freeman-Tukey.
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' @param weights the weights to be used (should sum to 4)
#' @param cc Optional continuity correction. Either "none" (default), "yates", "pearson", or "williams"
#' 
#' @returns 
#' Dataframe with:
#' \item{n}{the sample size}
#' \item{k}{the number of categories}
#' \item{statistic}{the chi-square statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{two-sided p-value}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of expected counts below 5}
#' \item{testUsed}{a description of the test used}
#' 
#' 
#' @details
#' The formula used is (Read, 1987, p. 271):
#' \deqn{FT\left(b_0, b_1, \dots, b_x\right) = \sum_{i=1}^k \left(\sum_{j=0}^x b_j\times\left(\sqrt{\frac{F_i}{E_i}}\right)^j\right)\times\left(\sqrt{F_i} - \sqrt{E_i}\right)^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(FT, df\right)}
#' With, if no expected counts provided:
#' \deqn{E_i = \frac{n}{k}}
#' else:
#' \deqn{E_i = n\times\frac{E_{p_i}}{n_p}}
#' \deqn{n_p = \sum_{i=1}^k E_{p_i}}
#' The sum of the \eqn{b_i} should be four, i.e.
#' \deqn{\sum_{i=0}^x = 4}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{k} the number of categories
#' \item \eqn{F_i} the (absolute) frequency of category i
#' \item \eqn{E_i} the expected frequency of category i
#' \item \eqn{E_{p_i}} the provided expected frequency of category i
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies
#' \item \eqn{n_p} the sum of all provided expected counts
#' \item \eqn{\chi^2\left(\dots\right)}	the chi-square cumulative density function
#' }:
#' 
#' The default weights are the ones used by Read \eqn{\left(\frac{4}{3}, \frac{8}{3}\right)}, which would
#' be the same as using a Cressie-Read power divergence with \eqn{\lambda = \frac{1}{2}}
#' 
#' The Yates continuity correction (cc="yates") is calculated using (Yates, 1934, p. 222):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i > E_i \\ F_i + 0.5 & \text{ if } F_i < E_i \\ F_i & \text{ if } F_i = E_i \end{cases}}
#' \deqn{G_Y=2\times\sum_{i=1}^{k}\left(F_i^\ast\times ln\left(\frac{F_i^\ast}{E_{i}}\right)\right)}
#' Where if \eqn{F_i^\ast = 0} then \eqn{F_i^\ast\times ln\left(\frac{F_i^\ast}{E_{i}}\right) = 0}
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' The formula is also used by McDonald (2014, p. 87)
#' 
#' @references 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Read, C. B. (1993). Freeman-Tukey chi-squared goodness-of-fit statistics. *Statistics & Probability Letters, 18*(4), 271–278. https://doi.org/10.1016/0167-7152(93)90015-B
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' ts_freeman_tukey_read(ex1)
#' 
#' #Example 2: pandas series with various settings
#' ex2 = df1['mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_freeman_tukey_read(ex2, expCounts=eCounts, cc="yates")
#' ts_freeman_tukey_read(ex2, expCounts=eCounts, cc="pearson")
#' ts_freeman_tukey_read(ex2, expCounts=eCounts, cc="williams")
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_freeman_tukey_read(ex3)
#' 
#' @export
ts_freeman_tukey_read <- function(data, expCounts=NULL, weights=c(4/3, 8/3), cc = c("none", "yates", "pearson", "williams")){
  
  data = na.omit(data)
  
  if (length(cc)>1) {cc="none"}
  
  #determine the observed counts
  if (is.null(expCounts)){
    freqTable<-table(data)
    n = sum(freqTable)
    k = dim(freqTable)
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in 1:k){
      freq[nrow(freq) + 1,] = c(names(freqTable)[i], freqTable[i])
    }
    nE = n
  }
  else{
    #if expected counts are given
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCounts[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    k = nrow(expCounts)
    nE = sum(expCounts[,2])
    freq[,2] = as.numeric(freq[,2])
    n = sum(freq[,2])
  }
  
  df = k - 1
  
  if (is.null(expCounts)){
    expC = rep(n/k,k)
  }
  else{
    expC = c()
    for (i in 1:k){
      for (j in 1:k){
        if (expCounts[i,1]==freq[j,1]){
          expC = append(expC, expCounts[i,2]/nE * n)
        }
      }
    }
  }
  
  minExp = min(expC)
  propBelow5 = sum(expC < 5)/k
  
  if (cc == "yates"){
    for (i in 1:k){
      if (as.numeric(freq[i, 2] > expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) - 0.5}
      else if (as.numeric(freq[i, 2] < expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) + 0.5}
    }
  }

  r = length(weights)
  FT = 0
  for (i in 1:k){
    ff = 0
    F = as.numeric(freq[i, 2])
    for (j in 1:r) {
      ff = ff + weights[j]*sqrt(F/expC[i])**(j - 1)
    }
    FT = FT + ff*(sqrt(F) - sqrt(expC[i])) ^ 2 
  }
  
  if (cc == "pearson"){
    FT = (n - 1) / n * FT}
  else if (cc == "williams"){
    FT = FT / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}
  
  pValue = pchisq(FT, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "Freeman-Tukey-Read test of goodness-of-fit"
  if (cc == "pearson"){
    testUsed = paste(testUsed, ", with E. Pearson continuity correction")}
  else if (cc == "williams"){
    testUsed = paste(testUsed, ", with Williams continuity correction")}
  else if (cc =="yates"){
    testUsed = paste(testUsed, ", with Yates continuity correction")
  }
  
  statistic = FT
  testResults <- data.frame(n, k, statistic, df, pValue, minExp, propBelow5, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test")
  
  return (testResults)
  
}