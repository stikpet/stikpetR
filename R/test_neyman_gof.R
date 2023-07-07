#' Neyman Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test, 
#' but also an exact multinomial, G-test, Freeman-Tukey, Mod-Log Likelihood, Cressie-Read, and 
#' Freeman-Tukey-Read test are possible.
#' 
#' @param data A vector or dataframe
#' @param expCounts Optional dataframe with the categories and expected counts 
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
#' The formula used is (Neyman, 1949, p. 250):
#' \deqn{\chi_{N}^{2}=\sum_{i=1}^{k}\frac{\left(O_{i}-E_{i}\right)^{2}}{O_{i}}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{N}^{2},df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^k F_i}
#' If no expected counts provided:
#' \deqn{E_i = \frac{n}{k}}
#' else:
#' \deqn{E_i = n\times\frac{E_{p_i}}{n_p}}
#' \deqn{n_p = \sum_{i=1}^k E_{p_i}}
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
#' }
#' 
#' The Yates correction (yates) is calculated using (Yates, 1934, p. 222):
#' \deqn{\chi_{NY}^2 = \sum_{i=1}^k \frac{\left(\left|F_i - E_i\right| - 0.5\right)^2}{O_i}}
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{NP}^2 = \chi_{N}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{NW}^2 = \frac{\chi_{N}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' 
#' @references 
#' Neyman, J. (1949). Contribution to the theory of the chi-square test. *Berkeley Symposium on Math. Stat, and Prob*, 239–273. https://doi.org/10.1525/9780520327016-030
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' ts_neyman_gof(ex1)
#' 
#' #Example 2: dataframe with various settings
#' ex2 = df1['mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_neyman_gof(ex2, expCounts=eCounts, cc="yates")
#' ts_neyman_gof(ex2, expCounts=eCounts, cc="pearson")
#' ts_neyman_gof(ex2, expCounts=eCounts, cc="williams")
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_neyman_gof(ex3)
#' 
#' @export
ts_neyman_gof <- function(data, expCounts=NULL, cc = c("none", "yates", "pearson", "williams")){
  
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
  
  #the degrees of freedom
  df = k - 1
  
  #the true expected counts
  if (is.null(expCounts)){
    #assume all to be equal
    expC = rep(n/k,k)
  }
  else{
    #check if categories match
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
  
  #calculate the chi-square value
  chiVal = 0
  for (i in 1:k){
    chiVal = chiVal + (as.numeric(freq[i, 2]) - expC[i]) ^ 2 / as.numeric(freq[i, 2])
  }
  
  if (cc == "pearson"){
    chiVal = (n - 1) / n * chiVal}
  else if (cc == "williams"){
    chiVal = chiVal / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}

  pValue = pchisq(chiVal, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "Neyman chi-square test of goodness-of-fit"
  if (cc == "pearson"){
    testUsed = paste0(testUsed, ", with E. Pearson continuity correction")}
  else if (cc == "williams"){
    testUsed = paste0(testUsed, ", with Williams continuity correction")}
  else if (cc == "yates"){
    testUsed = paste0(testUsed, ", with Yates continuity correction")}
  
  statistic=chiVal
  
  testResults <- data.frame(n, k, statistic, df, pValue, minExp, propBelow5, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test")
  
  return (testResults)
  
}