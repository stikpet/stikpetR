#' Pairwise Binary Test for Post-Hoc Analysis
#' 
#' @description 
#' This function will perform a one-sample binary test for each possible pair in the data. This could either be a binomial, Wald or score test.
#' 
#' The unadjusted p-values and Bonferroni adjusted p-values are both determined.
#' 
#' @param data dataframe with scores
#' @param test {"binomial", "score", "wald"}, optional test to use for each pair
#' @param expCount optional dataframe with categories and expected counts
#' @param mtc optional string. Any of the methods available in p_adjust() to correct for multiple tests
#' @param ... optional additional arguments for the specific test that are passed along.
#' 
#' @returns
#' a dataframe with:
#' \item{category 1}{the label of the first category}
#' \item{category 2}{the label of the second category}
#' \item{n1}{the sample size of the first category}
#' \item{n2}{the sample size of the second category}
#' \item{n pair}{the sample size of of the pair}
#' \item{obs. prop. 1}{the proportion in the sample of the first category}
#' \item{exp. prop. 1}{the expected proportion for the first category}
#' \item{statistic}{the test statistic}
#' \item{p-value}{the unadjusted significance}
#' \item{adj. p-value}{the adjusted significance}
#' \item{test}{description of the test used}
#' 
#' @section Before, After and Alternatives:
#' Before this an omnibus test might be helpful:
#' \code{\link{ts_freeman_tukey_gof}}, for Freeman-Tukey Test of Goodness-of-Fit. 
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_mod_log_likelihood_gof}}, for Mod-Log Likelihood Test of Goodness-of-Fit. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_pearson_gof}}, for Pearson Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test.
#' 
#' After this you might want to add an effect size measure:
#' \code{\link{es_post_hoc_gof}} for various effect sizes.
#' 
#' Alternative post-hoc tests:
#' \code{\link{ph_pairwise_gof}} for Pairwise Goodness-of-Fit Tests.
#' \code{\link{ph_residual_gof_bin}} for Residuals Tests using Binary Tests
#' \code{\link{ph_residual_gof_gof}} for Residuals Using Goodness-of-Fit Tests
#' 
#' The binary test that is performed on each pair:
#' \code{\link{ts_binomial_os}} for One-Sample Binomial Test.
#' \code{\link{ts_score_os}} for One-Sample Score Test.
#' \code{\link{ts_wald_os}} for One-Sample Wald Test.
#' 
#' More info on the adjustment for multiple testing:
#' \code{\link{p_adjust}}, various adjustment methods.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' # Examples: get data
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' gssDf <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = gssDf['mar1']
#' 
#' #Example 1 using default settings (one-sample binomial tests with equal-distance method)
#' ph_pairwise_bin(ex1)
#' 
#' #Example 2 using a score test with Yates correction:
#' ph_pairwise_bin(ex1, test="score", mtc='holm', cc='yates')
#' 
#' @export
ph_pairwise_bin <- function(data, test="binomial", expCount=NULL, mtc='bonferroni', ...){
  data = na.omit(data)
  
  #the sample size n
  n = length(data)
  
  #determine the observed counts
  if (is.null(expCount)){
    #generate frequency table
    freqTable<-table(data)
    
    #number of categories to use (k)
    k = dim(freqTable)
    
    #rewrite frequency table in long format
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in 1:k){
      freq[nrow(freq) + 1,] = c(names(freqTable)[i], freqTable[i])
    }
    
    #number of expected counts is simply sample size
    nE = n
  }
  
  else{
    #if expected counts are given
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCount[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    
    #number of categories to use (k)
    k = nrow(expCount)
    
    #sum of the given expected counts
    nE = sum(expCount[,2])
  }
  
  #the true expected counts
  if (is.null(expCount)){
    #assume all to be equal
    expC = rep(n/k,k)
  }
  else{
    #check if categories match
    expC = c()
    for (i in 1:k){
      for (j in 1:k){
        if (expCount[i,1]==freq[j,1]){
          expC = append(expC, expCount[i,2]/nE * n)
        }
      }
    }
  }
  
  columns=c("category 1", "category 2", "n1", "n2", "n pair", "obs. prop. 1", "exp. prop. 1", "statistic", "p-value", "adj. p-value", "test")
  res = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) = columns
  adjFactor = k*(k - 1)/2
  for (i in 1:(k -1)){
    for (j in (i+1):k){
      n1 = as.numeric(freq[i,2])
      n2 = as.numeric(freq[j,2])
      obP1 = n1/(n1 + n2)
      exP1 = expC[i]/(expC[i] + expC[j])
      codes = c(freq[i,1], freq[j,1])
      
      if (test=="binomial"){              
        testResults = ts_binomial_os(data, codes=codes, p0=exP1, ...)
        statistic = "n.a."
        sig = testResults[1,1]
        testDescription = testResults[1,2]}
      else{ 
        if (test=="wald"){
          testResults = ts_wald_os(data, codes=codes, p0=exP1, ...)}
        else if (test=="score"){
          testResults = ts_score_os(data, codes=codes, p0=exP1, ...)}
        
        statistic = testResults[1,2]
        sig = testResults[1,3]
        testDescription = testResults[1,4]}
      
      adjSig = sig*adjFactor
      if (adjSig > 1){adjSig = 1}
      
      res[nrow(res) + 1,] = c(freq[i,1], freq[j,1], n1, n2, n1+n2, obP1, exP1, statistic, sig, adjSig, testDescription)
    }
  }
  
  p_adj = p_adjust(as.numeric(res[, 9]), method=mtc)
  res[, 10] = p_adj
  
  return (res)
}