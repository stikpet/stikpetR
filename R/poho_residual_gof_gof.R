#' Post-Hoc Residuals Using GoF for GoF
#' 
#' @description 
#' This function will perform a goodness-of-fit test using each category and collapsing the other categories.
#' 
#' The unadjusted p-values and Bonferroni adjusted p-values are both determined.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/RQw09eg9Ojo) and test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/PostHocAfterGoF.html)
#' 
#' @param data dataframe with scores
#' @param test c("pearson", "freeman-tukey", "freeman-tukey-read", "g", "mod-log-g", "neyman", "powerdivergence", "multinomial") optional test to use
#' @param expCount optional dataframe with categories and expected counts
#' @param mtc optional string. Any of the methods available in p_adjust() to correct for multiple tests
#' @param ... optional additional parameters to be passed to the test
#' 
#' @returns
#' a dataframe with:
#' \item{category}{the label of the first category}
#' \item{obs. count}{the observed count of the category}
#' \item{exp. count}{the expected count of the category}
#' \item{statistic}{the chi-square test statistic}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the unadjusted significance}
#' \item{adj. p-value}{the adjusted significance}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of cells with an expected count below 5}
#' \item{test}{description of the test used}
#' 
#' In case of a multinomial test, the same columns except there are no *minExp* and *propBelow5* columns and:
#' \item{p obs}{instead of *statistic*, showing the probability of the observed sample table}
#' \item{n combs.}{instead of *df*, showing the number of possible tables}
#' 
#' @section Before, After and Alternatives:
#' Before this an omnibus test might be helpful, these are also the tests used on each category:
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
#' \code{\link{ph_pairwise_bin}} for Pairwise Binary Tests.
#' \code{\link{ph_pairwise_gof}} for Pairwise Goodness-of-Fit Tests.
#' \code{\link{ph_residual_gof_bin}} for Residuals Tests using Binary Tests
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
#' #Example 1 using default settings
#' ph_residual_gof_gof(ex1)
#' 
#' #Example 2 using a G test and Holm correction:
#' ph_residual_gof_gof(ex1, test="g", mtc='holm')
#' 
#' @export
ph_residual_gof_gof <- function(data, test="pearson", expCount=NULL, mtc='bonferroni', ...){
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
  }
  
  else{
    #if expected counts are given
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCount[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    
    #number of categories to use (k)
    k = nrow(expCount)    
  }
  
  n = sum(as.numeric(freq[,2]))  
  
  
  #the true expected counts
  if (is.null(expCount)){
    #assume all to be equal
    expC = rep(n/k,k)
    #number of expected counts is simply sample size
    nE = n  
  }
  else{
    #sum of the given expected counts
    nE = sum(expCount[,2])
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
  
  if (test=="multinomial"){  
    columns=c("category", "obs. count", "exp. count", "p obs", "n combs.","p-value", "adj. p-value", "test")
  }
  else{
    columns=c("category", "obs. count", "exp. count", "statistic", "df", "p-value", "adj. p-value", "minExp", "propBelow5", "test")
  }
  
  res = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) = columns
  adjFactor = k
  for (i in 1:k){
    cat = freq[i,1]
    n1 = as.numeric(freq[i,2])
    e1 = expC[i]
    tempA = rep(cat, n1)
    tempB = rep('all other', n - n1)
    tempData = c(tempA, tempB)
    exP = data.frame('category' = c(cat, 'all other'), 'exp count' = c(e1, n - e1)) 
    if (test=="pearson"){              
      testResult = ts_pearson_gof(tempData, expCounts=exP, ...)}
    else if (test=="freeman-tukey"){
      testResult = ts_freeman_tukey_gof(tempData, expCounts=exP, ...)}
    else if (test=="freeman-tukey-read"){
      testResult = ts_freeman_tukey_read(tempData, expCounts=exP, ...)}
    else if (test=="g"){
      testResult = ts_g_gof(tempData, expCounts=exP, ...)}
    else if (test=="mod-log-g"){
      testResult = ts_mod_log_likelihood_gof(tempData, expCounts=exP, ...)}
    else if (test=="neyman"){
      testResult = ts_neyman_gof(tempData, expCounts=exP, ...)}
    else if (test=="powerdivergence"){
      testResult = ts_powerdivergence_gof(tempData, expCounts=exP, ...)}
    else if (test=="multinomial"){
      testResult = ts_multinomial_gof(tempData, expCounts=exP, ...)}
    if (test=="multinomial"){ 
      pObs = testResult[1,1]
      nComb = testResult[1,2]
      pValue = testResult[1,3]
      adjpValue = pValue
      testDesc = testResult[1,4]
      res[nrow(res) + 1,] = c(cat, n1/n, e1/n, pObs, nComb, pValue, adjpValue, testDesc)
    }
    else {
      statistic = testResult[1,3]
      df = testResult[1,4]
      pValue = testResult[1,5]
      adjpValue = pValue
      minExp = testResult[1,6]
      propBelow5 = testResult[1,7]
      testDesc = testResult[1,8]
      res[nrow(res) + 1,] = c(cat, n1, e1, statistic, df, pValue, adjpValue, minExp, propBelow5, testDesc)
    }
  }
  
  p_adj = p_adjust(as.numeric(res[, 6]), method=mtc)
  res[, 7] = p_adj
  
  return (res)
}



