#' Pairwise Goodness-of-Fit Tests for Post-Hoc Analysis
#' 
#' @description 
#' This function will perform a goodness-of-fit test for each possible pair in the data. This could be any of the goodness-of-fit tests, e.g. a Pearson chi-square.
#' 
#' The unadjusted p-values and Bonferroni adjusted p-values are both determined.
#' 
#' @param data dataframe with scores
#' @param test {"pearson", "freeman-tukey", "freeman-tukey-read", "g", "mod-log-g", "neyman", "powerdivergence", "multinomial"}, optional test to use for each pair
#' @param expCount optional dataframe with categories and expected counts
#' @param ... optional additional arguments for the specific test that are passed along.
#' 
#' @returns
#' a dataframe with:
#' \item{category 1}{the label of the first category}
#' \item{category 2}{the label of the second category}
#' \item{n1}{the sample size of the first category}
#' \item{n2}{the sample size of the second category}
#' \item{obs. prop. 1}{the observed proportion of category 1}
#' \item{exp. prop. 1}{the expected proportion of category 1}
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
#' @seealso 
#' \code{\link{ts_multinomial_gof}}
#' \code{\link{ts_powerdivergence_gof}}
#' \code{\link{ts_neyman_gof}}
#' \code{\link{ts_mod_log_likelihood_gof}}
#' \code{\link{ts_g_gof}}
#' \code{\link{ts_freeman_tukey_read}}
#' \code{\link{ts_freeman_tukey_gof}}
#' \code{\link{ts_pearson_gof}}
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_pairwise_gof <- function(data, test="pearson", expCount=NULL, ...){
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
  
  if (test=="multinomial"){  
    columns=c("category 1", "category 2", "n1", "n2", "obs. prop. 1", "exp. prop. 1", "p obs", "n combs.","p-value", "adj. p-value", "test")
  }
  else{
    columns=c("category 1", "category 2", "n1", "n2", "obs. prop. 1", "exp. prop. 1", "statistic", "df", "p-value", "adj. p-value", "minExp", "propBelow5", "test")
  }
  
  res = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) = columns
  adjFactor = k*(k - 1)/2
  for (i in 1:(k -1)){
    for (j in (i+1):k){
      cat1 = freq[i,1]
      cat2 = freq[j,1]
      dataPair = c(data[data == cat1], data[data == cat2])
      
      n1 = as.numeric(freq[i,2])
      n2 = as.numeric(freq[j,2])
      
      #determine the expected proportions for the pair
      expected_proportion_1 = expC[i]/n
      expected_proportion_2 = expC[j]/n
      exp_count_1 = (n1 + n2)*(expected_proportion_1*1/(expected_proportion_1+expected_proportion_2))
      exp_count_2 = (n1 + n2)*(expected_proportion_2*1/(expected_proportion_1+expected_proportion_2))
      exP = data.frame('category' = c(cat1, cat2), 'exp count' = c(exp_count_1, exp_count_2))
      
      if (test=="pearson"){              
        testResult = ts_pearson_gof(dataPair, expCounts=exP, ...)}
      else if (test=="freeman-tukey"){
        testResult = ts_freeman_tukey_gof(dataPair, expCounts=exP, ...)}
      else if (test=="freeman-tukey-read"){
        testResult = ts_freeman_tukey_read(dataPair, expCounts=exP, ...)}
      else if (test=="g"){
        testResult = ts_g_gof(dataPair, expCounts=exP, ...)}
      else if (test=="mod-log-g"){
        testResult = ts_mod_log_likelihood_gof(dataPair, expCounts=exP, ...)}
      else if (test=="neyman"){
        testResult = ts_neyman_gof(dataPair, expCounts=exP, ...)}
      else if (test=="powerdivergence"){
        testResult = ts_powerdivergence_gof(dataPair, expCounts=exP, ...)}
      else if (test=="multinomial"){
        testResult = ts_multinomial_gof(dataPair, expCounts=exP, ...)}
      
      if (test=="multinomial"){ 
        pObs = testResult[1,1]
        nComb = testResult[1,2]
        pValue = testResult[1,3]
        adjpValue = pValue*adjFactor
        if (adjpValue > 1){
          adjpValue = 1}
        testDesc = testResult[1,4]
        res[nrow(res) + 1,] = c(cat1, cat2, n1, n2, n1/(n1+n2), exp_count_1/(n1+n2), pObs, nComb, pValue, adjpValue, testDesc)
      }
      else {
        statistic = testResult[1,3]
        df = testResult[1,4]
        pValue = testResult[1,5]
        adjpValue = pValue*adjFactor
        if (adjpValue > 1){
          adjpValue = 1}
        minExp = testResult[1,6]
        propBelow5 = testResult[1,7]
        testDesc = testResult[1,8]
        res[nrow(res) + 1,] = c(cat1, cat2, n1, n2, n1/(n1+n2), exp_count_1/(n1+n2), statistic, df, pValue, adjpValue, minExp, propBelow5, testDesc)
      }
    }
  }
  
  return (res)
}