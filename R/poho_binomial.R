#' Pairwise Binomial Test for Post-Hoc Analysis
#' 
#' @param data dataframe with scores
#' @param expCount optional pandas dataframe with categories and expected counts
#' @param twoSidedMethod method to use for determining two-sided p-value of binomial test
#' @param posthoc the correction to use, currently only "bonferroni" available
#' 
#' @returns
#' a dataframe with:
#' \item{category 1}{the label of the first category}
#' \item{category 2}{the label of the second category}
#' \item{n1}{the sample size of the first category}
#' \item{n2}{the sample size of the second category}
#' \item{obs. prop. cat. 1}{the proportion in the sample of the first category}
#' \item{exp. prop. cat. 1}{the expected proportion for the first category}
#' \item{p-value}{the unadjusted significance}
#' \item{adj. p-value}{the adjusted significance}
#' 
#' @description 
#' This function will perform a one-sample binomial test for each possible pair in the 
#' data. It makes use of the ts_binomial_os() function.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ph_binomial <- function(data, expCount=NULL, twoSidedMethod='eqdist', posthoc="bonferroni"){
  
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
  
  columns=c("category 1", "category 2", "n1", "n2", "obs. prop. cat. 1", "exp. prop. cat. 1", "p-value", "adj. p-value")
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
      sig = ts_binomial_os(data, codes=codes, p0=exP1, twoSidedMethod=twoSidedMethod)[1,1]
      
      if (posthoc=="bonferroni"){adjSig = sig*adjFactor}
      if (adjSig > 1){adjSig = 1}
      
      res[nrow(res) + 1,] = c(freq[i,1], freq[j,1], n1, n2, obP1, exP1, sig, adjSig)
    }
  }
  
  
  return (res)
}