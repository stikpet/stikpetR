#' Exact Multinomial Test of Goodness-of-Fit
#' 
#' @param data A vector with the data
#' @param expCount Optional dataframe with the categories and expected counts 
#' @return dataframe with probability of the observed data, number of combinations used, p-value and test used
#' 
#' @examples
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_multinomial_gof(data)
#' ts_multinomial_gof(data, eCounts)
#' 
#' @details 
#' The exact multinomial test of goodness of fit is done in four steps
#' 
#' Step 1: Determine the probability of the observed counts using 
#' the probability mass function of the multinomial distribution
#' 
#' Step 2: Determine all possible permutations with repetition 
#' that create a sum equal to the sample size over the k-categories.
#' 
#' Step 3: Determine the probability of each of these permutations using the probability mass function of the multinomial distribution.
#' 
#' Step 4: Sum all probabilities found in step 3 that are equal or less than the one found in step 1.
#' 
#' **Alternatives**
#' 
#' The *EMT* library has a similar function: *multinomial.test()*
#' 
#' The *XNomial* library has a similar function: *xmulti()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
ts_multinomial_gof <- function(data, expCount=NULL) {
  
  #determine the observed counts
  if (is.null(expCount)){
    observed <- c(table(data))}
  else {
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCount[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    
    observed <- setNames(as.integer(t(freq[2])), t(freq[1]))
    
  }
    
  n <- sum(observed)
  k <- length(observed)
  expProb =rep(1/k, k)
  
  pObs = dmultinom(sort(observed, decreasing=TRUE), size=n, expProb)
  counts <- seq(0, n, by = 1)
  
  kCounts <- matrix(0 ,nrow=n+1, ncol=k)
  for (i in 1:k){
    kCounts[,i] <- counts
  }
  
  all_perm <- merge(kCounts[,1], as.data.frame(kCounts[,2]),all=TRUE)
  all_perm <- all_perm[rowSums(all_perm) <= n,]
  for (i in 3:k){
    set_to_merge <- data.frame(dummyname=kCounts[,i])
    colnames(set_to_merge) <- paste0("k_",i)
    all_perm <- merge(all_perm, set_to_merge,all=TRUE)
    all_perm <- all_perm[rowSums(all_perm) <= n,]
  }
  
  all_perm <- all_perm[rowSums(all_perm) == n,]
  ncomb <- dim(all_perm)[1]
  pObsAll <- apply(all_perm, 1, function(x) dmultinom(sort(x, decreasing=TRUE), size=n, expProb))
  dfPs <- data.frame(pObsAll)
  pValue <- sum(dfPs[which(round(dfPs$pObsAll, digits=8) <= round(pObs, 8)),1])
  
  testUsed <- "one-sample multinomial exact goodness-of-fit test"
  testResults <- data.frame(pObs, ncomb, pValue, testUsed)
  
  return(testResults)                 
}