#' Exact Multinomial Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is a Pearson chi-square test, 
#' but also a G-test, Freeman-Tukey, Neyman, Mod-Log Likelihood and Cressie-Read test are possible.
#' 
#' McDonald (2014, p. 82) suggests to always use this exact test as long as the sample size is less than 1000 
#' (which was just picked as a nice round number, when n is very large the exact test becomes 
#' computational heavy even for computers).
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' 
#' @returns 
#' Dataframe with:
#' \item{pObs}{probability of the observed data}
#' \item{ncomb}{number of combinations used}
#' \item{pValue}{two-sided p-value}
#' \item{testUsed}{a description of the test used}
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
#' @section Alternatives:
#' 
#' The *EMT* library has a similar function: *multinomial.test()*
#' 
#' The *XNomial* library has a similar function: *xmulti()*
#' 
#' @references 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1[1:20, 'mar1']
#' ts_multinomial_gof(ex1)
#' 
#' #Example 2: dataframe with various settings
#' ex2 = df1[1:20, 'mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_multinomial_gof(ex2, expCounts=eCounts)
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_multinomial_gof(ex3)
#' 
#' @export
ts_multinomial_gof <- function(data, expCounts=NULL) {
  
  data = na.omit(data)
  
  #determine the observed counts
  if (is.null(expCounts)){observed <- c(table(data))}
  else {
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCounts[,1]){
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
  for (i in 1:k){kCounts[,i] <- counts}
  
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