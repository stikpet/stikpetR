#' Exact Multinomial Test of Goodness-of-Fit
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' @returns 
#' Dataframe with:
#' \item{pObs}{probability of the observed data}
#' \item{ncomb}{number of combinations used}
#' \item{pValue}{two-sided p-value}
#' \item{testUsed}{a description of the test used}
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
#' @examples
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_multinomial_gof(data)
#' ts_multinomial_gof(data, eCounts)
#' 
#' @seealso 
#' Alternative tests with a nominal variable:
#' \itemize{
#' \item \code{\link{ts_pearson_gof}} Pearson chi-square test of goodness-of-fit
#' \item \code{\link{ts_g_gof}} G / Likelihood Ratio / Wilks test of goodness-of-fit
#' \item \code{\link{ts_freeman_tukey_gof}} Freeman-Tukey test of goodness-of-fit
#' \item \code{\link{ts_neyman_gof}} Neyman test of goodness-of-fit
#' \item \code{\link{ts_mod_log_likelihood_gof}} mod-log likelihood test of goodness-of-fit
#' \item \code{\link{ts_cressie_read_gof}} Cressie-Read / Power Divergence test of goodness-of-fit
#' \item \code{\link{ts_freeman_tukey_read}} Freeman-Tukey-Read test of goodness-of-fit
#' }
#' 
#' Effect sizes that might be of interest (although they all require a chi-square test statistic):
#' \itemize{
#' \item \code{\link{es_cramer_v_gof}} Cramér's V for goodness-of-fit
#' \item \code{\link{es_cohen_w}} Cohen w
#' \item \code{\link{es_jbm_e}} Johnston-Berry-Mielke E
#' }
#' 
#' @references 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
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