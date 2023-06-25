#' Freeman-Tukey Test of Goodness-of-Fit
#' 
#' @param data A vector with the data
#' @param expCount Optional dataframe with the categories and expected counts 
#' @param cc Optional continuity correction (default is "none")
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
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test, 
#' but also an exact multinomial, G-test, Neyman, Mod-Log Likelihood, Cressie-Read, and 
#' Freeman-Tukey-Read test are possible.
#' 
#' @details 
#' The formula used is (Ayinde & Abidoye, 2010, p. 21):
#' \deqn{T^{2}=\sum_{i=1}^{k}\left(\sqrt{F_{i}} - \sqrt{E_{i}}\right)^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(T^{2},df\right)}
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
#' The test is attributed to Freeman and Tukey (1950), but couldn't really find it in there.
#' Another source often mentioned is Bishop et al. (2007)
#' 
#' The Pearson correction (cc="pearson") is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (cc="williams") is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' The formula is also used by McDonald (2014, p. 87)
#' 	
#' @examples  
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_freeman_tukey_gof(data)
#' ts_freeman_tukey_gof(data, cc="pearson")
#' ts_freeman_tukey_gof(data, cc="williams")
#' ts_freeman_tukey_gof(data, eCounts)
#' 
#' @seealso 
#' Alternative tests with a nominal variable:
#' \itemize{
#' \item \code{\link{ts_pearson_gof}} Pearson chi-square test of goodness-of-fit
#' \item \code{\link{ts_multinomial_gof}} exact multinomial test of goodness-of-fit
#' \item \code{\link{ts_g_gof}} G / Likelihood Ratio / Wilks test of goodness-of-fit
#' \item \code{\link{ts_neyman_gof}} Neyman test of goodness-of-fit
#' \item \code{\link{ts_mod_log_likelihood_gof}} mod-log likelihood test of goodness-of-fit
#' \item \code{\link{ts_cressie_read_gof}} Cressie-Read / Power Divergence test of goodness-of-fit
#' \item \code{\link{ts_freeman_tukey_read}} Freeman-Tukey-Read test of goodness-of-fit
#' }
#' 
#' Effect sizes that might be of interest:
#' \itemize{
#' \item \code{\link{es_cramer_v_gof}} Cramér's V for goodness-of-fit
#' \item \code{\link{es_cohen_w}} Cohen w
#' \item \code{\link{es_jbm_e}} Johnston-Berry-Mielke E
#' }
#' 
#' @references
#' Ayinde, K., & Abidoye, A. O. (2010). Simplified Freeman-Tukey test statistics for testing probabilities in contingency tables. *Science World Journal, 2*(2), 21–27. https://doi.org/10.4314/swj.v2i2.51730
#' 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations Related to the angular and the square root. *The Annals of Mathematical Statistics, 21*(4), 607–611. https://doi.org/10.1214/aoms/1177729756
#'  
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_freeman_tukey_gof <- function(data, expCount=NULL, cc = c("none", "pearson", "williams")){
  
  if (length(cc)>1) {cc="none"}
  
  #determine the observed counts
  if (is.null(expCount)){
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
    for (i in expCount[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    k = nrow(expCount)
    nE = sum(expCount[,2])
    n = sum(freq[,1])
  }
  
  #the degrees of freedom
  df = k - 1
  
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
  
  minExp = min(expC)
  propBelow5 = sum(expC < 5)/k
  
  #calculate the chi-square value
  T2 = 0
  for (i in 1:k){
    T2 = T2 + (sqrt(as.numeric(freq[i, 2])) - sqrt(expC[i])) ^ 2
  }
  T2 = 4*T2
    
  if (cc == "pearson"){
    T2 = (n - 1) / n * T2}
  else if (cc == "williams"){
    T2 = T2 / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}

  pValue = pchisq(T2, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "Freeman-Tukey test of goodness-of-fit"
  if (cc == "pearson"){
    testUsed = paste(testUsed, ", with E. Pearson continuity correction")}
  else if (cc == "williams"){
    testUsed = paste(testUsed, ", with Williams continuity correction")}
  
  statistic = T2
  testResults <- data.frame(n, k, statistic, df, pValue, minExp, propBelow5, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test")
  
  return (testResults)
  
}
