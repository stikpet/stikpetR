#' Cressie-Read Test of Goodness-of-Fit
#' 
#' @param data A vector with the data
#' @param expCount Optional dataframe with the categories and expected counts 
#' @param cc Optional continuity correction. Either "none" (default), "yates", "yates2", "pearson", or "williams"
#' @param lambda optional power to use in equation (see details)
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
#' but also an exact multinomial, G-test, Freeman-Tukey, Neyman, Mod-Log Likelihood, and 
#' Freeman-Tukey-Read test are possible.
#' 
#' @details 
#' The formula used is (Cressie & Read, 1984, p. 442):
#' \deqn{\chi_{C}^{2} = \begin{cases} 2\times\sum_{i=1}^{k}\left(F_{i}\times ln\left(\frac{F_{i}}{E_{i}}\right)\right) & \text{ if } \lambda=0 \\ 2\times\sum_{i=1}^{k}\sum_{j=1}^c\left(E_{i}\times ln\left(\frac{E_{i}}{F_{i}}\right)\right) & \text{ if } \lambda=-1 \\ \frac{2}{\lambda\times\left(\lambda + 1\right)} \times \sum_{i=1}^{k} F_i\times\left(\left(\frac{F_i}{E_i}\right)^{\lambda} - 1\right) & \text{ else } \end{cases}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{C}^{2},df\right)}
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
#' Cressie and Read (1984, p. 463) suggest to use \eqn{\lambda = \frac{2}{3}},  which
#' is therefor the default in this function.
#' 
#' Note that
#' 
#' The Yates correction (yates) is calculated using (Yates, 1934, p. 222):
#' \deqn{\chi_{CY}^2 = \sum_{i=1}^k \frac{\left(\left|F_i - E_i\right| - 0.5\right)^2}{E_i}}
#' 
#' In some cases the Yates correction is slightly changed to (yates2) (Allen, 1990, p. 523):
#' \deqn{\chi_{PY}^2 = \sum_{i=1}^k \frac{\max\left(0, \left(\left|F_i - E_i\right| - 0.5\right)\right)^2}{E_i}}
#' 
#' Note that the Yates correction is usually only considered if there are only two categories. Some also argue this correction is too conservative (see for details Haviland (1990)).
#' 
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{CP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{CW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' 
#'  
#' @seealso 
#' Alternative tests with a nominal variable:
#' \itemize{
#' \item \code{\link{ts_pearson_gof}} Pearson chi-square test of goodness-of-fit
#' \item \code{\link{ts_multinomial_gof}} exact multinomial test of goodness-of-fit
#' \item \code{\link{ts_g_gof}} G / Likelihood Ratio / Wilks test of goodness-of-fit
#' \item \code{\link{ts_freeman_tukey_gof}} Freeman-Tukey test of goodness-of-fit
#' \item \code{\link{ts_neyman_gof}} Neyman test of goodness-of-fit
#' \item \code{\link{ts_mod_log_likelihood_gof}} Mod-Log Likelihood test of goodness-of-fit
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
#' @section Alternatives:
#' 
#' The *MSCquartets* library has a *powerDivStat()* function, which can return the test statistics, 
#' based on given observed and expected counts.
#' obs = as.vector(unname(table(nomData)))
#' 
#' k = length(obs)
#' 
#' n = sum(obs)
#' 
#' exp = rep(1/k, k)
#' 
#' n*powerDivStat(obs/n, exp, lambda=2/3)
#' 
#' @references 
#' Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. *Journal of the Royal Statistical Society: Series B (Methodological), 46*(3), 440–464. doi:10.1111/j.2517-6161.1984.tb01318.x
#' 
#' Haviland, M. G. (1990). Yates’s correction for continuity and the analysis of 2 × 2 contingency tables. *Statistics in Medicine, 9*(4), 363–367. doi:10.1002/sim.4780090403
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. doi:10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. doi:10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. doi:10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_cressie_read_gof <- function(data, expCount=NULL, cc = c("none", "yates", "pearson", "williams"), lambda=2/3){
  
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
  chiVal = 0
  if (cc == "none" || cc == "pearson" || cc == "williams"){
    for (i in 1:k){
      obs  = as.numeric(freq[i, 2])
      ex = expC[i]
      
      if (lambda==-1) {
        chiVal = chiVal + ex*log(ex/obs)
      }
      
      else if (lambda!=0) {
        chiVal = chiVal + obs*((obs/ex)**(lambda) - 1)
      }
      else{
        chiVal = chiVal + obs*log(obs/ex)
      }
      
    }
    
    if (lambda==-1 || lambda==0) {
      chiVal = chiVal = 2*chiVal
    }
    
    else {
      chiVal = 2/(lambda*(lambda + 1)) * chiVal
    }

    if (cc == "pearson"){
      chiVal = (n - 1) / n * chiVal}
    else if (cc == "williams"){
      chiVal = chiVal / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}
  }
  else if (cc == "yates"){
    for (i in 1:k){
      chiVal = chiVal + (abs(as.numeric(freq[i, 2]) - expC[i]) - 0.5) ^ 2 / expC[i]}
  }
  else if (cc == "yates2"){
    for (i in 1:k){
      chiVal = chiVal + (max(0, abs(as.numeric(freq[i, 2]) - expC[i]) - 0.5)) ^ 2 / expC[i]}
  }
  pValue = pchisq(chiVal, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "Cressie-Read test of goodness-of-fit"
  if (lambda!=2/3) {
    testUsed = paste0(testUsed, " (lambda = ", lambda, ")")}
  
  if (cc == "pearson"){
    testUsed = paste0(testUsed, ", with E. Pearson continuity correction")}
  else if (cc == "williams"){
    testUsed = paste0(testUsed, ", with Williams continuity correction")}
  else if (cc == "yates"){
    testUsed = paste0(testUsed, ", with Yates continuity correction")}
  
  testResults <- data.frame(n, k, chiVal, df, pValue, minExp, propBelow5, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test")
  
  return (testResults)
  
}
