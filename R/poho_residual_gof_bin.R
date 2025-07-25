#' Post-Hoc Residuals Using Binary Tests for GoF
#' 
#' @description 
#' This function will perform a residuals post-hoc test for each of the categories in a nominal field. This could either be a z-test using the standardized residuals, the adjusted residuals, or any of the one-sample binary tests.
#' 
#' The unadjusted p-values and Bonferroni adjusted p-values are both determined.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/mioMqScuhZ8) and test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/PostHocAfterGoF.html)
#' 
#' @param data dataframe with scores
#' @param test c("adj-residual", "std-residual", "binomial", "wald", "score") optional test to use
#' @param expCount optional dataframe with categories and expected counts
#' @param mtc optional string. Any of the methods available in p_adjust() to correct for multiple tests
#' @param ... optional additional parameters to be passed to the test
#' 
#' @returns
#' a dataframe with:
#' \item{category}{the label of the first category}
#' \item{obs. count}{the observed count}
#' \item{exp. count}{the expected count}
#' \item{statistic}{the test statistic}
#' \item{p-value}{the unadjusted significance}
#' \item{adj. p-value}{the adjusted significance}
#' \item{test}{description of the test used}
#' 
#' @details
#' The formula used is for the adjusted residual test:
#' \deqn{z = \frac{F_i - E_i}{\sqrt{E_i\times\left(1 - \frac{E_i}{n}\right)}}}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' The formula used for the standardized residual test:
#' \deqn{z = \frac{F_i - E_i}{\sqrt{E_i}}}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \itemize{
#' \item \eqn{F_i}, the observed count for category $i$
#' \item \eqn{E_i}, the expected count for category $i$
#' \item \eqn{\Phi\left(\dots\right)}, the cumulative distribution function of the standard normal distribution
#' }
#' 
#' If no expected counts are provide it is assumed they are all equal for each category, i.e. \eqn{E_i = \frac{n}{k}}
#' 
#' The Bonferroni adjustment is calculated using:
#' \deqn{p_{adj} = \min \left(p \times k, 1\right)}
#' 
#' The other tests use the formula from the one-sample test variant, using the expected count/n as the expected proportion.
#' 
#' The adjusted residuals will gave the same result as using a one-sample score test. Some sources will also call these adjusted residuals as standardized residuals (Agresti, 2007, p. 38), and the standardized residuals used in this function as Pearson residuals (R, n.d.). Haberman (1973, p. 205) and Sharpe (2015, p. 3) are sources for the terminology used in this function. 
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
#' \code{\link{ph_pairwise_bin}} for Pairwise Binary Tests.
#' \code{\link{ph_pairwise_gof}} for Pairwise Goodness-of-Fit Tests.
#' \code{\link{ph_residual_gof_gof}} for Residuals Using Goodness-of-Fit Tests
#' 
#' The binary test that is performed on each category:
#' \code{\link{ts_binomial_os}} for One-Sample Binomial Test.
#' \code{\link{ts_score_os}} for One-Sample Score Test.
#' \code{\link{ts_wald_os}} for One-Sample Wald Test.
#' 
#' More info on the adjustment for multiple testing:
#' \code{\link{p_adjust}}, various adjustment methods.
#'   
#' @references 
#' Agresti, A. (2007). *An introduction to categorical data analysis* (2nd ed.). Wiley-Interscience.
#' Haberman, S. J. (1973). The analysis of residuals in cross-classified tables. *Biometrics, 29*(1), 205-220. doi:10.2307/2529686
#' R. (n.d.). Chisq.test \[Computer software\]. https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html
#' Sharpe, D. (2015). Your chi-square test is statistically significant: Now what? Practical Assessment, *Research & Evaluation, 20*(8), 1-10.
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
#' ph_residual_gof_bin(ex1)
#' 
#' #Example 2 using a binomial test and Holm correction:
#' ph_residual_gof_bin(ex1, test="binomial", mtc='holm')
#' 
#' @export
ph_residual_gof_bin <- function(data, test="std-residual", expCount=NULL, mtc='bonferroni', ...){
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
  
  columns=c("category", "obs. count", "exp. count", "statistic","p-value", "adj. p-value", "test")
  
  res = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) = columns
  adjFactor = k*(k - 1)/2
  for (i in 1:k){
    cat = freq[i,1]
    n1 = as.numeric(freq[i,2])
    e1 = expC[i]
    if (test=="std-residual"){
      statistic = (as.numeric(freq[i,2]) - expC[i])/(expC[i]**0.5)
      sig = 2*(1 - pnorm(abs(statistic)))
      testDescription = "standardized residuals z-test"
    }
    else if (test=="adj-residual"){
      statistic = (as.numeric(freq[i,2]) - expC[i])/((expC[i]*(1 - expC[i]/n))**0.5)
      sig = 2*(1 - pnorm(abs(statistic)))
      testDescription = "adjusted residuals z-test"
    }
    else {
      tempA = rep(cat, n1)
      tempB = rep('all other', n - n1)
      tempData = c(tempA, tempB)
      if (test=="binomial"){
        
        testResults = ts_binomial_os(tempData, p0 = e1/n, ...)
        statistic = "n.a."
        sig = testResults[1,1]
        testDescription = testResults[1,2]
      }
      else{
        if (test=="wald"){
          testResults = ts_wald_os(tempData, codes=c(cat, 'all other'), p0=e1/n, ...)}
        else if (test=="score"){
          testResults = ts_score_os(tempData, codes=c(cat, 'all other'), p0=e1/n, ...)}
        statistic = testResults[1,2]
        sig = testResults[1,3]
        testDescription = testResults[1,4]
      }
    }
    adjSig = min(sig*k, 1)
    res[i,] = c(cat, n1, e1, statistic, sig, adjSig, testDescription)
  }
  
  p_adj = p_adjust(as.numeric(res[, 5]), method=mtc)
  res[, 6] = p_adj
  
  return (res)
}



