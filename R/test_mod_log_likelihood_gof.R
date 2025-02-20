#' Mod-Log Likelihood Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test, 
#' but also an exact multinomial, G-test, Freeman-Tukey, Neyman, Cressie-Read, and 
#' Freeman-Tukey-Read test are possible.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/RP8zn0Dhr8g) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/ModLogLikelihood.html)
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' @param cc Optional continuity correction. Either "none" (default), "yates", "pearson", or "williams"
#' 
#' @returns 
#' Dataframe with:
#' \item{n}{the sample size}
#' \item{k}{the number of categories}
#' \item{statistic}{the chi-square statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{two-sided p-value}
#' \item{minExp}{the minimum expected count}
#' \item{percBelow5}{the percentage of expected counts below 5}
#' \item{test used}{a description of the test used}
#' 
#' @details 
#' The formula used (Cressie & Read, 1984, p. 441):
#' \deqn{\chi_{MLR}^2=2\times\sum_{i=1}^{k}\left(E_{i}\times \ln\left(\frac{E_{i}}{F_{i}}\right)\right)}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{MLR}^2,df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^k F_i}
#' If no expected counts provided:
#' \deqn{E_i = \frac{n}{k}}
#' else:
#' \deqn{E_i = n\times\frac{E_{p_i}}{n_p}}
#' \deqn{n_p = \sum_{i=1}^k E_{p_i}}
#' 
#' *Symbols used*:
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
#' Cressie and Read (1984) is not the original source, but the source where I found
#' the formula.
#' 
#' The Yates continuity correction (cc="yates") is calculated using (Yates, 1934, p. 222):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i > E_i \\ F_i + 0.5 & \text{ if } F_i < E_i \\ F_i & \text{ if } F_i = E_i \end{cases}}
#' In some cases the Yates correction is slightly changed to (yates2) (Allen, 1990, p. 523):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i - 0.5 > E_i \\ F_i + 0.5 & \text{ if } F_i + 0.5 < E_i \\ F_i & \text{ else } \end{cases}}
#' 
#' \deqn{\chi_{MLRY}^2=2\times\sum_{i=1}^{k}\left(F_i^\ast\times ln\left(\frac{F_i^\ast}{E_{i}}\right)\right)}
#' Where if \eqn{E_i = 0} then \eqn{E_i\times ln\left(\frac{E_i}{F_i^\ast}\right) = 0}
#' 
#' Note that the Yates correction is usually only considered if there are only two categories. Some also argue this correction is too conservative (see for details Haviland (1990)).
#' 
#' The Pearson correction (cc="pearson") is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{MLREP}^2 = \chi_{MLR}^2\times\frac{n - 1}{n}}
#' 
#' The Williams correction (cc="williams") is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{MLRW}^2 = \frac{\chi_{MLR}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' 
#' @section Before, After and Alternatives:
#' BBefore this an impression using a frequency table or a visualisation might be helpful:
#' \code{\link{tab_frequency}}, for a frequency table
#' \code{\link{vi_bar_simple}}, for Simple Bar Chart. 
#' \code{\link{vi_cleveland_dot_plot}}, for Cleveland Dot Plot.
#' \code{\link{vi_dot_plot}}, for Dot Plot.
#' \code{\link{vi_pareto_chart}}, for Pareto Chart.
#' \code{\link{vi_pie}}, for Pie Chart.
#' 
#' After this you might an effect size measure:
#' \code{\link{es_cohen_w}}, for Cohen w.
#' \code{\link{es_cramer_v_gof}},  for Cramer's V for Goodness-of-Fit.
#' \code{\link{es_fei}}, for Fei.
#' \code{\link{es_jbm_e}}, for Johnston-Berry-Mielke E.
#' 
#' or perform a post-hoc test:
#' \code{\link{ph_pairwise_bin}}, for Pairwise Binary Tests.
#' \code{\link{ph_pairwise_gof}}, for Pairwise Goodness-of-Fit Tests.
#' \code{\link{ph_residual_gof_bin}}, for Residuals Tests using Binary tests.
#' \code{\link{ph_residual_gof_gof}}, for Residuals Using Goodness-of-Fit Tests.
#' 
#' Alternative tests:
#' \code{\link{ts_pearson_gof}}, for Pearson Chi-Square Goodness-of-Fit Test. 
#' \code{\link{ts_freeman_tukey_gof}}, for Freeman-Tukey Test of Goodness-of-Fit. 
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_multinomial_gof}}, for Multinomial Goodness-of-Fit Test. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test. 
#' 
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
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' ts_mod_log_likelihood_gof(ex1)
#' 
#' #Example 2: dataframe with various settings
#' ex2 = df1['mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_mod_log_likelihood_gof(ex2, expCounts=eCounts, cc="yates")
#' ts_mod_log_likelihood_gof(ex2, expCounts=eCounts, cc="pearson")
#' ts_mod_log_likelihood_gof(ex2, expCounts=eCounts, cc="williams")
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", 
#' "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_mod_log_likelihood_gof(ex3)
#' 
#' @export
ts_mod_log_likelihood_gof <- function(data, expCounts=NULL, cc = c("none", "yates", "yates2", "pearson", "williams")){
  
  data = na.omit(data)
  
  if (length(cc)>1) {cc="none"}
  
  #determine the observed counts
  if (is.null(expCounts)){
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
    for (i in expCounts[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    k = nrow(expCounts)
    nE = sum(expCounts[,2])
    freq[,2] = as.numeric(freq[,2])
    n = sum(freq[,2])
  }
  
  #the degrees of freedom
  df = k - 1
  
  #the true expected counts
  if (is.null(expCounts)){
    #assume all to be equal
    expC = rep(n/k,k)
  }
  else{
    #check if categories match
    expC = c()
    for (i in 1:k){
      for (j in 1:k){
        if (expCounts[i,1]==freq[j,1]){
          expC = append(expC, expCounts[i,2]/nE * n)
        }
      }
    }
  }
  
  minExp = min(expC)
  propBelow5 = sum(expC < 5)/k
  
  if (cc == "yates"){
    for (i in 1:k){
      if (as.numeric(freq[i, 2] > expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) - 0.5}
      else if (as.numeric(freq[i, 2] < expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) + 0.5}
    }
  }
  else if (cc == "yates2"){
    for (i in 1:k){
      if (as.numeric(freq[i, 2]) - 0.5 > expC[i]){
        freq[i, 2] = as.numeric(freq[i, 2]) - 0.5}
      else if (as.numeric(freq[i, 2]) + 0.5 < expC[i]){
        freq[i, 2] = as.numeric(freq[i, 2]) + 0.5}
    }
  }
  
  #calculate the chi-square value
  chiVal = 0
  for (i in 1:k){
    obs = as.numeric(freq[i, 2])
    chiVal = chiVal + expC[i]*log(expC[i]/obs)
  }
  chiVal = 2*chiVal
  
  if (cc == "pearson"){
    chiVal = (n - 1) / n * chiVal}
  else if (cc == "williams"){
    chiVal = chiVal / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}
  
  
  pValue = pchisq(chiVal, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "mod-log likelihood goodness-of-fit test"
  if (cc == "pearson"){
    testUsed = paste0(testUsed, ", and Pearson correction")}
  else if (cc == "williams"){
    testUsed = paste0(testUsed, ", and Williams correction")}
  else if (cc == "yates" || cc=="yates2"){
    testUsed = paste0(testUsed, ", and Yates correction")}
  
  statistic = chiVal
  
  testResults <- data.frame(n, k, statistic, df, pValue, minExp, propBelow5*100, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "percBelow5", "test used")
  
  return (testResults)
  
}