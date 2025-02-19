#' Freeman-Tukey Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test, but also an exact multinomial, G-test, Neyman, Mod-Log Likelihood, Cressie-Read, and Freeman-Tukey-Read test are possible.
#' 
#' The Freeman-Tukey attempts to make the distribution more like a normal distribution by using a square root transformation.
#' 
#' Lawal (1984) continued some work from Larntz (1978) and compared the modified Freeman-Tukey, G-test and the Pearson chi-square test, and concluded that for small samples the Pearson test is preferred, while for large samples either the Pearson or G-test. Making this Freeman-Tukey test perhaps somewhat redundant....
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' @param cc Optional continuity correction. Either "none" (default), "yates", "yates2", "pearson", or "williams"
#' @param modified int, optional. indicate the use of the modified version. Either 0 (default = no modification), 1 or 2
#' 
#' @returns 
#' Dataframe with:
#' \item{n}{the sample size}
#' \item{k}{the number of categories}
#' \item{statistic}{the chi-square statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{two-sided p-value}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of expected counts below 5}
#' \item{test Used}{a description of the test used}
#' 
#' 
#' @details 
#' The formula used is (Ayinde & Abidoye, 2010, p. 21):
#' \deqn{\chi_{FT}^2=\sum_{i=1}^{k}\left(\sqrt{F_{i}} - \sqrt{E_{i}}\right)^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{FT}^2,df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^k F_i}
#' If no expected counts provided:
#' \deqn{E_i = \frac{n}{k}}
#' else:
#' \deqn{E_i = n\times\frac{E_{p_i}}{n_p}}
#' \deqn{n_p = \sum_{i=1}^k E_{p_i}}
#' 
#' A modified version uses another possible smoothing (Bishop, 1969, p. 284; Larntz, 1978, p.253):
#' \deqn{\chi_{MFT}^2 = \sum_{i=1}^{k}\left(\sqrt{F_{i}} + \sqrt{F_{i} + 1} - \sqrt{4\times E_{i} + 1}\right)^2}
#' 
#' Or slightly different (Read & Cressie, 1988, p. 82):
#' \deqn{\chi_{MFT}^2 = \sum_{i=1}^{k}\left(\sqrt{F_{i}} + \sqrt{F_{i} + 1} - \sqrt{4\times \left(E_{i} + 1\right)}\right)^2}
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
#' The test is attributed to Freeman and Tukey (1950), but couldn't really find it in there. Another source often mentioned is Bishop et al. (2007)
#' 
#' The Yates continuity correction (cc="yates") is calculated using (Yates, 1934, p. 222):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i > E_i \\ F_i + 0.5 & \text{ if } F_i < E_i \\ F_i & \text{ if } F_i = E_i \end{cases}}
#' 
#' In some cases the Yates correction is slightly changed to (yates2) (Allen, 1990, p. 523):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i - 0.5 > E_i \\ F_i + 0.5 & \text{ if } F_i + 0.5 < E_i \\ F_i & \text{ else } \end{cases}}
#' 
#' Note that the Yates correction is usually only considered if there are only two categories. Some also argue this correction is too conservative (see for details Haviland (1990)).
#' 
#' The Pearson correction (cc="pearson") is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{adj}^2 = \chi_{FT}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (cc="williams") is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{adj}^2 = \frac{\chi_{FT}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' The formula is also used by McDonald (2014, p. 87)
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
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_mod_log_likelihood_gof}}, for Mod-Log Likelihood Test of Goodness-of-Fit. 
#' \code{\link{ts_multinomial_gof}}, for Multinomial Goodness-of-Fit Test. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test. 
#' 
#' @references
#' Allen, A. O. (1990). *Probability, statistics, and queueing theory with computer science applications* (2nd ed.). Academic Press.
#' 
#' Ayinde, K., & Abidoye, A. O. (2010). Simplified Freeman-Tukey test statistics for testing probabilities in contingency tables. *Science World Journal, 2*(2), 21–27. doi:10.4314/swj.v2i2.51730
#' 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations Related to the angular and the square root. *The Annals of Mathematical Statistics, 21*(4), 607–611. doi:10.1214/aoms/1177729756
#' 
#' Haviland, M. G. (1990). Yates’s correction for continuity and the analysis of 2 × 2 contingency tables. *Statistics in Medicine, 9*(4), 363–367. doi:10.1002/sim.4780090403
#'  
#'  Larntz, K. (1978). Small-sample comparisons of exact levels for chi-squared goodness-of-fit statistics. *Journal of the American Statistical Association, 73*(362), 253–263. doi:10.1080/01621459.1978.10481567
#'  
#'  Lawal, H. B. (1984). Comparisons of the X 2 , Y 2 , Freeman-Tukey and Williams’s improved G 2 test statistics in small samples of one-way multinomials. *Biometrika, 71*(2), 415–418. doi:10.2307/2336263
#'  
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. doi:10.2307/2332518
#' 
#' Read, T. R. C., & Cressie, N. A. C. (1988). Goodness-of-fit statistics for discrete multivariate data. Springer-Verlag.
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. doi:10.2307/2335081
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' ts_freeman_tukey_gof(ex1)
#' 
#' #Example 2: Dataframe with various settings
#' ex2 = df1['mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_freeman_tukey_gof(ex2, expCounts=eCounts, cc="yates")
#' ts_freeman_tukey_gof(ex2, expCounts=eCounts, cc="pearson")
#' ts_freeman_tukey_gof(ex2, expCounts=eCounts, cc="williams")
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_freeman_tukey_gof(ex3)
#' ts_freeman_tukey_gof(ex3, expCount=eCounts)
#' 
#' @export
ts_freeman_tukey_gof <- function(data, expCounts=NULL, cc = c("none", "yates", "yates2", "pearson", "williams"), modified=0){
  
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
      if (as.numeric(freq[i, 2]) > expC[i]){
        freq[i, 2] = as.numeric(freq[i, 2]) - 0.5}
      else if (as.numeric(freq[i, 2]) < expC[i]){
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
  T2 = 0
  for (i in 1:k){
    if (modified==1){
      T2 = T2 + (sqrt(as.numeric(freq[i, 2])) + sqrt(as.numeric(freq[i, 2])+1) - sqrt(4*expC[i]+1)) ^ 2
    }
    else if (modified==2){
      T2 = T2 + (sqrt(as.numeric(freq[i, 2])) + sqrt(as.numeric(freq[i, 2])+1) - sqrt(4*(expC[i]+1))) ^ 2
    }
    else if (modified==0){
    T2 = T2 + (sqrt(as.numeric(freq[i, 2])) - sqrt(expC[i])) ^ 2
    }
  }
  if (modified==0){T2 = 4*T2}
    
  if (cc == "pearson"){
    T2 = (n - 1) / n * T2}
  else if (cc == "williams"){
    T2 = T2 / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}

  pValue = pchisq(T2, df, lower.tail = FALSE)
  
  #Which test was used
  if (modified){testUsed = "modified Freeman-Tukey"}
  else{testUsed = "Freeman-Tukey"}
  if (cc == "pearson"){
    testUsed = paste0(testUsed, ", and Pearson correction")}
  else if (cc == "williams"){
    testUsed = paste0(testUsed, ", and Williams correction")}
  else if (cc =="yates" || cc=="yates2"){
    testUsed = paste0(testUsed, ", and Yates correction")
  }
  
  statistic = T2
  testResults <- data.frame(n, k, statistic, df, pValue, minExp, propBelow5, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test used")
  
  return (testResults)
  
}
