#' G (Likelihood Ratio) Test of Goodness-of-Fit
#' 
#' @param data A vector with the data
#' @param expCount Optional dataframe with the categories and expected counts 
#' @param cc c(NULL, "yates", "pearson", "williams") Optional continuity correction (default is NULL)
#' @return dataframe with test statistic, degrees of freedom, p-value, minimum expected count, proportion of expected counts below 5, and test used
#' 
#' @examples
#' data <- c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_g_gof(data)
#' ts_g_gof(data, cc="yates")
#' ts_g_gof(data, cc="pearson")
#' ts_g_gof(data, cc="williams")
#' ts_g_gof(data, eCounts)
#' 
#' @details 
#' The formula used (Wilks, 1938, p. 62):
#' \deqn{G=2\times\sum_{i=1}^{k}\left(F_{i}\times ln\left(\frac{F_{i}}{E_{i}}\right)\right)}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(G,df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^k F_i}
#' If no expected counts provided:
#' \deqn{E_i = \frac{n}{k}}
#' else:
#' \deqn{E_i = n\times\frac{E_{p_i}}{n_p}}
#' \deqn{n_p = \sum_{i=1}^k E_{p_i}}
#' 
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
#' The term ‘Likelihood Ratio Goodness-of-Fit’ can for example be found in an article 
#' from Quine and Robinson (1985), the term ‘Wilks’s likelihood ratio test’ can also be found 
#' in Li and Babu (2019, p. 331), while the term G-test is found in Hoey (2012, p. 4)
#' 
#' The Yates continuity correction (cc="yates") is calculated using (Yates, 1934, p. 222):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i > E_i \\ F_i + 0.5 & \text{ if } F_i < E_i \\ F_i & \text{ if } F_i = E_i \end{cases}}
#' \deqn{G_Y=2\times\sum_{i=1}^{k}\left(F_i^\ast\times ln\left(\frac{F_i^\ast}{E_{i}}\right)\right)}
#' Where if \eqn{F_i^\ast = 0} then \eqn{F_i^\ast\times ln\left(\frac{F_i^\ast}{E_{i}}\right) = 0}
#' 
#' The Pearson correction (cc="pearson") is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{G_{P} = G\times\frac{n - 1}{n}}
#' 
#' The Williams correction (cc="williams") is calculated using (Williams, 1976, p. 36):
#' \deqn{G_{W} = \frac{G}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' The formula is also used by McDonald (2014, p. 87)
#' 
#' **Alternatives**
#' 
#' The library *DescTools* has a similar function *GTest()*
#' 
#' The library *RVAideMemoire* has a similar function *G.test()*
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Hoey, J. (2012). The two-way likelihood ratio (G) test and comparison to two-way chi squared test. 1–6. https://doi.org/10.48550/ARXIV.1206.4881
#' 
#' Li, B., & Babu, G. J. (2019). *A graduate course on statistical inference*. Springer.
#' 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Quine, M. P., & Robinson, J. (1985). Efficiencies of chi-square and likelihood Ratio goodness-of-fit tests. *The Annals of Statistics, 13*(2), 727–742. https://doi.org/10.1214/aos/1176349550
#' 
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. *The Annals of Mathematical Statistics, 9*(1), 60–62. https://doi.org/10.1214/aoms/1177732360
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#'  
#' @export
ts_g_gof <- function(data, expCount=NULL, cc = NULL){
  
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
  
  if (!is.null(cc) && cc == "yates"){
    for (i in 1:k){
      if (as.numeric(freq[i, 2] > expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) - 0.5}
      else if (as.numeric(freq[i, 2] < expC[i])){
        freq[i, 2] = as.numeric(freq[i, 2]) + 0.5}
    }
  }
  
  #calculate the chi-square value
  chiVal = 0
  for (i in 1:k){
    chiVal = chiVal + as.numeric(freq[i, 2])*log(as.numeric(freq[i, 2])/expC[i])
  }
  chiVal = chiVal*2
    
  if (!is.null(cc) && cc == "pearson"){
    chiVal = (n - 1) / n * chiVal}
  else if (!is.null(cc) && cc == "williams"){
    chiVal = chiVal / (1 + (k ^ 2 - 1) / (6 * n * (k - 1)))}

  pValue = pchisq(chiVal, df, lower.tail = FALSE)
  
  #Which test was used
  testUsed = "G test of goodness-of-fit"
  if (!is.null(cc) && cc == "pearson"){
    testUsed = paste(testUsed, ", with E. Pearson continuity correction")}
  else if (!is.null(cc) && cc == "williams"){
    testUsed = paste(testUsed, ", with Williams continuity correction")}
  else if (!is.null(cc) && cc == "yates"){
    testUsed = paste(testUsed, ", with Yates continuity correction")}
  
  statistic = chiVal
  
  testResults <- data.frame(statistic, df, pValue, minExp, propBelow5, testUsed)
  
  return (testResults)
  
}