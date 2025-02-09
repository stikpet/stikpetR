#' Power Divergence Goodness-of Fit Tests
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories are equal (the null hypothesis)
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test (\eqn{\chi^2}), but also an exact multinomial, G-test (\eqn{G^2}), Freeman-Tukey (\eqn{T^2}), Neyman (\eqn{NM^2}), Mod-Log Likelihood (\eqn{GM^2}), and Freeman-Tukey-Read test are possible.
#' 
#' Cressie and Read (1984, p. 463) noticed how the \eqn{\chi^2}, \eqn{G^2}, \eqn{T^2}, \eqn{NM^2} and \eqn{GM^2} can all be captured with one general formula. The additional variable lambda (\eqn{\lambda}) was then investigated, and they settled on a \eqn{\lambda} of 2/3.
#' 
#' By setting \eqn{\lambda} to different values, we get the different tests:
#' \itemize{
#' \item{\eqn{\lambda = 1}}{Pearson chi-square}
#' \item{\eqn{\lambda = 0}}{G/Wilks/Likelihood-Ratio}
#' \item{\eqn{\lambda = -\frac{1}{2}}}{Freeman-Tukey}
#' \item{\eqn{\lambda = -1}}{Mod-Log-Likelihood}
#' \item{\eqn{\lambda = -2}}{Neyman}
#' \item{\eqn{\lambda = \frac{2}{3}}}{Cressie-Read}
#' }
#' 
#' @param data A vector or dataframe with the data
#' @param expCounts Optional dataframe with the categories and expected counts
#' @param lambd Optional either name of test or specific value. Either "cressie-read" (default), "g", "mod-log", "freeman-tukey", or "neyman"
#' @param cc Optional continuity correction. Either "none" (default), "yates", "pearson", or "williams"
#' 
#' @returns 
#' Dataframe with:
#' \item{statistic}{the chi-square statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{two-sided p-value}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of expected counts below 5}
#' \item{testUsed}{a description of the test used}
#' 
#' 
#' @details 
#' The formula used is (Cressie & Read, 1984, p. 442):
#' \deqn{\chi_{C}^{2} = \begin{cases} 2\times\sum_{i=1}^{k}F_{i}\times ln\left(\frac{F_{i}}{E_{i}}\right) & \text{ if } \lambda=0 \\ 2\times\sum_{i=1}^{k} E_{i}\times ln\left(\frac{E_{i}}{F_{i}}\right) & \text{ if } \lambda=-1 \\ \frac{2}{\lambda\times\left(\lambda + 1\right)} \times \sum_{i=1}^{k} F_{i}\times\left(\left(\frac{F_{i}}{E_{i}}\right)^{\lambda} - 1\right) & \text{ else } \end{cases}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{C}^{2},df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \deqn{E_{i} = \frac{n}{k}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{k} the number of categories
#' \item \eqn{F_{i}} the observed count of category i
#' \item \eqn{E_{i}} the expected count of category i
#' \item \eqn{n} the sum of all counts
#' \item \eqn{\chi^2\left(\dots\right)}	the chi-square cumulative density function
#' }
#' 
#' Cressie and Read (1984, p. 463) suggest to use \eqn{\lambda = \frac{2}{3}},  which
#' is therefor the default in this function.
#' 
#' The **Pearson chi-square statistic** can be obtained by setting \eqn{\lambda = 1}. 
#' 
#' The **Freeman-Tukey test** will be same as setting lambda to \eqn{-\frac{1}{2}}. 
#' 
#' **Neyman test** will be same as setting lambda to \eqn{-2}.
#' 
#' The Yates continuity correction (cc="yates") is calculated using (Yates, 1934, p. 222):
#' \deqn{F_i^\ast  = \begin{cases} F_i - 0.5 & \text{ if } F_i > E_i \\ F_i + 0.5 & \text{ if } F_i < E_i \\ F_i & \text{ if } F_i = E_i \end{cases}}
#' 
#' Note that the Yates correction is usually only considered if there are only two categories. Some also argue this correction is too conservative (see for details Haviland (1990)).
#' 
#' The Pearson correction (cc="pearson") is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{adj}^2 = \chi_{C}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (cc="williams") is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{adj}^2 = \frac{\chi_{C}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{k^2 - 1}{6\times n\times df}}
#' The formula is also used by McDonald (2014, p. 87)
#' 
#' @references 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. *Journal of the Royal Statistical Society: Series B (Methodological), 46*(3), 440–464. doi:10.1111/j.2517-6161.1984.tb01318.x
#' 
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. *The Annals of Mathematical Statistics, 21*(4), 607–611. doi:10.1214/aoms/1177729756
#' 
#' Haviland, M. G. (1990). Yates’s correction for continuity and the analysis of 2 × 2 contingency tables. *Statistics in Medicine, 9*(4), 363–367. doi:10.1002/sim.4780090403
#' 
#' Neyman, J. (1949). Contribution to the theory of the chi-square test. Berkeley Symposium on Math. Stat, and Prob, 239–273. doi:10.1525/9780520327016-030
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. doi:10.2307/2332518
#' 
#' Pearson, K. (1900). On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. *Philosophical Magazine Series 5, 50*(302), 157–175. doi:10.1080/14786440009463897
#' 
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. *The Annals of Mathematical Statistics, 9*(1), 60–62. doi:10.1214/aoms/1177732360
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
#' ts_powerdivergence_gof(ex1)
#' 
#' #Example 2: dataframe with various settings
#' ex2 = df1['mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_powerdivergence_gof(ex2, expCounts=eCounts)
#' ts_powerdivergence_gof(ex2, expCounts=eCounts, cc="pearson")
#' ts_powerdivergence_gof(ex2, expCounts=eCounts, cc="williams")
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_powerdivergence_gof(ex3)
#' 
#' @export
ts_powerdivergence_gof <- function(data, expCounts=NULL, lambd=c("cressie-read", "g", "mod-log", "freeman-tukey", "neyman"), cc=c("none", "yates", "pearson", "williams")){
  
  #set defaults
  if (length(cc)>1) {cc="none"}
  if (length(lambd)>1) {lambd=2/3}
  
  #Set correction factor to 1 (no correction)
  corFactor = 1
  
  data = data.frame(data)
  data = na.omit(data)
  
  #Test Used
  if (lambd == 2/3 || lambd == "cressie-read"){
    lambd = 2/3
    testUsed = "Cressie-Read"}
  else if (lambd==0 || lambd == "g"){
    lambd=0
    testUsed = "likelihood-ratio"}    
  else if (lambd==-1 || lambd == "mod-log"){
    lambd=-1
    testUsed = "mod-log likelihood ratio"}    
  else if (lambd==1 || lambd=="pearson"){
    lambd=1
    testUsed = "Pearson chi-square"}
  else if (lambd==-0.5 || lambd=="freeman-tukey"){
    lambd=-0.5
    testUsed = "Freeman-Tukey"}        
  else if (lambd==-2 || lambd=="neyman"){
    lambd==-2
    testUsed = "Neyman"}
  else {
    testUsed = paste("power divergence with lambda = ", lambd)}
  
  #The test itself
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
  
  #set williams correction factor
  if (cc=="williams"){
    corFactor = 1/(1 + (k**2 - 1)/(6*n*df))
    testUsed = paste(testUsed, ", and Williams correction")}
  
  #adjust frequencies if Yates correction is requested
  if (cc=="yates"){
    adjFreq = freq
    for (i in 1:k){
      if (adjFreq[i,2] > expCounts[i,2]){
        adjFreq[i,2] = adjFreq[i,2] - 0.5}
      else if (adjFreq[i,2] < expCounts[i,2]){
        adjFreq[i,2] = adjFreq[i,2] + 0.5}
    }
    
    freq = adjFreq
    testUsed = paste(testUsed, ", and Yates correction")
  }
  
  freqs = as.numeric(freq[,2])
  #determine the test statistic
  if (lambd==0){
    ts = 2*sum(freqs*log(freqs/expC))}
  else if (lambd==-1){
    ts = 2*sum(expC*log(expC/freqs))}
  else{
    ts = 2*sum(freqs*((freqs/expC)**(lambd) - 1))/(lambd*(lambd + 1))}
  
  #set E.S. Pearson correction
  if (cc=="pearson"){
    corFactor = (n - 1)/n
    testUsed = paste(testUsed, ", and Pearson correction")}
  
  #Adjust test statistic
  ts = ts*corFactor
  
  #Determine p-value
  pVal = pchisq(ts, df, lower.tail = FALSE)
  
  #Check minimum expected counts
  #Cells with expected count less than 5
  nbelow = length(which(expCounts<5))
  #Number of cells
  ncells = k
  #As proportion
  pBelow = nbelow/ncells
  #the minimum expected count
  minExp = min(expC)
  
  
  #prepare results
  testResults = data.frame(n, k, ts, df, pVal, minExp=minExp, percBelow5=pBelow*100, testUsed)
  colnames(testResults)<-c("n", "k", "statistic", "df", "p-value", "minExp", "propBelow5", "test")
  
  return (testResults)
}