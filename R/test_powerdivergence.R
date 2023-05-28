#' Power Divergence Tests
#' 
#' @param var1 A vector with the data
#' @param var2 Optional vector with data for tests of independence 
#' @param expCounts Optional counts according to null hypothesis
#' @param lambd Optional either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3
#' @param corr Optional correction to be used.
#' @returns 
#' Dataframe with:
#' \item{statistic}{the chi-square statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{two-sided p-value}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of expected counts below 5}
#' \item{testUsed}{a description of the test used}
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all 
#' the categories are equal (the null hypothesis), or with two nominal variables to test if they are 
#' independent.
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square 
#' test (\eqn{\chi^2}), but also an exact multinomial, G-test (\eqn{G^2}), Freeman-Tukey (\eqn{T^2}), 
#' Neyman (\eqn{NM^2}), Mod-Log Likelihood (\eqn{GM^2}), and Freeman-Tukey-Read test are possible.
#' 
#' Cressie and Read (1984, p. 463) noticed how the \eqn{\chi^2}, \eqn{G^2}, \eqn{T^2}, \eqn{NM^2} 
#' and \eqn{GM^2} can all be captured with one general formula. The additional variable lambda 
#' (\eqn{\lambda}) was then investigated, and they settled on a \eqn{\lambda} of 2/3.
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
#' @details 
#' The formula used is (Cressie & Read, 1984, p. 442):
#' \deqn{\chi_{C}^{2} = \begin{cases} 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(F_{i,j}\times ln\left(\frac{F_{i,j}}{E_{i,j}}\right)\right) & \text{ if } \lambda=0 \\ 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(E_{i,j}\times ln\left(\frac{E_{i,j}}{F_{i,j}}\right)\right) & \text{ if } \lambda=-1 \\ \frac{2}{\lambda\times\left(\lambda + 1\right)} \times \sum_{i=1}^{r}\sum_{j=1}^{c} F_{i,j}\times\left(\left(\frac{F_{i,j}}{E_{i,j}}\right)^{\lambda} - 1\right) & \text{ else } \end{cases}}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{C}^{2},df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^r F_{i,j}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{r} the number of categories in the first variable (the number of rows)
#' \item \eqn{c} the number of categories in the second variable (the number of columns)
#' \item \eqn{F_{i,j}} the observed count in row i and column j
#' \item \eqn{E_{i,j}} the expected count in row i and column j
#' \item \eqn{R_i} the i-th row total
#' \item \eqn{C_j} the j-th column total
#' \item \eqn{n} the sum of all counts
#' \item \eqn{\chi^2\left(\dots\right)}	the chi-square cumulative density function
#' }
#' 
#' Cressie and Read (1984, p. 463) suggest to use \eqn{\lambda = \frac{2}{3}},  which
#' is therefor the default in this function.
#' 
#' The **Pearson chi-square statistic** can be obtained by setting \eqn{\lambda = 1}. Pearson's original 
#' formula is (Pearson, 1900, p. 165):
#' \deqn{\chi_{P}^2 = \sum_{i=1}^r \sum_{j=1}^c \frac{\left(F_{i,j} - E_{i,j}\right)^2}{E_{i,j}}}
#' 
#' The **Freeman-Tukey test** has as a formula (Bishop et al., 2007, p. 513):
#' \deqn{T^2 = 4\times\sum_{i=1}^r \sum_{j=1}^c \left(\sqrt{F_{i,j}} - \sqrt{E_{i,j}}\right)^2}
#' This will be same as setting lambda to \eqn{-\frac{1}{2}}. Note that the source for the formula is often quoted to be from Freeman and Tukey (1950) 
#' but couldn't really find it in that article.
#' 
#' **Neyman test** formula was very similar to Pearson's, but the observed and expected counts swapped (Neyman, 1949, p. 250):
#' \deqn{\chi_{N}^2 = \sum_{i=1}^r \sum_{j=1}^c \frac{\left(E_{i,j} - F_{i,j}\right)^2}{F_{i,j}}}
#' This will be same as setting lambda to \eqn{-2}.
#' 
#' The Yates correction (yates) is calculated using (Yates, 1934, p. 222):
#' 
#' Use instead of \eqn{F_{i,j}} the adjusted version defined by:
#' \deqn{F_{i,j}^\ast = \begin{cases} F_{i,j} - 0.5 & \text{ if } F_{i,j}>E_{i,j}  \\ F_{i,j} & \text{ if } F_{i,j}= E_{i,j}\\ F_{i,j} + 0.5 & \text{ if } F_{i,j}<E_{i,j} \end{cases}}
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right)-1\right) \times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_j}\right)-1\right)}{6\times n\times df}}
#' 
#' @examples  
#' nom1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male","male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' ts_powerdivergence(nom1)
#' ts_powerdivergence(nom1, nom2)
#' 
#' @references 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. *Journal of the Royal Statistical Society: Series B (Methodological), 46*(3), 440–464. https://doi.org/10.1111/j.2517-6161.1984.tb01318.x
#' 
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. *The Annals of Mathematical Statistics, 21*(4), 607–611. https://doi.org/10.1214/aoms/1177729756
#' 
#' Neyman, J. (1949). Contribution to the theory of the chi-square test. Berkeley Symposium on Math. Stat, and Prob, 239–273. https://doi.org/10.1525/9780520327016-030
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Pearson, K. (1900). On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. *Philosophical Magazine Series 5, 50*(302), 157–175. https://doi.org/10.1080/14786440009463897
#' 
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. *The Annals of Mathematical Statistics, 9*(1), 60–62. https://doi.org/10.1214/aoms/1177732360
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_powerdivergence <- function(var1, var2=NULL, expCounts=NULL, lambd=c("cressie-read", "g", "mod-log", "freeman-tukey", "neyman"), corr=c("none", "yates", "pearson", "williams")){
  
  #set defaults
  if (length(corr)>1) {corr="none"}
  if (length(lambd)>1) {lambd=2/3}
  
  #Set correction factor to 1 (no correction)
  corFactor = 1
  
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
    testUsed = "power divergence with lambda = " + "lambd"}
  
  #The test itself
  
  if (is.null(var2)){
    #only one variable, so goodness-of-fit version
    
    freqs = table(var1)
    k = length(freqs)
    n = sum(freqs)
    r = k
    c = 1
    
    #Determine expected counts if not provided
    if (is.null(expCounts)){
      expCounts = rep(n/k, k)}
    
    df = k - 1
    
    #set williams correction factor
    if (corr=="williams"){
      corFactor = 1/(1 + (k**2 - 1)/(6*n*df))
      testUsed = testUsed + ", and Williams correction"}
    
    #adjust frequencies if Yates correction is requested
    if (corr=="yates"){
      adjFreq = freqs
      for (i in 1:k){
        if (adjFreq[i] > expCounts[i]){
          adjFreq[i] = adjFreq[i] - 0.5}
        else if (adjFreq[i] < expCounts[i]){
          adjFreq[i] = adjFreq[i] + 0.5}
      }
      
      freqs = adjFreq
      testUsed = testUsed + ", and Yates correction"
    }
  }    
  else{
    #two variables, so use test of independence
    
    #create a cross table
    ct = table(var1, var2)
    rowTotals = rowSums(ct)
    colTotals = colSums(ct)
    
    n = sum(rowTotals)
    r = nrow(ct)
    c = ncol(ct)
    
    df = (r - 1)*(c - 1)
    
    #Determine expected counts if not provided
    if (is.null(expCounts)){
      expCounts = outer(rowTotals, colTotals, '*')/n}
    
    #set williams correction factor
    if (corr=="williams"){
      corFactor = 1/(1 + (n*sum(1/rowTotals) - 1)*(n*sum(1/colTotals) - 1)/(6*n*(r - 1)*(c - 1)))
      testUsed = testUsed + ", and Williams correction"}
    
    #adjust frequencies if Yates correction is requested
    if (corr=="yates"){
      adjFreq = ct
      for (i in 1:r){
        for (j in 1:c){
          if (adjFreq[i,j] > expCounts[i,j]){
            adjFreq[i,j] = adjFreq[i,j] - 0.5}
          else if (adjFreq[i,j] < expCounts[i,j]){
            adjFreq[i,j] = adjFreq[i,j] + 0.5}
        }
      }
      ct = adjFreq
      testUsed = testUsed + ", and Yates correction"
    }
    
    freqs = ct
  }
  
  #determine the test statistic
  if (lambd==0){
    ts = 2*sum(freqs*log(freqs/expCounts))}
  else if (lambd==-1){
    ts = 2*sum(expCounts*log(expCounts/freqs))}
  else{
    ts = 2*sum(freqs*((freqs/expCounts)**(lambd) - 1))/(lambd*(lambd + 1))}
  
  #set E.S. Pearson correction
  if (corr=="pearson"){
    corFactor = (n - 1)/n
    testUsed = testUsed + ", and Pearson correction"}
  
  #Adjust test statistic
  ts = ts*corFactor
  
  #Determine p-value
  pVal = pchisq(ts, df, lower.tail = FALSE)
  
  #Check minimum expected counts
  #Cells with expected count less than 5
  nbelow = length(which(expCounts<5))
  #Number of cells
  ncells = r*c
  #As proportion
  pBelow = nbelow/ncells
  #the minimum expected count
  minExp = min(expCounts)
  
  
  #prepare results
  testResults = data.frame(ts, df, pVal, minExp=minExp, percBelow5=pBelow*100, testUsed)        
  
  return (testResults)
}