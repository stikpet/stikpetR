#' Power Divergence Test of Independence
#' @description 
#' A test that can be used with two nominal variables to test if they are independent.
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
#' 
#' @param field1 list or dataframe with the first categorical field
#' @param field2 list or dataframe with the second categorical field
#' @param categories1 optional list with order and/or selection for categories of field1
#' @param categories2 optional list with order and/or selection for categories of field2
#' @param cc optional methdod for continuity correction. Either NULL (default), "yates", "pearson", "williams".
#' @param lambd Optional either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3
#' 
#' @returns
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{n rows}{number of categories used in first field}
#' \item{n col.}{number of categories used in second field}
#' \item{statistic}{the test statistic (chi-square value)}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value)}
#' \item{min. exp.}{the minimum expected count}
#' \item{prop. exp. below 5}{proportion of cells with expected count less than 5}
#' \item{test}{description of the test used}
#' 
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
ts_powerdivergence_ind <- function(field1, field2, categories1=NULL, categories2=NULL, cc=NULL, lambd=2/3){
  
  if (is.null(cc)){cc = "none"}
  
  #Test Used
  if (lambd == 2 / 3 || lambd == "cressie-read"){
    lambd = 2 / 3
    testUsed = "Cressie-Read test of independence"}    
  else if (lambd == 0 || lambd == "likelihood-ratio"){
    lambd = 0
    testUsed = "likelihood ratio test of independence"}        
  else if (lambd == -1 || lambd == "mod-log"){
    lambd = -1
    testUsed = "mod-log likelihood ratio test of independence"}        
  else if (lambd == 1 || lambd == "pearson"){
    lambd = 1
    testUsed = "Pearson chi-square test of independence"}    
  else if (lambd == -0.5 || lambd == "freeman-tukey"){
    lambd = -0.5
    testUsed = "Freeman-Tukey test of independence"}        
  else if (lambd == -2 || lambd == "neyman"){
    lambd = -2
    testUsed = "Neyman test of independence"}
  else{
    testUsed = paste("power divergence test of independence with lambda = ", lambd)}  
  
  if (cc == "yates"){
    testUsed = paste(testUsed ,", with Yates continuity correction")}
  
  #create the cross table
  ct = tab_cross(field1, field2, categories1, categories2, totals="include")
  
  #basic counts
  nrows = nrow(ct) - 1
  ncols =  ncol(ct) - 1
  n = ct[nrows+1, ncols+1]
  
  #determine the expected counts & chi-square value
  chi2Val = 0
  expMin = -1
  nExpBelow5 = 0    
  expC = data.frame()
  for (i in 1:nrows){
    for (j in 1:ncols){
      expC[i, j] = ct[nrows+1, j] * ct[i, ncols+1] / n
      
      #add or remove a half in case Yates correction
      if (cc=="yates"){
        if (ct[i,j] > expC[i,j]){
          ct[i,j] = ct[i,j] - 0.5}
        else if (ct[i,j] < expC[i,j]){
          ct[i,j] = ct[i,j] + 0.5}
      }
      
      if (lambd == 0){
        chi2Val = chi2Val + ct[i,j] * log(ct[i,j] / expC[i,j])}
      else if (lambd == -1){
        chi2Val = chi2Val + expC[i,j] * log(expC[i,j] / ct[i,j])}
      else{
        chi2Val = chi2Val + ct[i,j] * ((ct[i,j] / expC[i,j]) ^ lambd - 1)}
      
      
      
      #check if below 5
      if (expMin < 0 || expC[i,j] < expMin){
        expMin = expC[i,j]}
      if (expC[i,j] < 5){
        nExpBelow5 = nExpBelow5 + 1}
    }
  }
  
  if (lambd == 0){
    chi2Val = 2 * chi2Val}
  else if (lambd == -1){
    chi2Val = 2 * chi2Val}
  else{
    chi2Val = 2 / (lambd * (lambd + 1)) * chi2Val}
  
  nExpBelow5 = nExpBelow5/(nrows*ncols)
  
  #Degrees of freedom
  df = (nrows - 1)*(ncols - 1)
  
  #Williams and Pearson correction
  if (cc == "williams"){
    testUsed = paste(testUsed, ", with Williams continuity correction")
    rTotInv = 0
    for (i in 1:nrows){
      rTotInv = rTotInv + 1 / ct[i, ncols+1]}
    
    cTotInv = 0
    for (j in 1:ncols){
      cTotInv = cTotInv + 1 / ct[nrows+1, j]}
    
    q = 1 + (n * rTotInv - 1) * (n * cTotInv - 1) / (6 * n * df)
    chi2Val = chi2Val / q}
  else if (cc == "pearson"){
    testUsed = paste(testUsed, ", with E.S. Pearson continuity correction")
    chi2Val = chi2Val * (n - 1) / n}
  
  #The test
  pvalue = 1 - pchisq(chi2Val, df)
  
  #Prepare the results
  results <- data.frame(n, nrows, ncols, chi2Val, df, pvalue, expMin, nExpBelow5, testUsed)
  colnames(results)<-c("n", "n rows", "n col.", "statistic", "df", "p-value", "min. exp.", "prop. exp. below 5", "test")
  
  return (results)
  
}