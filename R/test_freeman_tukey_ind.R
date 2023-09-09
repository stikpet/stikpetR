#' Freeman-Tukey Test of Independence
#' 
#' @param field1 list or dataframe with the first categorical field
#' @param field2 list or dataframe with the second categorical field
#' @param categories1 optional list with order and/or selection for categories of field1
#' @param categories2 optional list with order and/or selection for categories of field2
#' @param cc optional methdod for continuity correction. Either NULL (default), "yates", "pearson", "williams".
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
#' @details 
#' The formula used is (Bishop et al., 2007, p. 513):
#' \deqn{T^2=4\times\sum_{i=1}^r \sum_{j=1}^c \left(\sqrt{F_{i,j}} - \sqrt{E_{i,j}}\right)^2}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(T^2,df\right)}
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
#' The test is attributed to Freeman and Tukey (1950), but couldn't really find it in there.
#' Ayinde and Abidoye (2010) also show the formula in more modern notation.
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right)-1\right) \times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_j}\right)-1\right)}{6\times n\times df}}
#' 	
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Ayinde, K., & Abidoye, A. O. (2010). Simplified Freeman-Tukey test statistics for testing probabilities in contingency tables. *Science World Journal, 2*(2), 21–27. https://doi.org/10.4314/swj.v2i2.51730
#' 
#' Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). *Discrete multivariate analysis*. Springer.
#' 
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations Related to the angular and the square root. *The Annals of Mathematical Statistics, 21*(4), 607–611. https://doi.org/10.1214/aoms/1177729756
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' @examples  
#' nom1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male","male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' ts_freeman_tukey_ind(nom1, nom2)
#' ts_freeman_tukey_ind(nom1, nom2, cc="pearson")
#' ts_freeman_tukey_ind(nom1, nom2, cc="williams")
#'  
#' @export
ts_freeman_tukey_ind <- function(field1, field2, categories1=NULL, categories2=NULL, cc=NULL){
  
  if (is.null(cc)){cc = "none"}
  
  testUsed = "Freeman-Tukey test of independence"
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
      
      chi2Val = chi2Val + (ct[i, j]**0.5 - expC[i, j]**0.5)**2
      
      #check if below 5
      if (expMin < 0 || expC[i,j] < expMin){
        expMin = expC[i,j]}
      if (expC[i,j] < 5){
        nExpBelow5 = nExpBelow5 + 1}
    }
  }
  
  chi2Val = 4*chi2Val  
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