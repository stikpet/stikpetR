#' Neyman Test of Independence
#' @description
#' This test is similar as a Pearson Chi-Square test of independence. If the significance of this test is below 0.05 (or another pre-defined threshold), the two nominal variables have a significant association.
#' 
#' The test compares the observed counts of the cross table with the so-called expected counts. The expected values are the number of respondents you would expect if the two variables would be independent. See the Pearson Chi-Square test of independence for more details on expected counts.
#' 
#' One problem though is that the test should only be used if not too many cells have a so-called expected count, of less than 5, and the minimum expected count is at least 1. So you will also have to check first if these conditions are met. Most often 'not too many cells' is fixed at no more than 20% of the cells. This is often referred to as 'Cochran conditions', after Cochran (1954, p. 420). Note that for example Fisher (1925, p. 83) is more strict, and finds that all cells should have an expected count of at least 5.
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
#' 
#' @details 
#' The formula used is (Neyman, 1949, p. 250):
#' \deqn{\chi_{N}^{2}=\sum_{i=1}^r \sum_{j=1}^c\frac{\left(F_{i,j}-E_{i,j}\right)^{2}}{F_{i,j}}}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{N}^{2},df\right)}
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
#' Cochran, W. G. (1954). Some methods for strengthening the common \eqn{\chi^2} tests. *Biometrics, 10*(4), 417. doi:10.2307/3001616
#' 
#' Fisher, R. A. (1925). *Statistical methods for research workers*. Oliver and Boyd.
#' 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Neyman, J. (1949). Contribution to the theory of the chi-square test. Berkeley Symposium on Math. Stat, and Prob, 239-273. doi:10.1525/9780520327016-030
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 x 2 table. *Biometrika, 34*(1/2), 139-167. doi:10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33-37. doi:10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217-235. doi:10.2307/2983604
#'  
#' @export
ts_neyman_ind <- function(field1, field2, categories1=NULL, categories2=NULL, cc=NULL){
  
  if (is.null(cc)){cc = "none"}
  
  testUsed = "Neyman test of independence"
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
      
      chi2Val = chi2Val + (ct[i, j] - expC[i, j])**2 / ct[i, j]
      
      #check if below 5
      if (expMin < 0 || expC[i,j] < expMin){
        expMin = expC[i,j]}
      if (expC[i,j] < 5){
        nExpBelow5 = nExpBelow5 + 1}
    }
  }
  
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



