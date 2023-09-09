#' Pearson Chi-Square Test of Independence
#' @description
#' To test if two nominal variables have an association, the most commonly used test is the Pearson chi-square test of independence (Pearson, 1900). If the significance of this test is below 0.05, the two nominal variables have a significant association.
#' 
#' The test compares the observed counts of the cross table with the so-called expected counts. The expected values are the number of respondents you would expect if the two variables would be independent.
#' 
#' If for example I had 50 male and 50 female respondents, and 50 agreed with a statement and 50 disagreed with the statement, the expected value for each combination (male-agree, female-agree, male-disagree, and female-disagree) would be 25.
#' 
#' Note that if in the survey the real results would be that all male disagreed, and all female would agree, there is a full dependency (i.e. gender fully decides if you agree or disagree), even though the row and column totals would still be 50. In essence the Pearson chi-square test, checks if your data is more toward the expected values (independence) or the full dependency one.
#' 
#' One problem though is that the Pearson chi-square test should only be used if not too many cells have a so-called expected count, of less than 5, and the minimum expected count is at least 1. So you will also have to check first if these conditions are met. Most often ‘not too many cells’ is fixed at no more than 20% of the cells. This is often referred to as 'Cochran conditions', after Cochran (1954, p. 420). Note that for example Fisher (1925, p. 83) is more strict, and finds that all cells should have an expected count of at least 5 .
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
#'@details
#' The formula used is (Pearson, 1900, p. 165):
#' \deqn{\chi_p^2 = \sum_{i=1}^r \sum_{j=1}^c \frac{\left(F_{i,j} - E_{i,j}\right)^2}{E_{i,j}}}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(\chi_p^2, df\right)}
#' 
#' With:
#' \deqn{E_{i,j} = \frac{R_i \times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^r F_{i,j}}
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j} = \sum_{i=1}^r R_i = \sum_{j=1}^c C_j}
#' 
#' Symbols:
#' \itemize{
#' \item \eqn{r}, the number of rows
#' \item \eqn{c}, the number of columns
#' \item \eqn{F_{i,j}}, the observed count in row i and column j.
#' \item \eqn{E_{i,j}}, the expected count in row i and column j.
#' \item \eqn{R_i}, the row total of row i 
#' \item \eqn{C_j}, the column total of column j
#' \item \eqn{n}, the overall total.
#' \item \eqn{df}, the degrees of freedom
#' }
#' 
#' The **Yates** correction uses \eqn{F_{i,j}'} instead of \eqn{F_{i,j}}, defined as (Yates, 1934, p. 222):
#' \deqn{F_{i,j}' = \begin{cases} F_{i,j}-\frac{1}{2} & \text{ if } F_{i,j}> E_{i,j} \\ F_{i,j}+\frac{1}{2} & \text{ if } F_{i,j}< E_{i,j} \\ F_{i,j} & \text{ if } F_{i,j}= E_{i,j} \end{cases}}
#' 
#' The **Williams** correction, adjusts the Pearson chi-square value:
#' \deqn{\chi_{wil}^2 = \frac{\chi_p^2}{q}}
#' 
#' With:
#' \deqn{q = 1+\frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right) - 1\right)\times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_i}\right) - 1\right)}{6\times n \times\left(r - 1\right)\times\left(c - 1\right)}}
#' 
#' The formula is probably from Williams (1976) but the one shown here is taken from McDonald (1976, p. 36).
#' 
#' The **Pearson** correction also adjusts the Pearson chi-square value with (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{epearson}^2 = \frac{n - 1}{n}\times \chi_p^2}
#' 
#' @references
#' Cochran, W. G. (1954). Some methods for strengthening the common \eqn{\chi^2} tests. *Biometrics, 10*(4), 417. doi:10.2307/3001616
#' 
#' Fisher, R. A. (1925). *Statistical methods for research workers*. Oliver and Boyd.
#' 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. doi:10.2307/2332518
#' 
#' Pearson, K. (1900). On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. *Philosophical Magazine Series 5, 50*(302), 157–175. doi:10.1080/14786440009463897
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. doi:10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. doi:10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_pearson_ind <- function(field1, field2, categories1=NULL, categories2=NULL, cc=NULL){
  
  if (is.null(cc)){cc = "none"}
  
  testUsed = "Pearson chi-square test of independence"
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
      
      chi2Val = chi2Val + (ct[i, j] - expC[i, j])**2 / expC[i, j]
      
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