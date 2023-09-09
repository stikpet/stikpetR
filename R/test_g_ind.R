#' G (Likelihood Ratio / Wilks) Test of Independence
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
#' The formula used (Wilks, 1938, p. 62):
#' \deqn{G=2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(F_{i,j}\times ln\left(\frac{F_{i,j}}{E_{i,j}}\right)\right)}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(G,df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^r F_{i,j}}
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
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. *The Annals of Mathematical Statistics, 9*(1), 60–62. https://doi.org/10.1214/aoms/1177732360
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @examples  
#' nom1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male","male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' ts_g_ind(nom1, nom2)
#' ts_g_ind(nom1, nom2, cc="yates")
#' ts_g_ind(nom1, nom2, cc="pearson")
#' ts_g_ind(nom1, nom2, cc="williams")
#'  
#' @export
ts_g_ind <- function(field1, field2, categories1=NULL, categories2=NULL, cc=NULL){
  
  if (is.null(cc)){cc = "none"}
  
  testUsed = "G test of independence"
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
      
      chi2Val = chi2Val + ct[i, j]*log(ct[i, j]/expC[i, j])
      
      #check if below 5
      if (expMin < 0 || expC[i,j] < expMin){
        expMin = expC[i,j]}
      if (expC[i,j] < 5){
        nExpBelow5 = nExpBelow5 + 1}
    }
  }
  chi2Val= 2*chi2Val
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
