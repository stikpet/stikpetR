#' Bhapkar Test
#' @description
#' If you are only interested if the overall distribution changed (i.e. if the percentages from each category changed or not), you can perform a marginal homogeneity test. There are two that seem to be quite popular for this, the Stuart-Maxwell test (Stuart, 1955; Maxwell, 1970), and the Bhapkar test (Bhapkar, 1961; 1966). According Uebersax (2006) (which also has a nice example) the Bhapkar one is preferred.
#' 
#' Simply put, a marginal homogeneity test, looks at the row vs column proportions. Since in a paired test, the options are the same, if the row and column proportions are the same, nothing changed between the two variables.
#' 
#' @param field1 vector, the first categorical field
#' @param field2 vector, the first categorical field
#' @param categories vector, optional, order and/or selection for categories of field1 and field2
#' 
#' @returns
#' Dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{the chi-squared value}
#' \item{df}{the degrees of freedom used in the test}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used is:
#' \deqn{\chi_{B}^2 = n\times d^{\prime} \times S^{-1} \times d}
#' 
#' With:
#' \deqn{S_{i,i} = p_{i,.} + p_{.,i} - 2\times p_{i,i} - \left(p_{i,.} - p_{.,i}\right)^2}
#' \deqn{S_{i,j} = -\left(p_{i,j} + p_{j,i}\right) - \left(p_{i,.} - p_{.,i}\right)\times\left(p_{j,.} - p_{.,j}\right)}
#' \deqn{d_i = p_{i,.} - p_{.,i}}
#' \deqn{p_{i,j} = \frac{F_{i,j}}{n}}
#' \deqn{d = \begin{bmatrix} d_1 \\ d_2 \\ \dots \\ d_{r-1} \ \end{bmatrix}}
#' \deqn{S = \begin{bmatrix} S_{1,1} & S_{1,2} & \dots & S_{1,c-1} \\ S_{2,1} & S_{2,2} & \dots & S_{2,c-1} \\  \dots & \dots & \dots & \dots \\ S_{r-1,1} & S_{r-1,2} & \dots & S_{r-1,c-1} \\ \end{bmatrix}}
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' 
#' The p-value (sig.):
#' \deqn{df = r - 1 = c - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_B^2, df\right)}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{r}, is the number of rows (categories in the first variable)
#' \item \eqn{c}, is the number of columns (categories in the second variable)
#' \item \eqn{n}, is the total number of scores
#' \item \eqn{F_{i,j}}, is the frequency (count) of scores equal to the i-th category in the first variable, and the j-th category in the second.
#' \item \eqn{p_{i,.}}, The sum of the proportions in row i
#' \item \eqn{p_{.,i}}, The sum of the proportions in column i
#' \item \eqn{d^{\prime}}, is the transpose of the d vector
#' \item \eqn{S^{-1}}, is the inverse of the S matrix.
#' \item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function for the chi-square distribution.
#' }
#' 
#' *Note*
#' 
#' \itemize{
#' \item The d vector and S matrix are one row (and column) less.
#' \item This test only differs from the Stuart-Maxwell test in the calculation of S
#' \item The test was introduced by Bhapkar (1961, 1966)
#' }
#' 
#' @references 
#' Bhapkar, V. P. (1961). Some tests for categorical data. *The Annals of Mathematical Statistics, 32*(1), 72-83. doi:10.1214/aoms/1177705140
#' 
#' Bhapkar, V. P. (1966). A note on the equivalence of two test criteria for hypotheses in categorical data. *Journal of the American Statistical Association, 61*(313), 228-235. doi:10.1080/01621459.1966.10502021
#' 
#' Maxwell, A. E. (1970). Comparing the classification of subjects by two independent judges. *The British Journal of Psychiatry, 116*(535), 651-655. doi:10.1192/bjp.116.535.651
#' 
#' Stuart, A. (1955). A test for homogeneity of the marginal distributions in a two-way classification. *Biometrika, 42*(3/4), 412-416. doi:10.2307/2333387
#' 
#' Uebersax, J. (2006, August 30). McNemar tests of marginal homogeneity. http://www.john-uebersax.com/stat/mcnemar.htm
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_bhapkar <- function(field1, field2, categories=NULL){
  
  #create the cross table
  ct = tab_cross(field1, field2, categories, categories, totals="include")
  
  #basic counts
  k = nrow(ct)-1
  n = ct[k+1,k+1]
  
  #STEP 1: Convert to percentages based on grand total
  p = ct/n
  
  #STEP 2: Determine the differences between the row and the column totals
  d = c()
  for (i in 1:(k-1)){
    d[i] = p[i,k+1] - p[k+1,i]}
  
  #STEP 3: Create the variance and covariance matrix
  #For values on the diagonal add the row and column p total,
  #subtract twice the cell p and then
  #subtract the squared difference between the row and column p total.
  
  #For values not on the diagonal add the mirrored cell p and then
  #add a minus sign, then subtract the product of
  #the difference of the row p total of the current cell and the column p total of the mirrored cell,
  #with the difference of the row p total of the mirrored cell and column p total of the current cell.
  S = matrix(0, nrow=k-1, ncol=k-1)
  for (i in 1:k-1) {
    for (j in 1:k-1) {
      if (i==j) {
        S[i,j] = p[i,k+1] + p[k+1,j] - 2*p[i,j] - (p[i,k+1] - p[k+1,j])**2
      }
      else{
        S[i,j] = -(p[i,j] + p[j,i])- (p[i,k+1] - p[k+1,i])*(p[j,k+1] - p[k+1,j])
      }
    }
    
  }
  
  Sinv = solve(S)
  
  dT = t(d)
  
  chi2Value = (n*dT %*% Sinv %*% d)[1]
  
  df = k - 1
  pValue = pchisq(chi2Value, df, lower.tail = FALSE)
  
  statistic=chi2Value
  res = data.frame(n, statistic, df, pValue)
  colnames(res) = c("n", "statistic", "df", "p-value")
  return(res)
  
}



