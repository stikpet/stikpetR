#' Freeman-Tukey Test of Independence
#' 
#' @param nom1 A vector with the data of the first variable
#' @param nom2 A vector with the data of the second variable
#' @param cc c(NULL, "pearson", or "williams") Optional continuity correction (default is NULL)
#' @return dataframe with test statistic, degrees of freedom, p-value, minimum expected count, proportion of expected counts below 5, and test used
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
ts_freeman_tukey_ind <- function(nom1, nom2, cc=NULL){
  
  datFr = na.omit(data.frame(nom1, nom2))
  
  obs = table(datFr)
  
  R = rowSums(obs)
  nr = length(R)
  C = colSums(obs)
  nc = length(C)
  n = sum(R)
  
  expCount = matrix(1, nrow=nr, ncol=nc)
  for (j in 1:nc) {
    for (i in 1:nr){
      expCount[i,j] = R[i]*C[j]/n
    }
  }
  
  chiVal = 4*sum((sqrt(obs) - sqrt(expCount))^2)
  
  if (is.null(cc)==FALSE){
    
    if (cc=="pearson") {
      #Pearson Correction
      chiVal = (n - 1)/n * chiVal
    }
    else if (cc=="williams"){
      #Williams Correction
      q = 1 + (n*sum(1/R)-1) * (n*sum(1/C)-1) / (6*n*(nr - 1)*(nc - 1))
      chiVal = chiVal/q
    }
  }  
  df = (nr - 1)*(nc - 1)
  
  pValue = 1 - pchisq(chiVal, df)
  
  minExp = min(expCount)
  propBelow5 = sum(expCount < 5)/(nr*nc)
  
  statistic = chiVal
  testResults <- data.frame(statistic, df, pValue, minExp, propBelow5)
  
  return(testResults)
  
}