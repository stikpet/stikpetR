#' Stuart-Maxwell / Marginal Homogeneity Test
#' 
#' @param nom1 the scores on the first variable
#' @param nom2 the scores on the second variable
#' @return dataframe with the test statistic, degrees of freedom, and p-value (sig.)
#' 
#' @details 
#' The formula used is:
#' \deqn{\chi_{SM}^2 = n\times d' \times S^{-1}\times d}
#' \deqn{df = r - 1 = c - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{SM}\right)}
#' 
#' With:
#' \deqn{S_{i,i} = p_{i,.} + p_{.,i} - 2\times p_{i,i}}
#' \deqn{S_{i,j} = -\left(p_{i,j} + p_{j,i}\right)}
#' \deqn{d_i = p_{i,.} - p_{.,i}}
#' \deqn{p_{i,j} = \frac{F_{i,j}}{n}}
#' \deqn{d = \begin{bmatrix} d_1 \\ d_2 \\ \dots \\ d_{r-1} \end{bmatrix}}
#' \deqn{S = \begin{bmatrix} S_{1,1} & S_{1,2} & \dots & S_{1,c-1} \\ S_{2,1} & S_{2,2} & \dots & S_{2,c-1} \\ \dots & \dots & \dots & \dots \\ S_{r-1,1} & S_{r-1,2} & \dots & S_{r-1,c-1} \\ \end{bmatrix}}
#' \deqn{n = \sum_{i=1}^r\sum_{j=1}^c F_{i,j}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{r} is the number of rows (categories in the first variable)
#' \item \eqn{c} is the number of columns (categories in the second variable)
#' \item \eqn{n} is the total number of scores
#' \item \eqn{F_{i,j}} is the frequency (count) of scores equal to the i-th category in the first variable, and the j-th category in the second.
#' \item \eqn{p_{i,.}} The sum of the proportions in row i
#' \item \eqn{p_{.,i}} The sum of the proportions in column i
#' \item \eqn{d'} is the transpose of the d vector
#' \item \eqn{S^{-1}} is the inverse of the S matrix.
#' \item \eqn{\chi^2\left(\dots, \dots\right)} is the cumulative distribution function of the chi-square distribution
#' }
#' 
#' Note that the d vector and S matrix are one row (and column) less.
#' 
#' This test only differs from the Bhapkar test in the calculation of S
#' 	
#' The procedure is found ins Maxwell (1970) which is founded by Stuart (1955).
#' 
#' @references 
#' Maxwell, A. E. (1970). Comparing the classification of subjects by two independent judges. *The British Journal of Psychiatry, 116*(535), 651–655. https://doi.org/10.1192/bjp.116.535.651
#' 
#' Stuart, A. (1955). A test for homogeneity of the marginal distributions in a two-way classification. *Biometrika, 42*(3/4), 412–416. https://doi.org/10.2307/2333387
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' nom1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
#' nom2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
#' ts_stuart_maxwell(nom1, nom2)
#' 
#' @export
ts_stuart_maxwell <- function(nom1, nom2){
  datFrame = na.omit(data.frame(nom1, nom2))
  
  ct = table(datFrame$nom1, datFrame$nom2)
  r = nrow(ct)
  c = ncol(ct)
  
  n = sum(ct)
  pct = ct/n
  RS = rowSums(pct)
  CS = colSums(pct)
  d = (RS - CS)[1:(r - 1)]
  
  S = matrix(0, nrow=r-1, ncol=c-1)
  for (i in 1:(r-1)) {
    for (j in 1:(c-1)) {
      if (i==j) {
        S[i,j] = RS[i] + CS[j] - 2*pct[i,j]
      }
      else{
        S[i,j] = -(pct[i,j] + pct[j,i])
      }
    }
    
  }
  
  Sinv = solve(S)
  
  dT = t(d)
  
  chi2Value = (n*dT %*% Sinv %*% d)[1]
  df = r - 1
  pValue = pchisq(chi2Value, df, lower.tail = FALSE)
  
  statistic = chi2Value  
  results = data.frame(statistic, df, pValue)
  
  return(results)
  
}