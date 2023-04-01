#' McNemar-Bowker Test
#' 
#' @param nom1 the scores on the first variable
#' @param nom2 the scores on the second variable
#' @return dataframe with the test statistic, degrees of freedom, and p-value (sig.)
#' 
#' @details 
#' The formula used is (Bowker, 1948, p. 573):
#' \deqn{\chi_{B}^2 = \sum_{i=1}^{r-1}\sum_{j=i+1}^c \frac{\left(F_{i,j} - F_{j,i}\right)^2}{F_{i,j} + F_{j,i}}}
#' \deqn{df = \frac{r \times \left(r - 1\right)}{2} = \frac{c \times \left(c - 1\right)}{2}}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{B}\right)}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{r} is the number of rows (categories in the first variable)
#' \item \eqn{c} is the number of columns (categories in the second variable)
#' \item \eqn{n} is the total number of scores
#' \item \eqn{F_{i,j}} is the frequency (count) of scores equal to the i-th category in the first variable, and the j-th category in the second.
#' }
#' 
#' This test is an extension of the McNemar (1947) test, which is only for 2x2 tables
#' 
#' @references 
#' Bowker, A. H. (1948). A test for symmetry in contingency tables. *Journal of the American Statistical Association, 43*(244), 572–574. https://doi.org/10.2307/2280710
#' 
#' McNemar, Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika, 12*(2), 153–157. https://doi.org/10.1007/BF02295996
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
#' ts_mcnemar_bowker(nom1, nom2)
#' 
#' @export
ts_mcnemar_bowker <- function(nom1, nom2){
  datFrame = na.omit(data.frame(nom1, nom2))
  
  ct = table(datFrame$nom1, datFrame$nom2)
  r = nrow(ct)
  c = ncol(ct)
  
  chi2Value = 0
  for (i in 1:r) {
    for (j in 1:c) {
      if (i>j) {
        chi2Value = chi2Value + (ct[i,j] - ct[j,i])**2/(ct[i,j] + ct[j,i])
      }
    }
  }
  
  df = r*(r - 1)/2
  pValue = pchisq(chi2Value, df, lower.tail = FALSE)
  
  statistic = chi2Value
  results = data.frame(statistic, df, pValue)
  
  return(results)
  
}