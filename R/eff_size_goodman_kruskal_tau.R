#' Goodman-Kruskal tau
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param categories1 optional, categories to use for field1
#' @param categories2 optional, categories to use for field2
#' 
#' @return dataframe with the effect size value, the test statistic, degrees of freedom, and p-value
#' 
#' @details 
#' The formula used is (Goodman & Kruskal, 1954, p. 759):
#' \deqn{\tau_{Y|X} = \frac{n\times\sum_{i=1}^r\sum_{j=1}^c\frac{F_{i,j}^2}{R_i} - \sum_{j=1}^c C_j^2}{n^2 - \sum_{j=1}^c C_j^2}}
#' \deqn{sig. = 1 - \chi^2\left(\left(n - 1\right)\times\left(c - 1\right)\times\tau_{Y|X}, df\right)}
#' \deqn{\tau_{X|Y} = \frac{n\times\sum_{i=1}^r\sum_{j=1}^c\frac{F_{i,j}^2}{C_i} - \sum_{i=1}^r R_i^2}{n^2 - \sum_{i=1}^r R_i^2}}
#' \deqn{sig. = 1 - \chi^2\left(\left(n - 1\right)\times\left(r - 1\right)\times\tau_{X|Y}, df\right)}
#' 
#' With:
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{F_{i,j}} the absolute frequency (observed count) from row i and column j.
#' \item \eqn{c} the number of columns
#' \item \eqn{r} the number of rows
#' \item \eqn{R_i} row total of row i, it can be calculated using \eqn{R_i=\sum_{j=1}^c F_{i,j}}
#' \item \eqn{C_j} column total of column j, it can be calculated using \eqn{C_j=\sum_{i=1}^r F_{i,j}}
#' \item \eqn{n} the total number of cases, it can be calculated in various ways, \eqn{n = \sum_{j=1}^c C_j =\sum_{i=1}^r R_i = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \item \eqn{\chi\left(\dots, \dots\right)} the cumulative density function of the chi-square distribution
#' }
#' 
#' Light and Margolin developed a R2 measure for categorical data, 
#' they proposed a test CATANOVA (Categorical Anova) for this measure. 
#' This was a chi-square test (p. 538). Sarndal (1974, p. 178) concluded that 
#' R2 from Light and Mangolin, was the same as Goodman-Kendal tau, and uses their 
#' test for tau. Margolin and Light (1974) reach the same conclusion and proof the 
#' equivelance.
#' 
#' @references 
#' Goodman, L. A., & Kruskal, W. H. (1954). Measures of association for cross classifications. *Journal of the American Statistical Association, 49*(268), 732–764. https://doi.org/10.2307/2281536
#' 
#' Light, R. J., & Margolin, B. H. (1971). An analysis of variance for categorical data. *Journal of the American Statistical Association, 66*(335), 534–544. https://doi.org/10.1080/01621459.1971.10482297
#' 
#' Margolin, B. H., & Light, R. J. (1974). An analysis of variance for categorical data, II: Small sample comparisons with chi square and other competitors. *Journal of the American Statistical Association, 69*(347), 755–764. https://doi.org/10.1080/01621459.1974.10480201
#' 
#' Särndal, C. E. (1974). A comparative study of association measures. *Psychometrika, 39*(2), 165–187. https://doi.org/10.1007/BF02291467
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' 
#' 
#' @export
es_goodman_kruskal_tau <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #remove if not in categories
  if (!is.null(categories1)){
    field1[! field1 %in% categories1] = NA
  }
  if (!is.null(categories2)){
    field2[! field2 %in% categories2] = NA
  }
  
  
  dFra = na.omit(data.frame(field2, field1))
  
  ct = table(dFra)
  
  r = nrow(ct)
  c = ncol(ct)
  
  Rs = unname(rowSums(ct))
  Cs = unname(colSums(ct))
  n = sum(Rs)
  
  t1 = 0
  tr1 = 0
  for (i in 1:r) {
    for (j in 1:c) {
      t1 = t1 + ct[i,j]**2/Rs[i]
      tr1 = tr1 + ct[i,j]**2/Cs[j]  
    }
    
  }
  
  t2 = sum(Cs**2)
  tr2 = sum(Rs**2)
  
  tau = (n*t1 - t2)/(n**2  - t2)
  tau2= (n*tr1 - tr2)/(n**2  - tr2)
  value = c(tau, tau2)
  
  chi2Value = (n - 1)*(c - 1)*tau
  chi2Value2 = (n - 1)*(r - 1)*tau2  
  statistic = c(chi2Value, chi2Value2)  
  
  df = (r - 1)*(c - 1)
  
  pValue = 1 - pchisq(chi2Value, df)
  pValue2 = 1 - pchisq(chi2Value2, df)
  p = c(pValue, pValue2)
  
  dependent = c("field1", "field2")
  results = data.frame(dependent, value, statistic, df, p)
  
  return(results)
  
}