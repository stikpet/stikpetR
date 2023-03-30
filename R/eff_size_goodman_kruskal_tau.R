#' Goodman-Kruskal tau
#' 
#' @param nom1 the scores on the first variable
#' @param nom2 the scores on the second variable
#' @param dir c("rows", "columns") optional to select which lambda to return (rows = first variable)
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
#' @examples 
#' nom1 <- c("Fully Disagree", "Disagree", "Fully agree", "Neither disagree nor agree", "Agree", "Agree", "Neither disagree nor agree", "Disagree", "Agree", "Agree", "Agree", "Agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Fully agree", "Fully agree", "Fully Disagree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Disagree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully Disagree", "Fully agree", "Agree", "Agree", "Neither disagree nor agree", "Disagree", "Neither disagree nor agree", "Fully agree", "Fully agree", "Disagree", "Disagree", "Neither disagree nor agree", "Disagree", "Agree", "Disagree", "Fully agree", "Fully agree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Fully Disagree")
#' nom2 <- c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' es_goodman_kruskal_tau(nom1, nom2, dir="rows")
#' es_goodman_kruskal_tau(nom1, nom2, dir="columns")
#' 
#' @export
es_goodman_kruskal_tau <- function(nom1, nom2, dir="rows"){
  
  if (dir=="columns") {
    dFra = na.omit(data.frame(nom2, nom1))
    
  }
  else {
    dFra = na.omit(data.frame(nom1, nom2))
    
  }
  ct = table(dFra)
  
  r = nrow(ct)
  c = ncol(ct)
  
  Rs = unname(rowSums(ct))
  Cs = unname(colSums(ct))
  n = sum(Rs)
  
  t1 = 0
  for (i in 1:r) {
    for (j in 1:c) {
      t1 = t1 + ct[i,j]**2/Rs[i]
    }
    
  }
  
  t2 = sum(Cs**2)
  
  tau = (n*t1 - t2)/(n**2  - t2)
  
  chi2Value = (n - 1)*(c - 1)*tau
  df = (r - 1)*(c - 1)
  
  pValue = 1 - pchisq(chi2Value, df)
  
  statistic = chi2Value
  results = data.frame(tau, statistic, df, pValue)
  
  return(results)
  
}