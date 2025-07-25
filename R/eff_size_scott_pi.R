#' Scott Pi
#' @description
#' An effect size meaure, that measures the how strongly two raters or variables, agree with each other. Full agreement would result in a pi of 1.
#' 
#' The measure is very similar to Cohen's kappa. The difference is with the calculation of the expected marginal proportions. Cohen's kappa uses a squared geometric mean, while Scott's pi uses squared arithmetic means. 
#' 
#' Scott developed this in criticism on Bennett-Alpert-Goldstein's S (see es_bag_s()).
#' 
#' @param field1 vector, the first categorical field
#' @param field2 vector, the first categorical field
#' @param categories vector, optional, order and/or selection for categories of field1 and field2
#' 
#' @returns
#' Dataframe with:
#' \item{Scott pi}{the Scott pi value}
#' \item{n}{the sample size}
#' \item{statistic}{the test statistic (z-value)}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Scott, 1955, p. 323):
#' \deqn{\pi = \frac{p_0 - p_e}{1 - p_e}}
#' 
#' With:
#' \deqn{P = \sum_{i=1}^r F_{i,i}}
#' \deqn{p_0 = \frac{P}{n}}
#' \deqn{p_e = \sum_{i=1}^r\left(\frac{R_i + C_i}{2\times n}\right)^2}
#' 
#' The asymptotic standard errors is calculated using (Scott, 1955, p. 325):
#' \deqn{ASE = \sqrt{\left(\frac{1}{1 - p_e}\right)^2\times\frac{p_0\times\left(1 - p_0\right)}{n - 1}}}
#' 
#' The p-value (significance) is then calculated using:
#' \deqn{z_{\pi} = \frac{\pi}{ASE}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(z_{\kappa}\right)\right)}
#' 
#' *Symbols used*
#' \itemize{
#'  \item \eqn{F_{i,j}}, the observed count in row i and column j.
#'  \item \eqn{r}, is the number of rows (categories in the first variable)
#'  \item \eqn{c}, is the number of columns (categories in the second variable)
#'  \item \eqn{n}, is the total number of scores
#'  \item \eqn{R_i}, the row total of row i.\eqn{R_i = \sum_{j=1}^c F_{i,j}}
#'  \item \eqn{C_j}, the column total of column j.\eqn{C_j = \sum_{i=1}^r F_{i,j}}
#' }
#' 
#' @references
#' Scott, W. A. (1955). Reliability of content analysis: The case of nominal scale coding. *The Public Opinion Quarterly, 19*(3), 321-325.
#'  
#' @export
es_scott_pi <- function(field1, field2, categories=NULL){
  #create the cross table
  ct = tab_cross(field1, field2, categories, categories, totals="include")
  
  #basic counts
  k = nrow(ct)-1
  n = ct[k+1, k+1]
  
  #determine p0 and pe
  p0 = 0
  pe = 0
  for (i in 1:k){
    p0 = p0 + ct[i, i]
    pe = pe + ((ct[i, k+1] + ct[k+1, i])/(2*n))**2
  }
  p0 = p0/n
  
  #Scott's Pi
  scottPi = (p0 - pe) / (1 - pe)
  
  #Test
  ase = ((1 / (1 - pe))**2 * p0 * (1 - p0) / (n - 1))**0.5
  z = scottPi / ase
  pValue = 2 * (1 - pnorm(abs(z))) 
  
  res = data.frame(scottPi, n, z, pValue)
  colnames(res) = c("Scott pi", "n", "statistic", "p-value")
  
  return(res)
  
}



