#' Post-Hoc Dunn Test (for Cochran Q test)
#' @description 
#' An adaptation from IBM SPSS on the Dunn test, so it can be used as a post-hoc test for a Cochran Q test.
#' 
#' @param data dataframe with the binary scores
#' @param success indicator for what is considered a success (default is first value found)
#' 
#' @returns 
#' A dataframe with:
#' \item{category 1}{label of first variable in comparison}
#' \item{category 2}{label of second variable in comparison}
#' \item{n suc. 1}{number of successes in first variable in comparison}
#' \item{n suc. 2*}{number of successes in second variable in comparison}
#' \item{statistic}{test statistic}
#' \item{z-value}{standardized test statistic (z-value)}
#' \item{p-value}{p-value of the z-value}
#' \item{adj. p-value}{Bonferroni corrected p-value}
#' 
#' @details 
#' The formula used (IBM, 2021, p. 814):
#' \deqn{z_{1,2} = \frac{\bar{d}_{1,2}}{SE}}
#' \deqn{sig. = 2\times\left(1-\Phi\left(\left|z_{1,2}\right|\right)\right)}
#' 
#' With:
#' \deqn{\bar{d}_{1,2} = \frac{ns_1 - ns_2}{n}}
#' \deqn{SE = \sqrt{2\times\frac{k\times\sum_{i=1}^n R_i - \sum_{i=1}^n R_i^2}{n^2\times k\times\left(k-1\right)}}}
#' \deqn{R_i = \sum_{j=1}^k s_{i,j}}
#' \deqn{ns_j = \sum_{i=1}^n s_{i,j}}    
#' \deqn{s_{i,j} = \begin{cases} 1 & \text{ if } x_{i,j}= \text{success} \\ 0 & \text{ if } x_{i,j} \neq \text{success} \end{cases}}
#' 
#' IBM SPSS mentions this is an adaptation from Dunn (1964), originally for the Kruskal-Wallis test.
#' 
#' The Bonferroni adjustment is done using:
#' \deqn{sig._{adj} = \min \left(sig. \times n_c, 1\right)}
#' \deqn{n_c = \frac{k\times\left(k-1\right)}{2}}
#' 
#' *Symbols used*
#' 
#' \itemize{
#' \item \eqn{x_{i,j}}, the score in row i and column j
#' \item \eqn{k}, the number of variables
#' \item \eqn{n}, the total number of cases used
#' \item \eqn{ns_j}, the total number of successes in column j
#' \item \eqn{R_i}, the total number of successes in row i
#' \item \eqn{\Phi\left(\dots\right)}, the standard normal cumulative distribution function.
#' \item \eqn{n_c}, the number of comparisons (pairs)
#' }
#' 
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics, 6*(3), 241–252. doi:10.1080/00401706.1964.10490181
#' 
#' IBM. (2021). IBM SPSS Statistics Algorithms. IBM.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_dunn_q <- function(data, success = NULL){
  dFr = na.omit(data)  
  k = ncol(dFr)
  n = nrow(dFr)
  
  if (is.null(success)){
    success=dFr[1,1]}
  
  #row successes
  R = rep(0, n)
  for (i in 1:n) {
    R[i] = sum(dFr[i,]==success)
  }
  rst = sum(R)
  rs2t = sum(R**2)
  
  #standard error
  se = (2 * (k * rst - rs2t) / (k * (k - 1) * n**2))**0.5
  
  #number of comparisons
  ncomp = k * (k - 1) / 2
  
  #the pairwise comparisons
  varNames = colnames(data)
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 8))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      cat1 = varNames[i]
      cat2 = varNames[j]
      selDf = dFr[c(cat1, cat2)]
      n1 = sum(selDf[cat1]==success)
      n2 = sum(selDf[cat2]==success)
      t = (n1 - n2)/n
      z = t/se
      pVal = 2 * (1 - pnorm(abs(z)))
      if (pVal*ncomp > 1){
        pAdj = 1}
      else{
        pAdj = pVal*ncomp}
      
      res[resRow, 1] = cat1
      res[resRow, 2] = cat2
      res[resRow, 3] = n1
      res[resRow, 4] = n2
      res[resRow, 5] = t
      res[resRow, 6] = z
      res[resRow, 7] = pVal
      res[resRow, 8] = pAdj
      
      resRow=resRow+1
    }
  }
  
  colnames(res) = c("category 1", "category 2", "n suc. 1", "n suc. 2", "statistic", "z-value", "p-value", "adj. p-value")
  return (res)
  
}