#' Goodman-Kruskal Gamma
#' 
#' @param scores the numeric scores of the first variable
#' @param groups the numeric scores of the second variable
#' @returns 
#' A dataframe with:
#' \item{g}{the Goodman-Kruskal Gamma value}
#' \item{ASE1}{the asymptotic standard error not assuming null hypothesis}
#' \item{ASE0}{the asymptotic standard error assuming null hypothesis}
#' \item{statistic}{the z-value used for the test}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Goodman & Kruskal, 1954, p. 749):
#' \deqn{\gamma = \frac{P - Q}{P + Q}}
#' \deqn{z_{\gamma} = \frac{\gamma}{ASE_0}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(z_{\gamma}\right)\right)}
#' With:
#' \deqn{P = \sum_{i=1}^r \sum_{j=1}^c P_{i,j}}
#' \deqn{Q = \sum_{i=1}^r \sum_{j=1}^c Q_{i,j}}
#' \deqn{P_{i,j} = F_{i,j}\times C_{i,j}}
#' \deqn{Q_{i,j} = F_{i,j}\times D_{i,j}}
#' \deqn{C_{i,j} = \sum_{h<i}\sum_{k<j} F_{h,k} + \sum_{h>i}\sum_{k>j} F_{h,k}}
#' \deqn{D_{i,j} = \sum_{h<i}\sum_{k>j} F_{h,k} + \sum_{h>i}\sum_{k<j} F_{h,k}}
#' \deqn{ASE_1 = \frac{4}{\left(P + Q\right)^2}\times\sqrt{\sum_{i=1}^r \sum_{j=1}^c F_{i,j}\times\left(Q\times C_{i,j} - P\times D_{i,j}\right)^2}}
#' \deqn{ASE_0 = \frac{2}{P + Q}\times\sqrt{\sum_{i=1}^r \sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2 - \frac{\left(P - Q\right)^2}{n}}}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{F_{i,j}} the count of scores equal to i in the first variable and j in the second
#' \item \eqn{n} the total sample size
#' \item \eqn{r} the number of unique categories in the first variable (number of rows)
#' \item \eqn{c} the number of unique categories in the second variable (number of columns)
#' \item \eqn{P} double the number of concordant pairs
#' \item \eqn{Q} double the number of discordant pairs
#' }
#' 
#' The formula for \eqn{ASE_1} can be found in Goodman and Kruskal (1963, p. 324) and also in
#' Brown and Benedetti (1977, p. 310). 
#' The formula for \eqn{ASE_0} was found in Brown and Benedetti (1977, p. 311) who also refer
#' to Goodman and Kruskal (1972).
#' 
#' Note that Kendall \eqn{\tau_a} is the same as Goodman-Kruskall gamma.
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' GoodmanKruskalGamma(table(ord1, ord2), conf.level=0.95)
#' 
#' *library(MESS)*
#' 
#' gkgamma(table(ord1, ord2))
#' 
#' *library(ryouready)*
#' 
#' ord.gamma(table(ord1, ord2))
#' 
#' @references 
#' Brown, M. B., & Benedetti, J. K. (1977). Sampling behavior of test for correlation in two-way contingency tables. *Journal of the American Statistical Association, 72*(358), 309–315. https://doi.org/10.2307/2286793
#' 
#' Goodman, L. A., & Kruskal, W. H. (1954). Measures of association for cross classifications. J*ournal of the American Statistical Association, 49*(268), 732–764. https://doi.org/10.2307/2281536
#' 
#' Goodman, L. A., & Kruskal, W. H. (1963). Measures of Association for Cross Classifications III: Approximate Sampling Theory. *Journal of the American Statistical Association, 58*(302), 310–364. https://doi.org/10.1080/01621459.1963.10500850
#' 
#' Goodman, L. A., & Kruskal, W. H. (1972). Measures of Association for Cross Classifications IV: Simplification of Asymptotic Variances. *Journal of the American Statistical Association, 67*(338), 415–421. https://doi.org/10.1080/01621459.1972.10482401
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' ord1 = c(5, 5, 5, 5, 3, 3, 4, 3, 4, 3, 4, 4, 4, 5, 3, 1, 3, 2, 1, 1, 1, 1, 2, 1, 2, 1, 2, 5, 3, 2, 1, 3, 4, 1, 1, 2, 3, 1, 1, 1, 3, 1, 1, 3, 4, 2, 2, 1, 1, 1, 2, 1, 1, 3)
#' ord2 = c(5, 5, 4, 5, 3, 3, 3, 3, 3, 5, 4, 3, 4, 5, 3, 2, 5, 2, 1, 1, 2, 2, 2, 3, 2, 3, 1, 5, 3, 3, 1, 3, 4, 1, 1, 2, 4, 1, 1, 3, 3, 1, 1, 3, 4, 2, 3, 2, 1, 2, 4, 2, 2, 4)
#' r_goodman_kruskal_gamma(ord1, ord2)
#' 
#' @export
r_goodman_kruskal_gamma <- function(ord1, ord2){
  
  ct = table(ord1, ord2)
  
  nr = nrow(ct)
  nc = ncol(ct)
  
  topLeft = matrix(0, nr, nc)
  for (iCell in 2:nr) {
    for (jCell in 2:nc){
      obs = 0
      for (i in 1:(iCell - 1)){
        for (j in 1:(jCell - 1)){
          obs = obs + ct[i,j]
        }
      }
      topLeft[iCell, jCell] = obs
    }
  }
  
  botRight = matrix(0, nr, nc)
  for (iCell in 1:(nr-1)) {
    for (jCell in 1:(nc-1)){
      obs = 0
      for (i in (iCell+1):nr){
        for (j in (jCell+1):nc){
          obs = obs + ct[i,j]
        }
      }
      botRight[iCell, jCell] = obs
    }
  }
  C = topLeft + botRight
  
  topRight = matrix(0, nr, nc)
  for (iCell in 2:nr) {
    for (jCell in 1:(nc-1)){
      obs = 0
      for (i in 1:(iCell - 1)){
        for (j in (jCell + 1):nc){
          obs = obs + ct[i,j]
        }
      }
      topRight[iCell, jCell] = obs
    }
  }
  
  botLeft = matrix(0, nr, nc)
  for (iCell in 1:(nr-1)) {
    for (jCell in 2:nc){
      obs = 0
      for (i in (iCell+1):nr){
        for (j in 1:(jCell - 1)){
          obs = obs + ct[i,j]
        }
      }
      botLeft[iCell, jCell] = obs
    }
  }
  D = topRight + botLeft
  
  P = sum(ct*C)
  Q = sum(ct*D)
  
  g = (P - Q)/(P + Q)
  
  n = sum(ct)
  ASE1 = 4*sqrt(sum(ct*(Q*C - P*D)**2))/(P + Q)**2
  ASE0 = 2*sqrt(sum(ct*(C - D)**2) - (P - Q)**2/n)/(P + Q)
  z = g/ASE0  

  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic = z
  results = data.frame(g, ASE1, ASE0, statistic, pValue)
  
  return(results)
  
}

