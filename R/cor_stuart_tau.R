#' Stuart Tau c / Kendall Tau c
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @param cc boolean to indicate the use of a continuity correction
#' @returns 
#' A dataframe with:
#' \item{tau}{the correlation coefficient}
#' \item{statistic}{the statistic from the test}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details
#' 
#' Tau looks at concordant pairs versus discordant pairs. These can be calculated using:
#' \deqn{n_c = \frac{P}{2}}
#' \deqn{n_d = \frac{Q}{2}}
#' If the scores are placed in a cross table we can use:
#' \deqn{P = \sum_{i=1}^r \sum_{j=1}^c P_{i,j}}
#' \deqn{Q = \sum_{i=1}^r \sum_{j=1}^c Q_{i,j}}
#' \deqn{P_{i,j} = F_{i,j}\times C_{i,j}}
#' \deqn{Q_{i,j} = F_{i,j}\times D_{i,j}}
#' \deqn{C_{i,j} = \sum_{h<i}\sum_{k<j} F_{h,k} + \sum_{h>i}\sum_{k>j} F_{h,k}}
#' \deqn{D_{i,j} = \sum_{h<i}\sum_{k>j} F_{h,k} + \sum_{h>i}\sum_{k<j} F_{h,k}}
#' 
#' Alternative, we don't have to use a cross table:
#' \deqn{P = \sum_{i=1}^n P_i}
#' \deqn{Q = \sum_{i=1}^n Q_i}
#' \deqn{P_i = \sum_{j=1}^n \begin{cases} 1 & \text{sign} \left(x_i - x_j\right) \times \text{sign} \left(y_i - y_j\right) = 1 \\ 0 & \text{ else}\end{cases}}
#' \deqn{Q_i = \sum_{j=1}^n \begin{cases} 1 & \text{sign} \left(x_i - x_j\right) \times \text{sign} \left(y_i - y_j\right) = -1 \\ 0 & \text{ else}\end{cases}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_c} the number of concordant pairs
#' \item \eqn{n_d} the number of discondant pairs
#' \item \eqn{n} is the number of pairs
#' \item \eqn{r} the number of categories in the first variable (i.e. number of rows)
#' \item \eqn{c} the number of categories in the second variable (i.e. number of columns)
#' \item \eqn{F_{h,k}} is the number of cases that scored h for the first variable, and k for the second
#' \item \eqn{P} is double the number of concordant pairs
#' \item \eqn{Q} is double the number of discordant pairs
#' }
#' 
#' The formula used is (Stuart, 1953, p. 107):
#' \deqn{\tau_c = \frac{P - Q}{n^2\times\frac{m - 1}{m}}}
#' With:
#' \deqn{m = \text{min}\left(r, c\right)}
#' 
#' **Testing**
#' 
#' The following normal approximation can be used (Brown & Benedetti, 1977, p. 311):
#' \deqn{z_{c} = \frac{\tau_c}{ASE}}
#' \deqn{ASE_0 = \frac{2\times m}{\left(m - 1\right)^2}\times\sqrt{\sum_{i=1}^r \sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2 - \frac{\left(P - Q\right)^2}{n}}}
#' 
#' The calculation of \eqn{\sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2}
#' can then also be accomplished using:
#' \deqn{\sum_{i=1}^n \left(P_i - Q_i\right)^2}
#' 
#' The continuity correction is applied as (Schaeffer & Levitt, p. 342):
#' \deqn{\tau_{cc} = \left|\tau\right| - \frac{2}{n\times\left(n - 1\right)}}
#' 
#' Note that this correction should actually be adjusted in case ties are present.
#' Hopefully this can be implemented in a future update.
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' StuartTauC(ord1, ord2)
#' 
#' *library(ryouready)*
#' 
#' ord.tau(table(ord1, ord2))
#' 
#' @examples 
#' ord1 = c(78,78,60,53,85,84,73,78,78,75,65,72,58,92,65)
#' ord2 = c(4,23,25,48,17,8,4,26,23,19,24,35,29,4,14)
#' 
#' r_stuart_tau(ord1, ord2)
#' r_stuart_tau(ord1, ord2, cc=TRUE)
#' 
#' @references 
#' Brown, M. B., & Benedetti, J. K. (1977). Sampling behavior of test for correlation in two-way contingency tables. *Journal of the American Statistical Association, 72*(358), 309–315. https://doi.org/10.2307/2286793
#' 
#' Schaeffer, M. S., & Levitt, E. E. (1956). Concerning Kendall’s tau, a nonparametric correlation coefficient. *Psychological Bulletin, 53*(4), 338–346. https://doi.org/10.1037/h0045013
#' 
#' Stuart, A. (1953). The estimation and comparison of strengths of association in contingency tables. *Biometrika, 40*(1/2), 105. https://doi.org/10.2307/2333101
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' 
#' @export
r_stuart_tau <- function(ord1, ord2, cc=FALSE){
  
  datFrame = na.omit(data.frame(ord1, ord2))
  ct = table(datFrame$ord1, datFrame$ord2)
  
  nr = nrow(ct)
  nc = ncol(ct)
  
  datFrame = na.omit(data.frame(ord1, ord2))
  n = length(datFrame$ord1)
  ASE0 = P = Q = 0
  for (i in 1:n) {
    pC = pD = 0
    for (j in 1:n){
      if (ord1[i] != datFrame$ord1[j] && datFrame$ord2[i] != datFrame$ord2[j]) {
        if(sign(datFrame$ord1[i] - datFrame$ord1[j]) == sign(datFrame$ord2[i] - datFrame$ord2[j])){
          P = P + 1
          pC = pC + 1
        }
        else {
          pD = pD + 1
          Q = Q + 1
        }
      }
    }    
    ASE0 = ASE0 + (pC - pD)**2
  }
  
  nco = P/2
  nd = Q/2
  
  m = min(nr, nc)
  tau = (P - Q)/(n**2*(m - 1)/m)
  
  tauTest = tau
  if(cc){
    tauTest = abs(tau) - 2/(n*(n - 1))
  }
  
  ASE0 = 4*m**2/((m - 1)**2*n**4)*(ASE0 - (P - Q)**2/n)
  
  z = tauTest/sqrt(ASE0)
  
  
  pValue = 2*(1 - pnorm(abs(z)))
  
  results = data.frame(tau, ASE0, z, pValue)
  
  return(results)
}