#' Kendall Tau (a and b)
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @param ver for tau a, or tau b
#' @param test which test to use (see details)
#' @param cc boolean to indicate the use of a continuity correction
#' @returns 
#' A dataframe with:
#' \item{tau}{the correlation coefficient}
#' \item{statistic}{the statistic from the test (only if applicable)}
#' \item{df}{the degrees of freedom (only if applicable)}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details
#' 
#' Kendall tau looks at concordant pairs versus discordant pairs. These can be calculated using:
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
#' #' *Symbols used:*
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
#' 
#' **Tau a**
#' 
#' The formula used for the value of \eqn{\tau_a} is (Kendall, 1938, p. 82):
#' \deqn{\tau_a = \frac{P - Q}{n\times\left(n - 1\right)}}
#' 
#' The formula can also be written as:
#' \deqn{\tau_a = \frac{n_c - n_d}{n_0}}
#' With:
#' \deqn{n_0 = \frac{n\times\left(n - 1\right)}{2}}
#' 
#' Note that Kendall \eqn{\tau_a} is the same as Goodman-Kruskall gamma.
#' 
#' **Tau b**
#' 
#' For \eqn{\tau_b} the formula used is (Kendall, 1945, p. 243):
#' \deqn{\tau_b = \frac{P - Q}{\sqrt{D_r \times D_c}}}
#' 
#' With:
#' \deqn{D_r = n^2 - \sum_{i=1}^r RS_i^2}
#' \deqn{D_c = n^2 - \sum_{j=1}^c CS_i^2}
#' \deqn{RS_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{CS_i = \sum_{i=1}^r F_{i,j}}
#' 
#' Note that \eqn{RS_i} are the frequencies of the scores in the first variable and
#' \eqn{CS_i} are the frequencies of the scores in the second variable.
#' 
#' Alternative the formula can be written as:
#' \deqn{\tau_b = \frac{n_c - n_d}{\sqrt{\left(n_0 - t_1\right)\times\left(n_0 - t_2\right)}}}
#' With:
#' \deqn{t_1 = \sum_{i=1}^n \frac{RS_i*\left(RS_i - 1\right)}{2}}
#' \deqn{t_2 = \sum_{i=1}^n \frac{CS_i*\left(CS_i - 1\right)}{2}}
#' 
#' **Testing**
#' 
#' For *Tau a* the following normal approximation can be used (Kendall, 1962, p. 51):
#' \deqn{z_{a} = \frac{3\times\left(\frac{P - Q}{2}\right)}{\sqrt{\frac{n\times\left(n - 1\right)\times\left(2\times n + 5\right)}{2}}}} 
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{a}\right|\right)\right)}
#' Or written with \eqn{\tau_a} (Schaeffer & Levitt, p. 341):
#' \deqn{z_{a} = \frac{\tau_a}{\sqrt{\frac{\sigma_{\tau_a}^2}{n}}}}
#' With:
#' \deqn{\sigma_{\tau_a}^2 = \frac{4\times n + 10}{9\times n\times\left(n - 1\right)}}
#' 
#' For *Tau b* an approximation can be used with:
#' \deqn{z_{b} = \frac{\tau_b}{ASE}}
#' 
#' For the equation of ASE two variations to choose from:
#' Using Brown and Benedetti (1977, p. 311):
#' \deqn{ASE_0 = 2\times\sqrt{\frac{\sum_{i=1}^r \sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2 - \frac{\left(P - Q\right)^2}{n}}{D_r\times D_c}}}
#' This is used when *test="bb"*
#' 
#' The calculation of \eqn{\sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2}
#' can then also be accomplished using:
#' \deqn{\sum_{i=1}^n \left(P_i - Q_i\right)^2}
#' 
#' Or a version from Kendall (1962, p. 55):
#' \deqn{ASE = 2\times\sqrt{v}}
#' With:
#' \deqn{v = \frac{v_0 - v_r - v_c}{18} + v_1 + v_2}
#' \deqn{v_0 = n\times\left(n - 1\right)\times\left(2\times n + 5\right)}
#' \deqn{v_r = \sum_{i=1}^r RS_i\times\left(RS_i - 1\right)\times\left(2\times RS_i + 5\right)}
#' \deqn{v_c = \sum_{j=1}^c CS_j\times\left(CS_j - 1\right)\times\left(2\times CS_j + 5\right)}
#' \deqn{v_1 = \frac{\left(\sum_{i=1}^r RS_i\times\left(RS_i - 1\right)\times\left(RS_i - 2\right)\right) \times \left(\sum_{j=1}^c CS_j\times\left(CS_j - 1\right)\times\left(CS_j - 2\right)\right)}{9\times n\times\left(n - 1\right)\times\left(n - 2\right)}}
#' This is used when *test="kendall"* and the default
#' 
#' For the AS 71 algorithm the test statistic is given by:
#' \deqn{S = \binom{n}{2}\times\left|\tau\right| = \frac{n\times\left(n - 1\right)}{2}\times\left|\tau\right|}
#' 
#' The absolute value of \eqn{\tau} is used, since AS 71 only provides the upper-tail probabilities.
#' 
#' For the Kendall method the test statistic is the number of concordant pairs \eqn{n_c}.
#' 
#' See the documentation of *di_kendall_tau()* for more info on the AS 71 and Kendall algorithms.
#' 
#' The continuity correction is applied as (Schaeffer & Levitt, p. 342):
#' \deqn{\tau_{cc} = \left|\tau\right| - \frac{2}{n\times\left(n - 1\right)}}
#' or depending on the approximation used (Kendall, 1961, p. 54):
#' \deqn{S_{cc} = \left|S\right| - 1}
#' Where:
#' \deqn{S = n_c - n_d}
#' 
#' Note that this correction should actually be adjusted in case ties are present.
#' Hopefully this can be implemented in a future update.
#' 
#' **Alternatives**
#' 
#' *R's stats library*
#' 
#' cor.test(ord1, ord2, method="kendall", exact=FALSE)
#' 
#' cor.test(ord1, ord2, method="kendall", exact=TRUE)
#' 
#' cor.test(ord1, ord2, method="kendall", exact=FALSE, continuity = TRUE)
#' 
#' *library(ryouready)*
#' 
#' ord.tau(table(ord1, ord2))
#' 
#' @examples 
#' ord1 = c(5, 8, 6, 3, 2, 9)
#' ord2 = c(2, 1, 4, 5, 7, 8)
#' r_kendall_tau(ord1, ord2, ver="a")
#' r_kendall_tau(ord1, ord2, ver="a", cc=TRUE)
#' r_kendall_tau(ord1, ord2, ver="b", test="kendall")
#' r_kendall_tau(ord1, ord2, ver="b", test="bb")
#' r_kendall_tau(ord1, ord2, ver="b", test="bb", cc=TRUE)
#' r_kendall_tau(ord1, ord2, ver="b", test="kendall", cc=TRUE)
#' r_kendall_tau(ord1, ord2, ver="b", test="exact-AS71")
#' r_kendall_tau(ord1, ord2, ver="b", test="exact-kendall")
#' 
#' ord1 = c(5, 8, 6, 3, 2, 9)
#' ord2 = c(8, 7, 4, 1, 3, 6)
#' r_kendall_tau(ord1, ord2, ver="b")
#' 
#' ord1 = c(5, 8, 6, 3, 2, 9, 11, 13, 1)
#' ord2 = c(8, 7, 4, 1, 3, 6, 5, 11, 2)
#' r_kendall_tau(ord1, ord2, ver="b")
#' 
#' ord1 = c(5, 3, 3, 4, 3, 4, 3, 4, 4, 4, 5, 3, 1, 3, 2)
#' ord2 = c(5, 3, 3, 3, 3, 3, 5, 4, 3, 4, 5, 3, 2, 5, 2)
#' r_kendall_tau(ord1, ord2, ver="b")
#' 
#' @references 
#' Brown, M. B., & Benedetti, J. K. (1977). Sampling behavior of test for correlation in two-way contingency tables. *Journal of the American Statistical Association, 72*(358), 309–315. https://doi.org/10.2307/2286793
#' 
#' Kendall, M. G. (1938). A new measure of rank correlation. *Biometrika, 30*(1–2), 81–93. https://doi.org/10.1093/biomet/30.1-2.81
#' 
#' Kendall, M. G. (1945). The treatment of ties in ranking problems. *Biometrika, 33*(3), 239–251. https://doi.org/10.1093/biomet/33.3.239
#' 
#' Kendall, M. G. (1962). *Rank correlation methods* (3rd ed.). Charles Griffin.
#' 
#' Schaeffer, M. S., & Levitt, E. E. (1956). Concerning Kendall’s tau, a nonparametric correlation coefficient. *Psychological Bulletin, 53*(4), 338–346. https://doi.org/10.1037/h0045013
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
r_kendall_tau <- function(ord1, ord2, ver=c("a", "b"), test=c("kendall", "bb", "exact-AS71", "exact-kendall"), cc=FALSE){
  
  if (length(ver)>1) {
    ver="b"
  }
  
  if (length(test)>1) {
    test="kendall"
  }
  
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
  
  if (ver=="a") {
    tau = (nco - nd)/(n * (n - 1)/2)
    
    tauTest = tau
    if(cc){
      tauTest = abs(tau) - 2/(n*(n - 1))
    }
    
    z = 3*tauTest/sqrt(((4*n + 10)/(n - 1))/n)
    pValue = 2*(1 - pnorm(abs(z)))
    
    statistic = z
    if (tau < 0) {
      statistic = -abs(z)
    }

  }
  else if (ver=="b"){
    
    t1 = table(datFrame$ord1)
    t2 = table(datFrame$ord2)
    
    #check for ties:
    if (max(t1)>1 || max(t2)>1) {
      if (test=="AS71" || test=="kendall") {
        test="kendall"
        warning("Ties present, switch method to Kendall-Gibbons") 
      }
    }

    Dr = n**2 - sum(t1**2)
    Dc = n**2 - sum(t2**2)
    
    tau = (P - Q)/sqrt(Dr*Dc)
    tauTest = tau
    if (test=="bb") {
      ASE0 = 2*sqrt((ASE0 - (P - Q)**2/n)/(Dr*Dc))
      
      if(cc){
        tauTest = abs(tau) - 2/(n*(n - 1))
      }
      
      z = tauTest/(ASE0)
      pValue = 2*(1 - pnorm(abs(z)))
      statistic = z
      if (tau < 0) {
        statistic = -abs(z)
      }
    }
    else if (test=="kendall") {
      #Kendall (1962, p. 55)
      v0 = n*(n - 1)*(2*n + 5)
      vt1 = sum(t1*(t1 - 1)*(2*t1 + 5))
      vt2 = sum(t2*(t2 - 1)*(2*t2 + 5))
      v1 = sum(t1*(t1 - 1))*sum(t2*(t2 - 1))/(2*n*(n - 1))
      v2 = sum(t1*(t1 - 1)*(t1 - 2))*sum(t2*(t2-1)*(t2-2))/(9*n*(n-1)*(n-2))
      v = (v0 - vt1 - vt2)/18 + v1 + v2
      
      if (cc) {
        z = (abs(nco - nd)-1)/sqrt(v)
      }
      else {
        z = (nco - nd)/sqrt(v)
      }
      pValue = 2*(1 - pnorm(abs(z)))
      statistic = z
      if (tau < 0) {
        statistic = -abs(z)
      }
    }
    
    else if (test=="exact-AS71"){
      S = choose(n, 2)*abs(tau)
      pValue = di_kendall_tau(n, tau, method="AS71")
      statistic = S
    }
    else if (test=="exact-kendall"){
      pValue = di_kendall_tau(n, tau, method="kendall")
      statistic = nco
    }
  }
  
  results = data.frame(tau, statistic, pValue)
  
  return(results)
}
