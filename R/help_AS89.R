#' Helper Function - Algorithm AS 89
#' 
#' @param n the number of scores (should be equal in both variables)
#' @param IS the test statistic (see details)
#' @returns 
#' \item{pValue}{the two-sided significance (p-value)}
#' 
#' @description 
#' Algorithm AS 89 (Best & Roberts, 1975) is for upper-tail probabilities
#' 
#' 
#' @details 
#' The test statistic \eqn{S} defined as:
#' \deqn{S = \sum_{i = 1}^n d_i^2 = \sum_{i=1}^n \left(r_{x_i} - r_{y_i}\right)^2}
#' Which if there are no ties is equal to:
#' \deqn{S = \frac{\left(n^3 - n\right)\times\left(1 - r_s\right)}{6}}
#' 
#' @references 
#' Best, D. J., & Roberts, D. E. (1975). Algorithm AS 89: The upper tail probabilities of Spearman's rho. *Applied Statistics, 24*(3), 377-379. https://doi.org/10.2307/2347111
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' n = 12
#' S = 8
#' he_AS89(n, S)
#' 
#' @export
he_AS89 <- function(n, IS){
  
  L = rep(0, 6)
  
  #Test admissibility of arguments and initialize
  PRHO = 1
  if (n <= 1 || IS <= 0) {
    return (PRHO)
  }
  
  PRHO = 0
  if (IS > n*(n*n - 1)/3) {
    return (PRHO)
  }
  
  JS = IS
  if (JS != 2*(JS/2)) {
    JS = JS + 1
  }
  
  if (n > 6) {
    #GOTO 6
    #6
    B = 1/n
    X = (6*(JS - 1)*B/(1/(B*B)-1)-1)*sqrt(1/B - 1)
    Y = X*X
    U = X*B*(0.2274 + B*(0.2531 + 0.1745*B) + Y*(-0.0758 + B*(0.1033 + 0.3932*B) - Y*B*(0.0879 + 0.0151*B - Y*(0.0072 - 0.0831*B + Y*B*(0.0131 - 0.00046*Y)))))
    PRHO = U/exp(Y/2) + pnorm(X, lower.tail = FALSE)
    if (PRHO < 0) {
      PRHO = 0
    }
    else if (PRHO>1){
      PRHO = 0
    }
    
    return(PRHO)
  }
  
  #Exact evaluation of probability
  
  NFAC = 1
  
  for (I in 1:n) {
    NFAC = NFAC*I
    L[I] = I
  }
  #1 continue
  
  PRHO = 1/NFAC
  if (JS == (n*(n*n - 1)/3)) {
    return(PRHO)
  }
  IFR = 0
  
  for (M in 1:NFAC) {
    ISE = 0
    for (I in 1:n) {
      ISE = ISE + (I - L[I])**2
    }
    #2 continue
    
    if (JS <= ISE) {
      IFR = IFR + 1
    }
    N1 = n
    
    #3
    LOOP3 = TRUE
    while (LOOP3) {
      LOOP3 = FALSE
      MT = L[1]
      NN = N1 - 1
      for (I in 1:NN) {
        L[I] = L[I + 1]
      }
      #4 continue
      
      L[N1] = MT
      if (L[N1] == N1 && N1 == 2) {
        N1 = N1 - 1
        
        if (M != NFAC) {
          LOOP3 = TRUE
        }
      }
    }
  }
  PRHO = IFR/NFAC
  return(PRHO)
  
}



