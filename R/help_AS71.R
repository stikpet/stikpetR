#' Helper Function - Algorithm AS 71
#' @param S the test statistic
#' @param N the sample size
#' @returns 
#' \item{pValue}{upper tail p-value of Kendall tau Distribution}
#' @details 
#' Algorithm AS 71 (Best & Gipps, 1974) uses as test statistic:
#' \deqn{S = \binom{n}{2}\times\left|\tau\right| = \frac{n\times\left(n - 1\right)}{2}\times\left|\tau\right|}
#' 
#' The Fortran code was translated to R by myself.
#' 
#' @references 
#' Best, D. J., & Gipps, P. G. (1974). Algorithm AS 71: The upper tail probabilities of Kendall’s tau. *Applied Statistics, 23*(1), 98–100. https://doi.org/10.2307/2347062
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
he_AS71 <- function(S, N){
  #Upper tail probabilities of tau
  H = rep(0, 15)
  L = matrix(0, nrow=2, ncol=15)
  
  #Check validity of IS and N values
  PRTAUS = 1
  if (N < 1) {
    return (PRTAUS)
  }
  M = N*(N - 1) / 2 - abs(S)
  if ((M < 0) || (M > (M/2)*2)) {
    return (PRTAUS)
  }
  if (M == 0 && S <= 0) {
    return (PRTAUS)
  }
  if (N > 8) {
    #GOTO 7
    
    #7
    #Calculation of Tchebycheff-Hermite polynomials
    X = (S - 1)/sqrt((6 + N*(5-N*(3 + 2*N)))/(-18))
    H[1] = X
    H[2] = X*X - 1
    
    for (I in 3:15) {
      #8
      H[I] = X*H[I-1] - (I - 1)*H[I-2]
    }
    
    #Probabilities calculated by modified Edgeworth series for n greater than 8
    
    R = 1/N
    SC = R*(H[3]*(-9E-2 + R*(4.5E-2 + R*(-5.325E-1 + R*5.06E-1))) + 
              R*(H[5]*(3.6735E-2 + R*(-3.6735E-2 + R*3.214E-1))
                 + H[7]*(4.05E-3 + R*(-2.3336E-2 + R*7.787E-2))
                 + R*(H[9]*(-3.3061E-3 - R*6.5166E-3)
                      + H[11]*(-1.215E-4 + R*2.5927E-3)
                      +R*(H[13]*1.4878E-4
                          + H[15]*2.7338E-6))))
    
    PRTAUS = pnorm(X, lower.tail=FALSE) + SC*0.398942*exp(-0.5*X*X)
    
    if (PRTAUS < 0) {
      PRTAUS = 0
    }
    if (PRTAUS > 1) {
      PRTAUS = 1
    }
    
    return(PRTAUS)
    
  }
  
  #Probabilities calculated by recurrence relation for N less than 9
  
  if (S < 0) {
    M = M - 2
  }
  
  IM = M /2 + 1
  L[1,1] = 1
  L[2,1] = 1
  
  if (IM >= 2) {
    for (I in 2:IM) {
      L[1, I] = 0
      L[2, I] = 0
    }
  }
  
  #2
  IL = 1
  I = 1
  M = 1
  J = 1
  JJ = 2
  
  #3
  GOTO3 = TRUE
  while (GOTO3) {
    if (I==N){
      GOTO3 = FALSE
    }
    else {
      GOTO3 = FALSE
      IL = IL + I
      I = I + 1
      M = M*I
      J = 3 - J
      JJ = 3 - JJ
      IN = 1
      IO = 0
      K = min(IM, IL)
      
      #4
      GOTO4=TRUE
      while(GOTO4){
        GOTO4=FALSE
        IN = IN + 1
        
        if (IN > K) {
          GOTO3 = TRUE
        }
        else {
          L[JJ, IN] = L[JJ, IN - 1] + L[J, IN]
          
          if (IN <= I) {
            GOTO4 = TRUE
          }
          else{
            IO = IO + 1
            L[JJ, IN] = L[JJ, IN] - L[J, IO]
            GOTO4 = TRUE
          }
        }
        
      }
    }
  }
  
  #5
  K = 0
  for (I in 1:IM) {
    #6
    K = K + L[JJ, I]
  }
  PRTAUS = K/M
  
  if (S < 0) {
    PRTAUS = 1 - PRTAUS
  }
  
  return(PRTAUS)
  
}