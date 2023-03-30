#' Helper Function - Kendall Algorithm 
#' @param n the sample size
#' @param c the number of concordant pairs
#' @returns 
#' \item{pValue}{upper tail p-value of Kendall tau Distribution}
#' @details 
#' An algorithm found at https://github.com/scipy/scipy/blob/v1.10.1/scipy/stats/_mstats_basic.py#L774-L898
#' was adapted. This refers to Kendall (1970), and uses the helper function *he_kendall(n, C)*.
#' Where \eqn{C = n_c}, i.e. the number of concordant pairs.
#' This algorithm already returns a two-tailed result.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
he_kendall <- function(n, c){
  
  c = as.integer(min(c, floor(n*(n - 1)/2) - c))
  
  if (n==1) {
    prob = 1
  }
  else if (n==2) {
    prob = 1
  }
  else if (c==0) {
    if (n < 171) {
      prob = 2/factorial(n)
    }
    else{
      prob=0
    }
  }
  else if (c==1) {
    if (n < 172) {
      prob = 2/factorial(n-1)
    }
    else{
      prob=0
    }
  }
  else if (4*c==n*(n-1)){
    prob = 1
    
  }
  else if (n < 171) {
    new = rep(0, c+1)
    new[1:2] = 1
    for (j in 3:n){
      new = cumsum(new)
      if (j <= c) {
        new[(j+1):(c+1)] = new[(j+1):(c+1)] - new[1:(c-j+1)]
      }
    }
    prob = 2*sum(new)/factorial(n)
  }
  
  else{
    new = rep(0, c+1)
    new[1:2] = 1
    for (j in 3:n){
      new = cumsum(new)/j
      if (j <= c) {
        new[(j+1):(c+1)] = new[(j+1):(c+1)] - new[1:(c-j+1)]
      }
    }
    prob = sum(new)
  }
  
  return(prob)
}