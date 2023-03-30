#' Helper Function for Permutations
#' @description 
#' This function was posted by Museful (2013). It creates all possible permutations.
#' 
#' @param n the number of scores
#' @returns 
#' \item{A}{all possible permutations of integers from 1 to n}
#' 
#' @references 
#' Museful. (2013, November 25). Answer to “Generating all distinct permutations of a list in R.” Stack Overflow. https://stackoverflow.com/a/20199902/12149706
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' he_permutations(5)
#' 
#' @export
he_permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}
