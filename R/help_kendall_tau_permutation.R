#' Helper Function - Kendall Tau Permutation Test 
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @returns 
#' \item{pValue}{upper tail p-value of Kendall tau Distribution}
#' @details 
#' Uses a permutation test to calculate the probability. It is an adaption to 
#' the code that was posted by cuttlefish44 (2016) online.
#' 
#' The exact distribution is calculated using the following steps:
#' \enumerate{
#' \item Determine all possible permutations of the scores in the first variable
#' \item Determine for each permutation the Kendall tau with the second variable
#' \item Count how often the Spearman rho is above the Kendall tau between the original two variables
#' \item Divide the results by \eqn{n!}
#' }
#' 
#' @examples 
#' ord1 = c(5, 8, 6, 3, 2, 9)
#' ord2 = c(2, 1, 4, 5, 8, 7)
#' he_tau_permutation(ord1, ord2)
#' 
#' @references 
#' cuttlefish44. (2016, September 16). Answer to “Different methods for finding spearman’s coefficient produce diff p-values depending on presence of tied values.” Cross Validated. https://stats.stackexchange.com/a/235380/190640
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @export
he_tau_permutation <- function(ord1, ord2){
  n = length(ord1)
  rx = rank(ord1)
  ry = rank(ord2)
  tau = cor(rx, ry, method="kendall") 
  print(tau)
  permutation <- he_permutations(n)
  # the function to calculate tau between an argument and ord2
  f.tau <- function(a) cor(rank(a), ry, method="kendall") 
  
  x2.all.permutation <- matrix(ord1[permutation], ncol=n)
  x2.all.permutation.tau <- apply(x2.all.permutation, 1, function(a) f.tau(a))
  
  pValue = sum(x2.all.permutation.tau > tau) / factorial(n)
  
  return(pValue)
}