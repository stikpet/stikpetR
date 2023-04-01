#' Spearman Exact Distribution using Permutations
#' @description 
#' This code in this function was posted by cuttlefish44 (2016)
#' It performs a two-tailed exact test for Spearman rho.
#' @param ord1 the numeric scores of the first variable
#' @param ord2 the numeric scores of the second variable
#' @returns 
#' \item{pValue}{the two-tailed p-value}
#' 
#' @details 
#' 
#' The exact distribution is calculated using the following steps:
#' \enumerate{
#' \item Determine all possible permutations of the scores in the first variable
#' \item Determine for each permutation the Spearman rho with the second variable
#' \item Count how often the Spearman rho is above the Spearman rho between the original two variables
#' \item Divide the results by \eqn{n!}
#' }
#' 
#' @examples 
#' ord1 = c(5, 3, 3, 4, 3, 4, 3)
#' ord2 = c(5, 3, 3, 3, 3, 3, 5)
#' 
#' he_spearman_permutation(ord1, ord2)
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
#' 
#' @export
he_spearman_permutation <- function(ord1, ord2){
  n = length(ord1)
  rx = rank(ord1)
  ry = rank(ord2)
  rho = cor(rx, ry) 
  
  permutation <- he_permutations(n)
  # the function to calculate rho between an argument and ord2
  f.rho <- function(a) 1 - 6 * sum((rank(a) - rank(ord2))^2) / (n^3 - n)
  
  x2.all.permutation <- matrix(ord1[permutation], ncol=n)
  x2.all.permutation.rho <- apply(x2.all.permutation, 1, function(a) f.rho(a))
  
  pValue = sum(x2.all.permutation.rho > rho) / factorial(n) * 2
  
  return(pValue)
}