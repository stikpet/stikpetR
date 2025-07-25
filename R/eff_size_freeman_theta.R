#' Freeman Theta
#' @description 
#' According to Jacobson (1972, p. 42), this is the only measure for nominal-ordinal data, and is a modification of Somers d.
#' 
#' It can range from 0 to 1, with 0 indicating no influence of the catField on the scores of the ordField, and a 1 a perfect relationship.
#' 
#' Alternatives could be eta-squared and epsilon-squared.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' 
#' @returns
#' theta, float. The Freeman Theta value
#' 
#' @details 
#' The formula used is (Freeman, 1965, p. 116):
#' \deqn{\theta = \frac{D}{T}}
#' 
#' With:
#' \deqn{D = \sum D_{x,y}}
#' \deqn{D_{x,y} = \left|f_a - f_b\right|}
#' \deqn{f_a = \sum_{i=1}^{n_{lvl} - 1}\left(F_{x,i}\times\sum_{j=i+1}^{n_{lvl}} F_{y,j}\right)}
#' \deqn{f_b = \sum_{i=2}^{n_{lvl}}\left(F_{x,i}\times\sum_{j=1}^{i-1} F_{y,j}\right)}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{F_{x,i}}, from category x, the number of cases with level i.
#' \item \eqn{n_{lvl}}, the number of levels.
#' \item \eqn{n_i}, the total number of cases from category i
#' }
#' 
#' @references
#' Freeman, L. C. (1965). *Elementary applied statistics: For students in behavioral science*. Wiley.
#' 
#' Jacobson, P. E. (1972). Applying measures of association to nominal-ordinal data. *The Pacific Sociological Review, 15*(1), 41-60. doi:10.2307/1388286
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_freeman_theta <- function(catField, ordField, categories=NULL, levels=NULL){
  #create the cross table    
  ct = tab_cross(catField, ordField, order1=categories, order2=levels, totals="include")
  
  #basic counts
  k = nrow(ct)-1
  nlvl = ncol(ct)-1
  
  d = 0
  t = 0
  for (x in 1:(k - 1)){
    for (y in (x + 1):k){
      fb = 0
      for (i in 2:nlvl){
        fs = 0
        for (j in 1:(i-1)){
          fs = fs + ct[y, j]
        }
        fb = fb + ct[x, i] * fs
      }
      
      fa = 0
      for (i in 1:(nlvl - 1)){
        fs = 0
        for (j in (i + 1):nlvl){
          fs = fs + ct[y, j]
        }
        fa = fa + ct[x, i] * fs
      }
      
      d = d + abs(fa - fb)
      t = t + ct[x, nlvl+1] * ct[y, nlvl+1]
    }
  }
  
  theta = d / t
  
  return (theta)
}



