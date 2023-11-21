#' Berry-Johnston-Mielke R
#' @description 
#' A chance-corrected version of eta-squared, as an effect size measure for a Cochran Q test.
#' 
#' @param data dataframe with the scores
#' @param success indicator for what is considered a success (default is 1)
#' 
#' @returns 
#' \item{R}{the effect size measure}
#' 
#' @details 
#' The formula used (Berry et al., 2007 pp. 1237, 1239):
#' \deqn{R = 1 - \frac{\delta}{\mu_{\delta}}}
#' With:
#' \deqn{\mu_{\delta} = \frac{2}{n\times\left(n - 1\right)} \times \left(\sum_{i=1}^n p_i\right) \times \left(n - \sum_{i=1}^n p_i\right) - \sum_{i=1}^n p_i\times\left(1 - p_i\right)}
#' \deqn{\delta = \frac{1}{k\times\binom{n}{k}} \times \sum_{c=1}^k\sum_{i=1}^{n - 1}\sum_{j=i+1}^n \left|x_{i,c} - x_{j,c}\right|}
#' \deqn{p_i = \frac{\sum_{j=1}^k}{k}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{n} the number of rows
#' \item \eqn{k} the number of columns
#' \item \eqn{x_{i,j}} the score in row i and column j
#' }
#' 
#' The function actually uses for:
#' \deqn{\sum_{c=1}^k\sum_{i=1}^{n - 1}\sum_{j=i+1}^n \left|x_{i,c} - x_{j,c}\right| = \sum_{j=1}^k C_j \times\left(n - \sum_{j=1}^k C_j\right)}
#' With:
#' \deqn{C_j = \sum_{i = 1}^n x_{i,j}}
#' 
#' The original article has in the equation for \eqn{\mu_{\delta}} the first factor written as
#' \eqn{\frac{2}{k\times\left(k - 1\right)}}. In personal communication with one of the authors 
#' Alexis (2014) indicated this was wrong and \eqn{n} should be used.
#' 
#' @references 
#' Alexis. (2014, September 7). Answer to “Effect size of Cochran’s Q.” Cross Validated. https://stats.stackexchange.com/a/114649
#' 
#' Berry, K. J., Johnston, J. E., & Mielke, P. W. (2007). An alternative measure of effect size for Cochran’s Q test for related proportions. *Perceptual and Motor Skills, 104*(3_suppl), 1236–1242. doi:10.2466/pms.104.4.1236-1242
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_jbm_r <- function(data, success=NULL){
  dFr = na.omit(data)
  k = ncol(dFr)
  n = nrow(dFr)
  
  if (is.null(success)){
    success=dFr[1,1]}
  
  dFr <- dFr == success
  
  nColSuc = colSums(dFr)
  nFail = n - nColSuc
  nProd = nColSuc*nFail
  f1 = 1/(choose(n, 2)*k)
  delta = f1*sum(nProd)
  
  pRowSuc = rowSums(dFr)/k
  pComb = pRowSuc*(1 - pRowSuc)
  npSuc = sum(pRowSuc)
  npComb = sum(pComb)
  f2 = npSuc*(n - npSuc) - npComb
  mu = 2/(n*(n - 1))*f2
  
  R = 1 - delta/mu
  
  return(R)  
  
}