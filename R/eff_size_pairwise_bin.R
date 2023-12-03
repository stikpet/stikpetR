#' Binary Effect Size for Pairwise Test
#' @description 
#' When using a pairwise post-hoc test for a single nominal variable, the pair has become binary. This function then can determine the effect size for each pair.
#' 
#' Options are to use Cohen g, Cohen h', or the Alternative Ratio.
#' 
#' @param data, list with the data
#' @param expCounts dataframe, optional. The categories and expected counts
#' @param es, string, optional. effect size to use.Either "coheng" (default), "cohenh", "ar"
#' 
#' @returns
#' dataframe with 
#' \item{cat1}{label of first category in pair}
#' \item{cat2}{label of second category in pair}
#' \item{n1}{number of cases in first category}
#' \item{n2}{number of cases in second category}
#' \item{var}{followed by the effect size value}
#' 
#' @details
#' If expected counts are provided, for Cohen h' and the Alternative Ratio these expected counts are converted to expected proportions
#' 
#' See the separate functions of each effect size for more details.
#' 
#' * \code{\link{es_cohen_g}} for Cohen g
#' * \code{\link{es_cohen_h_os}} for Cohen h'
#' * \code{\link{es_alt_ratio}} for the Alternative Ratio
#' 
#' @seealso 
#' \code{\link{ph_binomial}}, performs a pairwise binomial test
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_pairwise_bin <- function(data, expCounts=NULL, es="coheng"){
  
  Fi = table(data)
  categ = names(Fi)
  res = as.data.frame(t(combn(categ, 2)))
  nPairs = nrow(results)
  res["n1"] = as.integer(Fi[res[,1]])
  res["n2"] = as.integer(Fi[res[,2]])
  
  if (es=="coheng"){
    res["Cohen g"] = res['n1'] / (res['n1']+res['n2']) - 0.5
  }
  
  else if (es=="cohenh"){
    pi = res['n1'] / (res['n1']+res['n2'])
    phi1 = 2*asin(pi**0.5)
    
    if (is.null(expCounts)){
      phic = 2*asin(0.5**0.5)
      res["Cohen h\'"] = phi1 - phic            
    }
    else {
      n1E = c()
      n2E = c()
      for (i in 1:nPairs){
        n1E = c(n1E, expCounts[expCounts[,1]==res[i,1], 2])
        n2E = c(n2E, expCounts[expCounts[,1]==res[i,2], 2])
      }
      pic = n1E / (n1E + n2E)
      phic = 2*asin(pic**0.5)
      h = phi1 - phic
      res["Cohen h\'"] = h
    }
  }
  
  else if (es=="ar"){
    pi1 = res['n1'] / (res['n1']+res['n2'])
    pi2 = res['n2'] / (res['n1']+res['n2'])
    if (is.null(expCounts)){
      ar1 = pi1/0.5
      ar2 = pi2/0.5
    }
    else {
      n1E = c()
      n2E = c()
      for (i in 1:nPairs){
        n1E = c(n1E, expCounts[expCounts[,1]==res[i,1], 2])
        n2E = c(n2E, expCounts[expCounts[,1]==res[i,2], 2])
      }
      pic1 = n1E / (n1E + n2E)
      pic2 = n2E / (n1E + n2E)
      ar1 = pi1/pic1
      ar2 = pi2/pic2
      
    }
    res['AR 1'] = ar1
    res['AR 2'] = ar2
    
  }
  
  
  return(res)
}