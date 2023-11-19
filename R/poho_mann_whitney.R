#' Post-Hoc Pairwise Mann-Whitney U Test
#' @description
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. A simple Bonferroni adjustment is also made for the multiple testing.
#' 
#' Other post-hoc tests that could be considered are Dunn, Nemenyi, Steel-Dwass, Conover-Iman, or pairwise Mood-Median.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param method string, optional. The use of the exact distribution (only if there are no ties) or the normal. Either "exact", "appr"
#' @param cc boolean, optional.Use of continuity correction. Default is True.
#' 
#' @returns
#' A dataframe with:
#' \item{category 1}{one of the two categories being compared}
#' \item{category 2}{second of the two categories being compared}
#' \item{statistic}{the test statistic (z-value)}
#' \item{p-value}{the p-value (significance)}
#' \item{adj. p-value}{the Bonferroni adjusted p-value}
#' 
#' @details
#' This function selects each possible pair of categories and then simply runs a Mann-Whitney U test, using only those two categories.
#' 
#' See ts_mann_whitney() for details of the calculations.
#' 
#' The Bonferroni adjustment is simply:
#' \deqn{p_{adj} = \min \left(p \times n_{comp}, 1\right)}
#' \deqn{n_{comp} = \frac{k\times\left(k-1\right)}{2}}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#'  \item \eqn{n_{comp}}, number of comparisons (pairs)
#'  \item \eqn{k}, number of categories
#' }
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_mann_whitney <- function(catField, ordField, categories=NULL, levels=NULL, method="appr", cc=TRUE){
  
  #Create the cross table
  ct = tab_cross(ordField, catField, order1=levels, order2=categories, totals="include")  
  
  #basic counts
  k = ncol(ct)-1
  nLvl = nrow(ct)-1  
  n = ct[nLvl+1, k+1]
  
  #number of pairs
  ncomp = k*(k-1)/2
  
  colNames = colnames(ct)
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 6))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      selCats = c(colNames[i], colNames[j])
      res[resRow,1] = colNames[i]
      res[resRow,2] = colNames[j]
      mwu = ts_mann_whitney(catField, ordField, selCats, levels, method, cc)
      
      res[resRow,3] = mwu[1, 4]
      res[resRow,4] = mwu[1, 5]
      res[resRow,5] = res[resRow,4] * ncomp
      if (res[resRow,5] > 1){
        res[resRow,5] = 1}
      res[resRow,6] = mwu[1, 6]
      
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("category 1","category 2","statistic","p-value","adj. p-value","test")
  return (res)
}