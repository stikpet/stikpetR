#' Post-Hoc Pairwise Mood-Median Test
#' @description
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. A simple Bonferroni adjustment is also made for the multiple testing.
#' 
#' Other post-hoc tests that could be considered are Dunn, Nemenyi, Steel-Dwass, Conover-Iman, or pairwise Mann-Whitney U test.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param test string, optional. The test of independence to use. Options are "pearson" (default), "fisher", "freeman-tukey", "g", "mod-log", "neyman", "power"
#' @param cc string, optional. Method for continuity correction. Options are NULL (default), "yates", "pearson", "williams"
#' @param lambd float or string, optional. Either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3. Only applies to Power Divergence test. Options besides a float are float, "cressie-read", "likelihood-ratio", "mod-log", "pearson", "freeman-tukey", "neyman"
#' 
#' @returns
#' A dataframe with:
#' \item{category 1}{one of the two categories being compared}
#' \item{category 2}{second of the two categories being compared}
#' \item{statistic}{the test statistic}
#' \item{df}{he degrees of freedom, if applicable}
#' \item{p-value}{the p-value (significance)}
#' \item{adj. p-value}{the Bonferroni adjusted p-value}
#' 
#' @details
#' This function selects each possible pair of categories and then simply runs a Mood-Median test, using only those two categories.
#' 
#' See ts_mood_median() for details of the calculations.
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
#' The R function pairwiseMedianTest from the rcompanion package, as well as the pairwise.mood.medtest function from the RVAideMemoire package, produce the same result, but will apply the Pearson continuity correction.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_mood_median <- function(catField, ordField, categories=NULL, levels=NULL, test="pearson", cc=NULL, lambd=2/3){
  
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
  res = data.frame(matrix(nrow = ncomp, ncol = 7))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      selCats = c(colNames[i], colNames[j])
      res[resRow,1] = colNames[i]
      res[resRow,2] = colNames[j]
      
      tstRes = ts_mood_median(catField, ordField, selCats, levels, test, cc, lambd)
      if (test=="fisher"){
        res[resRow, 3] = NA
        res[resRow, 4] = NA
        res[resRow, 5] = tstRes
        res[resRow, 6] = tstRes
        res[resRow, 7] = "Fisher exact"}
      else{
        res[resRow, 3] = tstRes[1,4]
        res[resRow, 4] = tstRes[1,5]
        res[resRow, 5] = tstRes[1,6]
        res[resRow, 6] = tstRes[1,6]
        res[resRow, 7] = tstRes[1,9]}
      
      
      res[resRow,6] = res[resRow,5] * ncomp
      if (res[resRow,6] > 1){
        res[resRow,6] = 1}
      
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("category 1","category 2","statistic", "df", "p-value","adj. p-value","test")
  return (res)
}