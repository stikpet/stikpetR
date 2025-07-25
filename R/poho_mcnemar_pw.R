#' Post-Hoc McNemar Test - Pairwise
#' @description
#' After a (McNemar-)Bowker test a post-hoc test can potentially locate where the changes occured. This can be done use a McNemar test, which is the Bowker test but for 2x2 tables.
#' 
#' There are two variations, one is to simply compare each possible pair of categories (pairwise comparison), or compare each category with all other categories (collapsed comparison). This function is for the pairwise version, see **ph_mcnemar_co()** for the collapsed version.
#' 
#' Instead of using the McNemar test it is also possible to use the binomial test, which will be used if exact is set to True.
#' 
#' @param field1 vector, the first categorical field
#' @param field2 vector, the first categorical field
#' @param categories vector, optional, order and/or selection for categories of field1 and field2
#' @param exact boolean, optional, use of exact binomial distribution (default is False)
#' @param cc boolean, optional, use of a continuity correction (default is False)
#' 
#' @returns
#' Dataframe with:
#' \item{field1}{the first category compared to the second}
#' \item{field2}{the second category compared to the first}
#' \item{n}{the sample size}
#' \item{statistic}{the chi-squared value (if applicable)}
#' \item{df}{the degrees of freedom used in the test (if applicable)}
#' \item{p-value}{the significance (p-value)}
#' \item{adj. p-value}{the Bonferroni adjusted p-value}
#' 
#' @details 
#' The formula used is (McNemar, 1947, p. 156):
#' \deqn{\chi_{M}^2 = \frac{\left(F_{1,2} - F_{2,1}\right)^2}{F_{1,2} + F_{2,1}}}
#' \deqn{df = 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_M^2, df\right)}
#' If a continuity correction is applied the formula changes to:
#' \deqn{\chi_{M*}^2 = \frac{\left(\left|F_{1,2} - F_{2,1}\right| - 1\right)^2}{F_{1,2} + F_{2,1}}}
#' 
#' The formula used for the binomial test is:
#' \deqn{sig. = 2\times\text{Bin}\left(F_{1,2} + F_{2,1}, \min\left(F_{1,2}, F_{2,1}\right), 0.5\right)}
#' 
#' The formula used for the binomial test with a mid-p correction:
#' \deqn{sig. = 2\times\text{Bin}\left(F_{1,2} + F_{2,1}, \min\left(F_{1,2}, F_{2,1}\right), 0.5\right) -\text{bin}\left(F_{1,2} + F_{2,1}, \min\left(F_{1,2}, F_{2,1}\right), 0.5\right)}
#' 
#' The number of pairwise tests \eqn{n_{comp}} ) is:
#' \deqn{n_{comp} = \frac{k\times \left(k - 1\right)}{2}}
#' 
#' The adjusted p-value is then determined using a Bonferroni correction:
#' \deqn{sig._{adj} = \begin{cases} sig. \times n_{comp} & \text{ if } sig. \times n_{comp} \leq = 1 \\ 1 & \text{ if } sig. \times n_{comp} > 1 \end{cases}}
#' 
#' *Symbols used*
#' \itemize{
#'  \item \eqn{F_{1,2}}, the observed count of cases that scored category 1 on the first variable, and category 2 on the second.
#'  \item \eqn{F_{2,1}}, the observed count of cases that scored category 2 on the first variable, and category 1 on the second.
#'  \item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function for the chi-square distribution.
#'  \item \eqn{\text{Bin}\left(\dots\right)}, the cumulative distribution function for the binomial distribution.
#'  \item \eqn{\text{bin}\left(\dots\right)}, the probability mass function for the binomial distribution.
#' }
#' 
#' @references
#' McNemar, Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika, 12*(2), 153-157. doi:10.1007/BF02295996
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_mcnemar_pw <- function(field1, field2, categories=NULL, exact=FALSE, cc=FALSE){
  ct = tab_cross(field1, field2, categories, categories, totals="include")  
  
  #basic counts
  k = nrow(ct)-1
  n = ct[k+1, k+1]
  
  #number of pairs
  ncomp = k*(k-1)/2
  
  rowNames = row.names(ct)
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 7))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      cats = c(rowNames[i], rowNames[j])
      res[resRow,1] = rowNames[i]
      res[resRow,2] = rowNames[j]
      
      if (exact){
        tab = tab_cross(field1, field2, order1=cats, order2=cats)
        n = tab[1,2] + tab[2,1]
        minCount = min(tab[1,2], tab[2,1])
        pVal = pbinom(minCount, n,0.5)*2
        if (cc){
          pVal = pVal - dbinom(minCount, n,0.5)}
        
        stat = NA
        df = NA
      }
      
      else{
        resMc = ts_mcnemar_bowker(field1, field2, categories=cats, cc=cc)
        n = resMc[1,1]
        stat = resMc[1,2]
        df = resMc[1,3]
        pVal = resMc[1,4]
      }
      
      res[resRow,3] = n
      res[resRow,4] = stat
      res[resRow,5] = df
      res[resRow,6] = pVal
      res[resRow,7] = res[resRow,6] * ncomp
      if (res[resRow,7] > 1){
        res[resRow,7] = 1}
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("field1", "field2", "n", "statistic", "df", "p-value", "adj. p-value")
  return (res)
}



