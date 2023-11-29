#' Fisher Exact test
#' 
#' @description
#' Perhaps the most commonly used test when you have two binary variables is the Fisher (Exact) Test (Fisher, 1922, 1950). It tests if "the relative proportions of one variable are independent of the second variable; in other words, the proportions at one variable are the same for different values of the second variable" (McDonald, 2014, p. 77).
#' 
#' Note that for a 2x2 table there are quite a lot of different tests. Upton (1982) discusses 24 of them. For larger tables a Fisher-Freeman-Halton Exact Test could be used.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with order for categories of field1
#' @param categories2 : optional list with order for categories of field2
#' 
#' @returns 
#' pval : the two-sided p-value (sig.)
#' 
#' @details
#' The formula used is from Fisher (1950, p. 96):
#' \deqn{p = \sum_{i=a_{min}}^{a_{max}}\begin{cases} p_i & \text{if } p_i \leq p_s \\  0 & \text{ else } \end{cases}}
#' 
#' With:
#' \deqn{p_x = \frac{\binom{R_1}{x}\times \binom{n - R_1}{C_1-x}}{\binom{n}{C_1}}}
#' \deqn{a_{min} = \max\left(0, C_1 + R_1 - n\right)}
#' \deqn{a_{max} = \min\left(R_1, C_1\right)}
#' \deqn{\binom{x}{y}=\frac{x!}{y!\times\left(x-y\right)!}}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{p_s}, the probability of sample cross table, i.e. p_x with x being the upper-left cell of the the cross table from the sample data.
#' \item \eqn{R_1}, is the total of the first row, 
#' \item \eqn{C_1} the total of the first column. 
#' \item \eqn{n}, is the total sample size.    
#' }
#' 
#' The reason for the minimum value of 'a', is first that it cannot be negative, since these are counts. So 0 would be the lowest ever possible. However, once 'a' is set, and the totals are fixed, all other values should also be positive (or zero). The value for 'b' will be if 'a' is 0, it will simply be R1 - a. The value for 'c' is also no issue, this is simply C1 - a. However 'd' might be negative, even if a = 0. The value for 'd' is n - R1 - c. Since c = C1 - a, we get d = n - R1 - C1 + a. But this could be negative if R1 + C1 > n. So, 'a' must be at least C1 + R1 - n.
#' 
#' The maximum for 'a' is simply the minimum of either it's row total, or column total.
#' 
#' Note that \eqn{p_x} is the probability mass function of a hypergeometric distribution.
#' 
#'  
#' @seealso 
#' \code{\link{tab_cross}}, to create a cross-table.
#' 
#' @references
#' Fisher, R. A. (1922). On the Interpretation of \eqn{\chi^2} from Contingency Tables, and the Calculation of P. *Journal of the Royal Statistical Society, 85*(1), 87–94. https://doi.org/10.2307/2340521
#' 
#' Fisher, R. A. (1950). *Statistical methods for research workers* (11th rev.). Oliver and Boyd.
#' 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' Upton, G. J. G. (1982). A comparison of alternative tests for the 2 x 2 comparative trial. *Journal of the Royal Statistical Society. Series A (General), 145*(1), 86–105. https://doi.org/10.2307/2981423
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ts_fisher(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#'  
#' @export
ts_fisher <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  Rs = rowSums(ct)
  Cs = colSums(ct)
  n = sum(Rs)
  
  amin = max(0, Cs[1] + Rs[1] - n)
  amax = min(Rs[1], Cs[1])
  
  pSample = dhyper(ct[1,1], Cs[1], Cs[2], Rs[1])
  
  pVal = 0
  for (i in amin:amax) {
    pFori = dhyper(i, Cs[1], Cs[2], Rs[1])
    
    if (pFori <= pSample) {
      pVal = pVal + pFori
    }
    
  }
  
  return(pVal)
}