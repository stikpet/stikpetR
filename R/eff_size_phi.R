#' Pearson/Yule Phi Coefficient / Cole C2 / Mean Square Contingency
#' 
#' @description
#' After performing chi-square test the question of the effect size comes up. An obvious candidate to use in a measure of effect size is the test statistic, the \eqn{\chi^2}. One of the earliest and often mentioned measure uses this: the phi coefficient (or mean square contingency). Both Yule (1912, p. 596) and Pearson (1900, p. 12) mention this measure, and Cole (1949, p. 415) refers to it as Cole C2. It is also the same as Cohen's w (Cohen, 1988, p. 216), but Cohen does not specify it to be only for 2x2 tables.
#' 
#' It is interesting that this gives the same result, as if you would assign a 0 and 1 to each of the two variables categories, and calculate the regular correlation coefficient.
#' 
#' Pearson (1904, p. 6) calls the squared value (i.e. not taking the square root) the Mean Square Contingency.
#' 
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with order for categories of field1
#' @param categories2 : optional list with order for categories of field2
#' 
#' @return phi coefficient
#' 
#' @details
#' 
#' The formula used is (Pearson, 1900, p. 12):
#' \deqn{\phi = \frac{a\times d - b\times c}{\sqrt{R_1\times R_2 \times C_1 \times C_2}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' The formula is also sometimes expressed with a \eqn{\chi^2} value (Pearson, 1904, p.6; Cohen, 1988, p. 216):
#' 
#' \deqn{\phi = \sqrt{\frac{\chi^2}{n}}}
#' 
#' Note that Cohen w did not limit the size of the table, but uses the same formula.
#' 
#' @seealso
#' \code{\link{th_cohen_w}}, rules of thumb for Cohen w
#'
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411-424. https://doi.org/10.2307/1932444
#' 
#' Pearson, K. (1900). Mathematical Contributions to the Theory of Evolution. VII. On the Correlation of Characters not Quantitatively Measurable. *Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character*, 195, 1-405.
#' 
#' Pearson, K. (1904). *Contributions to the Mathematical Theory of Evolution. XIII. On the theory of contingency and its relation to association and normal correlation*. Dulau and Co.
#' 
#' Yule, G. U. (1912). On the methods of measuring association between two attributes. *Journal of the Royal Statistical Society, 75*(6), 579-652. https://doi.org/10.2307/2340126
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_phi(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_phi <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #the row totals
  rowTots <- margin.table(ct, 1)
  R1 <- as.numeric(unname(rowTots[1]))
  R2 <- as.numeric(unname(rowTots[2]))
  
  #the column totals
  colTots <- margin.table(ct, 2)
  C1 <- as.numeric(unname(colTots[1]))
  C2 <- as.numeric(unname(colTots[2]))
  
  phi =(a*d - b*c)/(R1*R2*C1*C2)**0.5
  
  return(phi)
  
}



