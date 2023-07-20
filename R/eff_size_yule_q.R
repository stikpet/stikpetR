#' Yule Q / Cole C4 / Pearson Q2
#' 
#' @description
#' Yule Q as well as Yule Y measure how much bigger the diagonal top-left to bottom-right is, than top-right to bottom-left. If the two variables are paired it can be seen as the difference between the pairs that are in agreement with the ones that are in disagreement over the total number of pairs.
#' 
#' It ranges from -1 to 1.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Yule Q
#' 
#' @details
#' Yule's Q (1900) is also Cole's C4 (1949, p. 415) and Pearson's Q2 (1900, p. 15).
#' 
#' The formula used (Yule, 1900, p. 272):
#' \deqn{Q = \frac{a\times d - b\times c}{a\times d + b\times c}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' }
#' 
#' @section Alternatives:
#' 
#' The *psych* library has a function that also shows Yule Q: *Yule()*
#' 
#' @seealso
#' \code{\link{th_yule_q}}, rules of thumb for Yule Q
#' 
#' \code{\link{es_convert}}, to convert Yule Q to an Odds Ratio or Yule Y.
#'
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London*, 195, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' Yule, G. U. (1900). On the association of attributes in statistics: With illustrations from the material of the childhood society, &c. *Philosophical Transactions of the Royal Society of London*, 194, 257–319. https://doi.org/10.1098/rsta.1900.0019
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_yule_q(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_yule_q <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  q = (a*d - b*c)/(a*d + b*c)
  
  return(q)
  
}