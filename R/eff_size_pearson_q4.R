#' Pearson Q4
#' 
#' @description
#' An approximation for the tetrachoric correlation coefficient. 
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Pearson Q4
#' 
#' @details
#' The formula used (Pearson, 1900, p. 16):
#' \deqn{Q_4 = \sin\left(\frac{\pi}{2} \times \frac{1}{1 + \frac{2\times b\times c\times n}{\left(a\times d -b\times c\right) \times \left(b+c\right)}}\right)}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{n} the sum of all counts
#' }
#' 
#' Note that Pearson Q2 is the same as Yule Q and Pearson Q3 is the same as Yule r.
#' 
#' @references 
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London, 195*, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_pearson_q4(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_pearson_q4 <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  n = a+b+c+d
  
  q4 = sin(pi/2 * 1/(1+2*b*c*n/((a*d-b*c)*(b+c))))
  
  return(q4)
  
}