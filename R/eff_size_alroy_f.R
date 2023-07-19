#' Alroy's Forbes Adjustment
#' 
#' @description
#' A measure of association between two binary variables.
#' 
#' An adjustment to the Forbes coefficient.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Alroy F
#' 
#' @details 
#' The formula used (Alroy, 2015):
#' \deqn{F' = \frac{a\times\left(n' + \sqrt{n'}\right)}{a\times\left(n' + \sqrt{n'}\right) + \frac{3}{2}\times b\times c}}
#' With:
#' \deqn{n' = a + b + c}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' @references 
#' Alroy, J. (2015). A new twist on a very old binary similarity coefficient. *Ecology, 96*(2), 575–586. https://doi.org/10.1890/14-0471.1
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_alroy_f(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_alroy_f <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  nAdj = a + b + c
  
  fadj = a*(nAdj + sqrt(nAdj))/(a*(nAdj + sqrt(nAdj)) + 3/2*b*c)
  
  return(fadj)
}