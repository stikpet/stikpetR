#' Cole C1
#' 
#' @description 
#' A measure of association between two binary variables.
#' 
#' This is an adjustment of the Forbes coefficient (\code{\link{es_forbes}}), to convert the scale between -1 and 1, by simply subtracting 1 from Forbes F.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Cole C1
#' 
#' @details 
#' The formula used (Cole, 1949, p. 415):
#' \deqn{C_1 = \frac{a\times d - b\times c}{\left(a+b\right)\times\left(a+c\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' It is an adjustment of Forbes F so it would range from -1 to 1. 
#' 
#' Note that Cole C2 is the phi coefficient, Cole C3 is McEwen-Michael coefficient, Cole C4 is Yule Q, and Cole C6 is Yule r.
#' 
#'
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_cole_c1(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_cole_c1 <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  c1 = (a*d - b*c)/((a+b) * (a+c))
  
  return(c1)
}