#' Digby H
#' 
#' @description
#' A measure of association between two binary variables.
#' 
#' Yule Q and Yule Y can each be written in the format of:
#' \deqn{\frac{OR^x - 1}{OR^x + 1}}
#' With OR being the Odds Ratio. For Yule Q the \eqn{x=1} and for Yule Y \eqn{x=0.5}. Digby (1983, p. 754) showed that Yule’s Q consistently overestimates the association, while Yule’s Y underestimates it It seems that a better approximation might be somewhere between 0.5 and 1 as the power to use on the Odds Ratio. Digby proposed to use 0.75
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Digby H
#' 
#' @details 
#' The formula used (Digby, 1983, pp. 753-754):
#' \deqn{H = \frac{\left(a\times d\right)^{3/4} - \left(b\times c\right)^{3/4}}{\left(a\times d\right)^{3/4} + \left(b\times c\right)^{3/4}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' }
#' 
#' @references 
#' Digby, P. G. N. (1983). Approximating the tetrachoric correlation coefficient. *Biometrics, 39*(3), 753–757. https://doi.org/10.2307/2531104
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_digby_h(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_digby_h <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  h = ((a*d)**0.75 - (b*c)**0.75) / ((a*d)**0.75 + (b*c)**0.75)
  
  return(h)
}