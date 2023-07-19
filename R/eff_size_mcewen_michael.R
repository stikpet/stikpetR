#' McEwen-Michael Coefficient / Cole C3
#' 
#' @description 
#' A measure of association between two binary variables.
#' 
#' A problem with measures that use the Forbes coefficient or Odds Ratio is that if only one cell is very large compared to the others, or 0, the association will be quite large (close to -1 or 1). Michael and McEwen attempted to overcome this by adding a correction.
#' 
#' Note that Cole (1949, p. 415) refers to this as C3. Cole (1949, p. 417) critized this approach and proposed some alternatives himself.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return McEwen-Michael Coefficient
#' 
#' @details 
#' The formula used (Michael, 1920, p. 57)):
#' \deqn{\frac{a\times d - b\times c}{\left(\frac{a + d}{2}\right)^2 + \left(\frac{b + c}{2}\right)^2}}
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' Michael, E. L. (1920). Marine Ecology and the coefficient of association: A plea in behalf of quantitative biology. *Journal of Ecology, 8*(1), 54–59. https://doi.org/10.2307/2255213
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_mcewen_michael(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_mcewen_michael <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  mm = (a*d - b*c)/(((a+d)/2)**2 + ((b+c)/2)**2)
  
  return(mm)
}