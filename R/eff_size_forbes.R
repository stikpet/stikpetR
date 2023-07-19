#' Forbes Coefficient
#' 
#' @description
#' An effect size measure for two binary variables.
#' 
#' A measure for the association between two binary variables. If all values are equal in the cross table, there is no association. Forbes uses that if all values are the same, then:
#' \deqn{1 = \frac{\left(a+b\right)\times\left(a+c\right)}{n}}
#' 
#' While the above equation would result in 0 or 2 if there is a perfect association (i.e. b or c is 0).
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Forbes Coefficient
#' 
#' @details 
#' 
#' The formula used (Forbes, 1907, p. 279):
#' \deqn{F = \frac{n\times\min\left(a, d\right)}{C_1\times R_1}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_1} the sum of counts in the 1st row 
#' \item \eqn{C_1} the sum of counts in the 1st column 
#' }
#' 
#' @references 
#' Forbes, S. A. (1907). On the local distribution of certain Illinois fishes: An essay in statistical ecology. *Illinois Natural History Survey Bulletin, 7*(8), 273–303.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_forbes(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_forbes <- function(field1, field2, categories1=NULL, categories2=NULL){
  
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
  
  #grand total
  n <- sum(colTots)
  
  forb = n*min(a,d)/(C1 * R1)
  
  return(forb)
}