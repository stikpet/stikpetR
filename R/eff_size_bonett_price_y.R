#' Bonett and Price Y*
#' 
#' @description
#' A measure of association between two binary variables.
#' 
#' Yule Q and Yule Y can each be written in the format of:
#' \deqn{\frac{OR^x - 1}{OR^x + 1}}
#' With OR being the Odds Ratio. For Yule Q the \eqn{x=1} and for Yule Y \eqn{x=0.5}. Digby (1983, p. 754) showed that Yule’s Q consistently overestimates the association, while Yule’s Y underestimates it It seems that a better approximation might be somewhere between 0.5 and 1 as the power to use on the Odds Ratio. Bonett and Price derived a formula to determine the optimal value for x in each situation.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Bonett and Price Y*
#' 
#' @details 
#' The formula used (Bonett and Price, 2007, pp. 433-434):
#' \deqn{Y* = \frac{\hat{\omega}^x-1}{\hat{\omega}^x+1}}
#' With:
#' \deqn{x = \frac{1}{2}-\left(\frac{1}{2}-p_{min}\right)^2}
#' \deqn{p_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)}{n}}
#' \deqn{\hat{\omega} = \frac{\left(a+0.1\right)\times\left(d+0.1\right)}{\left(b+0.1\right)\times\left(c+0.1\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' Note that \eqn{\hat{\omega}} is a biased corrected version of the Odds Ratio
#' 
#' @references 
#' Bonett, D. G., & Price, R. M. (2007). Statistical inference for generalized yule coefficients in 2 × 2 contingency tables. *Sociological Methods & Research, 35*(3), 429–446. https://doi.org/10.1177/0049124106292358
#'
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_bonett_price_y(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_bonett_price_y <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #the row totals
  rowTots <- rowSums(ct)
  
  #the column totals
  colTots <- colSums(ct)

  #grand total
  n <- sum(colTots)
  
  #smallest proportion of rows and columns
  pMin = min(min(rowTots),  min(colTots))/n
  
  #determine power to use
  pwr = 1/2 - (1/2 - pMin)^2
  
  #odds ratio
  OR = a*d/(b*c)
  #adjust the odds ratio
  ORadj = (a+0.1)*(d+0.1)/((b+0.1)*(c+0.1))
  
  ybp = (ORadj^pwr - 1)/(ORadj^pwr + 1)
  
  return(ybp)
}