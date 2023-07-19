#' Yule's Y / Coefficient of Colligation
#' 
#' @description
#' A measure of association between two binary variables. This method is based on the Odds Ratio. It converts the sample data cross table to a similar one but having the diagonal values the same, while keeping the same odds ratio as the original.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Yule's Y
#' 
#' @details 
#' The formula used (Yule, 1912, p. 592):
#' \deqn{Y = \frac{\sqrt{a\times d} - \sqrt{b\times c}}{\sqrt{a\times d} + \sqrt{b\times c}}}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' Yule Y can be converted to Yule Q for which rules-of-thumb (*th_yule_q(q)*) can be used, or converted to an Odds Ratio using also with rules of thumb available via *th_odds_ratio(or)*
#' 
#' @seealso 
#' \code{\link{es_convert}}, to convert Yule Y to Yule Q or an Odds Ratio.
#'
#' @references 
#' Yule, G. U. (1912). On the methods of measuring association between two attributes. *Journal of the Royal Statistical Society, 75*(6), 579–652. https://doi.org/10.2307/2340126
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_yule_y(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_yule_y <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #Yule's Y (1912, p. 592) is a further adaptation of Yule's Q into:
  y = ((a*d)**0.5 - (b*c)**0.5) / ((a*d)**0.5 + (b*c)**0.5)
  
  return(y)
}