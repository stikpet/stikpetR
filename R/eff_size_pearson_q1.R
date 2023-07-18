#' Pearson Q1
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Pearson Q1
#' 
#' @details
#' 
#' The formula used (Pearson, 1900, p. 15):
#' \deqn{Q_1 = \sin\left(\frac{\pi}{2} \times \frac{a\times d - b\times c}{\left(a+b\right)\times\left(b+d\right)}\right)}
#' 
#' With:
#' The cross table arranged such that \eqn{a\times d > b\times c}, \eqn{a>d}, and \eqn{c>b}.
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' }
#' 
#' Note that Pearson (1900) stated: "Q1 was found of little service" (p. 16).
#' 
#' Note that Pearson Q2 is the same as Yule Q and Pearson Q3 is the same as Yule r.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London, 195*, 1–405. https://doi.org/10.1098/rsta.1900.0022
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_pearson_q1(bin1, bin2)
#' 
#' @export
es_pearson_q1 <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #Pearson requires for Q1 that ad>bc, a>d, and c>b
  sw = -1
  runs=0
  while(runs < 2){
    #swop the rows
    if (a*d < b*c || a < d || c < b){
      at = a
      a = c
      c = at
      bt = b
      b = d
      d = bt
      sw = -sw
    }
    
    #swop columns
    if (a*d < b*c || a < d || c < b){
      at = a
      a = b
      b = at
      ct = c
      c = d
      d = ct
      sw = -sw
    }
    
    #swop the rows again
    if (a*d < b*c || a < d || c < b){
      at = a
      a = c
      c = at
      bt = b
      b = d
      d = bt
      sw = -sw
    }
    
    #pivot the table
    if (a*d < b*c || a < d || c < b){
      bt = b
      b = c
      c = bt
      sw = -sw
    }
    
    runs = runs + 1    
  }
  
  q1 = sw*sin(pi/2 * (c+d)*(a+c)/((a+b) * (b+d)))
  
  return(q1)
  
}
