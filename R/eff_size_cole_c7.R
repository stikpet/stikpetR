#' Cole C7 / Coefficient of Interspecific Association
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @return Cole C7
#' 
#' @details 
#' The formula used is (Cole, 1949, pp. 420-421). 
#' \deqn{C_7 = \begin{cases} \frac{a\times d - b\times c}{R_1\times C_2} & \text{ if } a\times d\geq b\times c \\ \frac{a\times d - b\times c}{R_1 \times C_1} & \text{ if } a\times d < b\times c \text{ and } a\leq d  \\ \frac{a\times d - b\times c}{R_2\times C_2} & \text{ if } a\times d < b\times c \text{ and } a >  d \end{cases}}
#' 
#' #' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' Cole refers to this as the 'coefficient of interspecific association' 
#' 
#' Note that Cole C2 is the phi coefficient, Cole C3 is McEwen-Michael coefficient, Cole C4 is Yule Q, and Cole C6 is Yule r.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Cole, L. C. (1949). The measurement of interspecific associaton. *Ecology, 30*(4), 411–424. https://doi.org/10.2307/1932444
#' 
#' @examples 
#' bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' es_cole_c7(bin1, bin2)
#' 
#' @export
es_cole_c7 <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  if (a*d >= b*c){
    c7 = (a*d-b*c)/((a+b)*(b+d))
  }else if(a<=d){
    c7 = (a*d-b*c)/((a+b)*(a+c))
  }else{
    c7 = (a*d-b*c)/((b+d)*(c+d))}
  
  return(c7)
}