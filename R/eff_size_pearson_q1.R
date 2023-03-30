#' Pearson Q1
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
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
es_pearson_q1 <- function(var1, var2){
  
  data = data.frame(var1, var2)
  
  #remove missing values
  data = na.omit(data)
  
  #Create a cross table first
  dataTable = table(data)

  #store the individual cells
  a = dataTable[1,1]
  b = dataTable[1,2]
  c = dataTable[2,1]
  d = dataTable[2,2]
  
  #Pearson requires for Q1 that ad>bc, a>d, and c>b
  #The original
  ap <- a
  bp <- b
  cp <- c
  dp <- d
  
  sw = -1
  
  #Check once, if not swop two columns
  nSwops = 0
  while(ap*dp<bp*cp || ap < dp || cp < bp){
    #after three swoppings rotate table
    if (nSwops==3){
      at <- ap
      ap <- bp
      bp <- dp
      dp <- cp
      cp <- at
      
      nSwops=0
      
      sw <- -sw
    }
    
    else{
      at <- ap
    
      if (ap*dp<bp*cp) {
        #swop the columns
        ap <- bp
        bp <- at
        
        ct <- cp
        cp <- dp
        dp <- ct
      } else{
        #swop the rows
        ap <- cp
        cp <- at
        
        bt <- bp
        bp <- dp
        dp <- bt
      }
      
      sw <- -sw
      nSwops = nSwops+1
    }
  }
  
  q1 = sw*sin(pi/2 * (cp+dp)*(ap+cp)/((ap+bp) * (bp+dp)))
  
  return(q1)
  
}
