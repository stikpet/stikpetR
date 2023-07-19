#' Bonett and Price rho
#' @description
#' An approximation for the tetrachoric correlation coefficient. 
#' 
#' @param var1 A vector with the binary data from the first variable
#' @param var2 A vector with the binary data from the second variable
#' @param version c(1, 2) Optional either version 1 for rho, or version 2 (default) for rho-hat
#' @return Bonett and Price r
#' 
#' @details
#' An approximation for the tetrachoric correlation coefficient.
#' 
#' Formula for version 1 is (Bonett & Price, 2005, p. 216):
#' \deqn{\rho^* = \cos\left(\frac{\pi}{1+\omega^c}\right)}
#' With:
#' \deqn{\omega = OR = \frac{a\times d}{b\times c}}
#' \deqn{c = \frac{1-\frac{\left|R_1-C_1\right|}{5\times n} - \left(\frac{1}{2}-p_{min}\right)^2}{2}}
#' \deqn{p_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)}{n}}
#' 
#' Formula for version 2 is  (Bonett & Price, 2005, p. 216):
#' \deqn{\hat{\rho}^* = \cos\left(\frac{\pi}{1+\hat{\omega}^{\hat{c}}}\right)}
#' with:
#' \deqn{\hat{\omega} = \frac{\left(a+\frac{1}{2}\right)\times \left(d+\frac{1}{2}\right)}{\left(b+\frac{1}{2}\right)\times \left(c+\frac{1}{2}\right)}}
#' \deqn{\hat{c} = \frac{1-\frac{\left|R_1-C_1\right|}{5\times \left(n+2\right)} - \left(\frac{1}{2}-\hat{p}_{min}\right)^2}{2}}
#' \deqn{\hat{p}_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)+1}{n+2}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' \item \eqn{n} the sum of all counts
#' }
#' 
#' @references 
#' Bonett, D. G., & Price, R. M. (2005). Inferential methods for the tetrachoric correlation coefficient. *Journal of Educational and Behavioral Statistics, 30*(2), 213–225. https://doi.org/10.3102/10769986030002213
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_bonett_price_r(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_bonett_price_r <- function(field1, field2, categories1=NULL, categories2=NULL, version=2){
  
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
  
  n = a+b+c+d
  
  if(version==1){
    pMin = min(min(rowTots),  min(colTots))/n
    cBP = (1 - abs(R1 - C1)/(5*n) - (0.5 - pMin)**2)/2
    omg = a*d/(b*c) 
    r = cos(pi / (1 + omg**cBP))}
  else if(version==2){
    pMin2 = (min(min(rowTots),  min(colTots)) + 1)/(n+2)
    cBP2 = (1 - abs(R1 - C1)/(5*(n+2)) - (0.5 - pMin2)**2)/2
    omg2 = (a+0.5)*(d+0.5)/((b+0.5)*(c+0.5)) 
    r = cos(pi  / (1 + omg2**cBP2))}

  return(r)
  
}