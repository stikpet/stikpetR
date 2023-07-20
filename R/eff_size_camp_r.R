#' Camp r
#' 
#' @description
#' An approximation for the tetrachoric correlation coefficient. 
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' @param method : optional method to use. Either "cureton" (default), "camp1", or "camp2"
#' 
#' @return Camp r
#' 
#' @details
#' Camp (1934, pp. 309) describes the following steps for the calculation:
#' Step 1: If total of column 1 (C1) is less than column 2 (C2), swop the two columns
#' 
#' Step 2: Calculate \eqn{p = \frac{C1}{n}}, \eqn{p_1 = \frac{a}{n}}, and \eqn{p_2 = \frac{c}{C2}}
#' 
#' Step 3: Determine \eqn{z_1}, \eqn{z_2}, \eqn{z_3} as the normal deviate 
#' corresponding to the area \eqn{p_1}, \eqn{p_2}, \eqn{p} resp. (inverse standard normal cumulative distribution)
#' 
#' Step 4: Determine y the normal ordinate corresponding to \eqn{z_3} (the height of the normal distribution)
#' 
#' Step 5: Calculate \eqn{m = \frac{p\times\left(1-p\right)\times\left(z_1 + z_2\right)}{y}}
#' 
#' Step 6: Find phi in a table of phi values
#' 
#' Camp suggested for a very basic approximation to simply use \eqn{\phi=1}.
#' 
#' For a better approximation Camp made the following table:
#' 
#' | p | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
#' |---|-----|-----|-----|-----|-----|
#' |phi|0.637|0.63|0.62|0.60|0.56|
#' 
#' Cureton (1968, p. 241) expanded on this table and produced:
#' 
#' | p | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
#' |-----|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
#' | 0.5 | 0.637 | 0.636 | 0.636 | 0.635 | 0.635 | 0.634 | 0.634 | 0.633 | 0.633 | 0.632 | 0.631 |
#' | 0.6 | 0.631 | 0.631 | 0.630 | 0.629 | 0.628 | 0.627 | 0.626 | 0.625 | 0.624 | 0.622 | 0.621 |
#' | 0.7 | 0.621 | 0.620 | 0.618 | 0.616 | 0.614 | 0.612 | 0.610 | 0.608 | 0.606 | 0.603 | 0.600 |
#' | 0.8 | 0.600 | 0.597 | 0.594 | 0.591 | 0.587 | 0.583 | 0.579 | 0.574 | 0.569 | 0.564 | 0.559 |
#' 
#' Step 7: Calculate \eqn{r_t = \frac{m}{\sqrt{1+\phi\times m^2}}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell of the cross table 
#' \item \eqn{b} the count in the top-right cell of the cross table 
#' \item \eqn{c} the count in the bottom-left cell of the cross table 
#' \item \eqn{d} the count in the bottom-right cell of the cross table 
#' \item \eqn{R_i} the sum of counts in the i-th row 
#' \item \eqn{C_i} the sum of counts in the i-th column 
#' }
#' 
#' Cureton (1968) describes quite a few shortcomings with this approximation, and circumstances when it might be appropriate.
#' 
#' @references 
#' Camp, B. H. (1934). *Mathematical part of elementary statistics*. D.C. Heath and Company, London.
#' 
#' Cureton, E. E. (1968). Tetrachoric correlation by the Camp approximation. *Educational and Psychological Measurement, 28*(2), 239–244. https://doi.org/10.1177/001316446802800202
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_camp_r(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_camp_r <- function(field1, field2, categories1=NULL, categories2=NULL, method=c("cureton", "camp1", "camp2")){
  
  if (length(method)>1){method="cureton"}
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #the column totals
  colTots <- margin.table(ct, 2)
  C1 <- as.numeric(unname(colTots[1]))
  C2 <- as.numeric(unname(colTots[2]))
  
  #grand total
  n <- sum(colTots)
  
  #create duplicates of original
  C1a = C1
  C2a = C2
  ca=a
  cb=b
  cc=c
  cd=d
  
  #step 1: if C1 < C2 swop the two columns
  switch=1
  if (C1 < C2){
    a2 = ca
    ca = cb
    cb = a2
    c2 = cc
    cc = cd
    cd = c2
    C1a = ca + cc
    C2a = cb + cd
    switch=-1
  }
  
  #step 2: determine three proportions
  p1 = ca/C1a
  p2 = cd/C2a
  p = C1a/n
  
  #step 3: determine the corresponding z-values
  z1 = qnorm(p1)
  z2 = qnorm(p2)
  z3 = qnorm(p)
  
  #step 4: determine the height of the normal distribution
  y = dnorm(z3)
  
  #step 5: calculate m
  m = p*(1-p)*(z1+z2)/y
  
  #step 6: find phi in table of phi-values
  if (method=="camp1"){
    phi = 1
  }
  else if (method=="camp2"){
    phiValues = data.frame("0.5" = 0.637, "0.6" = 0.63, "0.7" = 0.62, "0.8" = 0.6, "0.9" = 0.56)
    phi = phiValues[paste("X" ,floor(p*10)/10, sep="")]
    
  }
  else {
    phiValues = data.frame("0.5" = c(0.637, 0.636, 0.636, 0.635, 0.635, 0.634, 0.634, 0.633, 0.633, 0.632, 0.631), "0.6" = c(0.631, 0.631, 0.63, 0.629, 0.628, 0.627, 0.626, 0.625, 0.624, 0.622, 0.621), "0.7" = c(0.621, 0.62, 0.618, 0.616, 0.614, 0.612, 0.61, 0.608, 0.606, 0.603, 0.6), "0.8" = c(0.6, 0.597, 0.594, 0.591, 0.587, 0.583, 0.579, 0.574, 0.569, 0.564, 0.559))
    phi = phiValues[round(p*100,0)-floor(p*10)*10+1, paste("X" ,floor(p*10)/10, sep="")]
  }
  
  #step 7: calculate r
  r = switch * m/(1+phi*m**2)**0.5
  
  return(r)
}