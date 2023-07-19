#' Becker and Clogg rho
#' 
#' @description
#' An approximation for the tetrachoric correlation coefficient. 
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' @param version c(1, 2) optional, which version to calculate (see details), default ver=1
#' 
#' @return Becker and Clogg r
#' 
#' @details 
#' Version 1 will calculate:
#' \deqn{\rho^* = \frac{g-1}{g+1}}
#' 
#' Version 2 will calculate:
#' \deqn{\rho^{**} = \frac{OR^{13.3/\Delta} - 1}{OR^{13.3/\Delta} + 1}}
#' 
#' With:
#' \deqn{g=e^{12.4\times\phi - 24.6\times\phi^3}}
#' \deqn{\phi = \frac{\ln\left(OR\right)}{\Delta}}
#' \deqn{OR=\frac{\left(\frac{a}{c}\right)}{\left(\frac{b}{d}\right)} = \frac{a\times d}{b\times c}}
#' \deqn{\Delta = \left(\mu_{R1} - \mu_{R2}\right) \times \left(v_{C1} - v_{C2}\right)}
#' \deqn{\mu_{R1} = \frac{-e^{-\frac{t_r^2}{2}}}{p_{R1}}, \mu_{R2} = \frac{e^{-\frac{t_r^2}{2}}}{p_{R2}}}
#' \deqn{v_{C1} = \frac{-e^{-\frac{t_c^2}{2}}}{p_{C1}}, v_{C2} = \frac{e^{-\frac{t_c^2}{2}}}{p_{C2}}}
#' \deqn{t_r = \Phi^{-1}\left(p_{R1}\right), t_c = \Phi^{-1}\left(p_{C1}\right)}
#' \deqn{p_{x} = \frac{x}{n}}
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
#' \item \eqn{\Phi^{-1}\left(x\right)} for the the inverse standard normal cumulative distribution function
#' }
#' 
#' These formulas can be found in Becker and Clogg (1988, pp. 410-412)
#' 
#' @references 
#' Becker, M. P., & Clogg, C. C. (1988). A note on approximating correlations from Odds Ratios. *Sociological Methods & Research, 16*(3), 407–424. https://doi.org/10.1177/0049124188016003003
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_becker_clogg_r(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_becker_clogg_r <- function(field1, field2, categories1=NULL, categories2=NULL, version=1){
  
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
  
  
  pR1 = R1/n
  pR2 = R2/n
  pC1 = C1/n
  pC2 = C2/n
  
  tr = qnorm(pR1)
  tc = qnorm(pC1)
  
  mR1 = -exp(-tr**2/2) / pR1
  mR2 = exp(-tr**2/2) / pR2
  
  vC1 = -exp(-tc**2/2) / pC1
  vC2 = exp(-tc**2/2) / pC2
  
  delta = (mR1 - mR2)*(vC1 - vC2)
  
  OR = a*d/(b*c)
  
  if (version==2) {
    rt = (OR**(13.3/delta) - 1) / (OR**(13.3/delta) + 1)}
  else if (version==1){
    phiBC = log(OR) / delta
    
    g = exp(12.4*phiBC - 24.6*phiBC**3)
    
    rt = (g - 1)/(g + 1)}
  
  return(rt)
  
}