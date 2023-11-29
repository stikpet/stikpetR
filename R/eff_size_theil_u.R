#' Theil U / Uncertainty Coefficient
#' @description 
#' Theil U is a measure of nominal association. According to Wikipedia: "given Y, what fraction of the bits of X can we predict? In this case we can think of X as containing the total information, and of Y as allowing one to predict part of such information." (2022).
#' 
#' The term Theil U can also refer to two completely different measures, often used in forecasting and sometimes referred to as index of inequality.
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param categories1 optional, categories to use for field1
#' @param categories2 optional, categories to use for field2
#' 
#' @returns 
#' dataframe with 
#' \item{dependent}{the field used as dependent variable}
#' \item{n}{the sample size}
#' \item{value}{the Theil U value}
#' \item{ASE_0}{the asymptotic standard error assuming the null hypothesis}
#' \item{ASE_1}{the asymptotic standard error assuming the alternative hypothesis}
#' \item{statistic}{the z-value}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used:
#' \deqn{U_{Y|X} = \frac{H_X + H_Y - H_{XY}}{H_Y}}
#' \deqn{U_{X|Y} = \frac{H_X + H_Y - H_{XY}}{H_X}}
#' \deqn{U = 2\times\frac{H_X + H_Y - H_{XY}}{H_X + H_Y}}
#' 
#' With:
#' \deqn{H_X = -\sum_{i=1}^r \frac{R_i}{n}\ln\left(\frac{R_i}{n}\right)}
#' \deqn{H_Y = -\sum_{j=1}^c \frac{C_j}{n}\ln\left(\frac{C_j}{n}\right)}
#' \deqn{H_X = -\sum_{i=1}^r\sum_{j=1}^c \frac{F_{i,j}}{n}\ln\left(\frac{F_{i,j}}{n}\right), \text{ for }F_{i,j}>0}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{F_{i,j}} the absolute frequency (observed count) from row i and column j.
#' \item \eqn{c} the number of columns
#' \item \eqn{r} the number of rows
#' \item \eqn{R_i} row total of row i, it can be calculated using \eqn{R_i=\sum_{j=1}^c F_{i,j}}
#' \item \eqn{C_j} column total of column j, it can be calculated using \eqn{C_j=\sum_{i=1}^r F_{i,j}}
#' \item \eqn{n} the total number of cases, it can be calculated in various ways, \eqn{n = \sum_{j=1}^c C_j =\sum_{i=1}^r R_i = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' }
#' 
#' The asymptotic standard erros are calculated using:
#' \deqn{ASE\left(U_{Y|X}\right)_1 = \frac{\sqrt{\sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(H_Y \times \ln\left(\frac{F_{i,j}}{R_i}\right) + \left(H_X - H_XY\right)\times\ln\left(\frac{C_j}{n}\right)\right)^2}}{n\times H_Y^2}}
#' \deqn{ASE\left(U_{X|Y}\right)_1 = \frac{\sqrt{\sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(H_X \times \ln\left(\frac{F_{i,j}}{C_j}\right) + \left(H_Y - H_XY\right)\times\ln\left(\frac{R_i}{n}\right)\right)^2}}{n\times H_X^2}}
#' \deqn{ASE\left(U\right)_1 = \frac{\sqrt{\sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(H_XY \times \ln\left(\frac{R_i \times C_j}{n^2}\right) - \left(H_X + H_Y\right)\times\ln\left(\frac{F_{i,j}}{n}\right)\right)^2}}{n\times \left(H_X +H_Y\right)^2}}
#' 
#' \deqn{ASE\left(U_{Y|X}\right)_0 = \frac{\sqrt{P - n\times\left(H_X + H_Y - H_{XY}\right)^2}}{n\times H_Y}}
#' \deqn{ASE\left(U_{X|Y}\right)_0 = \frac{\sqrt{P - n\times\left(H_X + H_Y - H_{XY}\right)^2}}{n\times H_X}}
#' \deqn{ASE\left(U\right)_0 = \frac{2\times\sqrt{P - n\times\left(H_X + H_Y - H_{XY}\right)^2}}{n\times \left(H_X + H_Y\right)}}
#' 
#' With:
#' \deqn{P = \sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(\ln\left(\frac{R_i\times C_j}{n\times F_{i,j}}\right)\right)}
#' 
#' The test statistic is:
#' \deqn{T_i = \frac{U_i}{ASE\left(U_i\right)_0}}
#' 
#' The formula’s were taken from SPSS 15 Algorithms (2006, p. 117), 
#' unclear what the original source is, probably Theil (1970) or Theil (1972)
#' 
#' @references 
#' SPSS. (2006). SPSS 15.0 algorithms.
#' 
#' Theil, H. (1970). On the estimation of relationships involving qualitative variables. *American Journal of Sociology, 76*(1), 103–154. doi:10.1086/224909
#' 
#' Theil, H. (1972). *Statistical decomposition analysis: With applications in the social and administrative sciences* (Vol. 14). North-Holland Pub. Co.; American Elsevier Pub. Co.
#' 
#' Wikipedia. (2022). Uncertainty coefficient. In Wikipedia. https://en.wikipedia.org/w/index.php?title=Uncertainty_coefficient&oldid=1099636947#Definition
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_theil_u <- function(field1, field2, categories1=NULL, categories2=NULL){
  #remove if not in categories
  if (!is.null(categories1)){
    field1[! field1 %in% categories1] = NA
  }
  if (!is.null(categories2)){
    field2[! field2 %in% categories2] = NA
  }
  
  dFra = na.omit(data.frame(field1, field2))
  
  ct = table(dFra)
  
  r = nrow(ct)
  c = ncol(ct)
  
  Rs = unname(rowSums(ct))
  Cs = unname(colSums(ct))
  n = sum(Rs)
  
  hx = -1*sum(Rs/n*log(Rs/n))
  hy = -1*sum(Cs/n*log(Cs/n))
  
  hxy = 0
  for (i in 1:r) {
    for (j in 1:c) {
      if (ct[i,j]>0) {
        hxy = hxy + ct[i,j]/n*log(ct[i,j]/n)
      }
    }
  }
  hxy = -hxy
  
  num = hx + hy - hxy
  
  U = 2*num/(hx + hy)   
  UY = num/hy
  UX = num/hx
  
  
  ASE_1 = 0
  ASE_1Y = 0
  ASE_1X = 0    
  p = 0
  for (i in 1:r) {
    for (j in 1:c) {
      if (ct[i,j]>0) {
        ASE_1 = ASE_1 + ct[i,j]*(hxy*log(Rs[i]*Cs[j]/n**2) - (hx + hy)*log(ct[i,j]/n))**2  
        ASE_1Y = ASE_1Y + ct[i,j]*(hy*log(ct[i,j]/Rs[i]) + (hx - hxy)*log(Cs[j]/n))**2
        ASE_1X = ASE_1X + ct[i,j]*(hx*log(ct[i,j]/Cs[j]) + (hy - hxy)*log(Rs[i]/n))**2  
        p = p + ct[i,j]*log(Rs[i]*Cs[j]/(n*ct[i,j]))**2
        
      }
    }
  }
  ASE_1Y = sqrt(ASE_1Y)/(n*hy**2)
  ASE_0Y = sqrt(p - n*num**2)/(n*hy)
  
  ASE_1X = sqrt(ASE_1X)/(n*hx**2)
  ASE_0X = sqrt(p - n*num**2)/(n*hx)
  
  ASE_1 = 2*sqrt(ASE_1)/(n*(hx+hy)**2)
  ASE_0 = 2*sqrt(p - n*num**2)/(n*(hx+hy))
  
  #results
  
  value = c(U, UX, UY)
  ase0 = c(ASE_0, ASE_0X, ASE_0Y)
  ase1 = c(ASE_1, ASE_1X, ASE_1Y)
  dependent = c("symmetric", "field1", "field2")
  statistic = value/ase0
  pvalue = 2*(1-pnorm(abs(statistic))) 
  results = data.frame(dependent, n, value, ase0, ase1, statistic, pvalue)
  colnames(results)<-c("dependent", "n", "value", "ASE_0", "ASE_1", "statistic", "p-value")
  
  return(results)
  
}