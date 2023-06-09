#' Theil U / Uncertainty Coefficient
#' 
#' @param var1 the scores on the first variable
#' @param var2 the scores on the second variable
#' @param dir c("both", "rows", "columns") optional to select which U to return (rows = first variable)
#' @return dataframe with the effect size value, the asymptotic standard errors (assuming null and alternative)
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
#' Theil, H. (1970). On the estimation of relationships involving qualitative variables. *American Journal of Sociology, 76*(1), 103–154. https://doi.org/10.1086/224909
#' 
#' Theil, H. (1972). *Statistical decomposition analysis: With applications in the social and administrative sciences* (Vol. 14). North-Holland Pub. Co.; American Elsevier Pub. Co.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' nom1 <- c("Fully Disagree", "Disagree", "Fully agree", "Neither disagree nor agree", "Agree", "Agree", "Neither disagree nor agree", "Disagree", "Agree", "Agree", "Agree", "Agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Neither disagree nor agree", "Fully agree", "Fully agree", "Fully Disagree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Disagree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully Disagree", "Fully agree", "Agree", "Agree", "Neither disagree nor agree", "Disagree", "Neither disagree nor agree", "Fully agree", "Fully agree", "Disagree", "Disagree", "Neither disagree nor agree", "Disagree", "Agree", "Disagree", "Fully agree", "Fully agree", "Disagree", "Agree", "Disagree", "Neither disagree nor agree", "Fully Disagree")
#' nom2 <- c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' es_theil_u(nom1, nom2, dir="both")
#' es_theil_u(nom1, nom2, dir="rows")
#' es_theil_u(nom1, nom2, dir="columns")
#' 
#' @export
es_theil_u <- function(var1, var2, dir="both"){
  #rows = for var1
  #columns = for var2
  
  
  dFra = na.omit(data.frame(var1, var2))
  
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
  
  if (dir=="columns") {
    U = num/hy
    
    ASE_1 = 0
    p = 0
    for (i in 1:r) {
      for (j in 1:c) {
        if (ct[i,j]>0) {
          ASE_1 = ASE_1 + ct[i,j]*(hy*log(ct[i,j]/Rs[i]) + (hx - hxy)*log(Cs[j]/n))**2
          p = p + ct[i,j]*log(Rs[i]*Cs[j]/(n*ct[i,j]))**2
          
        }
      }
    }
    ASE_1 = sqrt(ASE_1)/(n*hy**2)
    ASE_0 = sqrt(p - n*num**2)/(n*hy)

  }
  else if (dir=="rows") {
    U = num/hx
    
    ASE_1 = 0
    p = 0
    for (i in 1:r) {
      for (j in 1:c) {
        if (ct[i,j]>0) {
          ASE_1 = ASE_1 + ct[i,j]*(hx*log(ct[i,j]/Cs[j]) + (hy - hxy)*log(Rs[i]/n))**2

          p = p + ct[i,j]*log(Rs[i]*Cs[j]/(n*ct[i,j]))**2
          
        }
      }
    }
    ASE_1 = sqrt(ASE_1)/(n*hx**2)

    ASE_0 = sqrt(p - n*num**2)/(n*hx)

  }
  else {
    U = 2*num/(hx + hy)  
    
    ASE_1 = 0
    p = 0
    for (i in 1:r) {
      for (j in 1:c) {
        if (ct[i,j]>0) {
          ASE_1 = ASE_1 + ct[i,j]*(hxy*log(Rs[i]*Cs[j]/n**2) - (hx + hy)*log(ct[i,j]/n))**2
          
          p = p + ct[i,j]*log(Rs[i]*Cs[j]/(n*ct[i,j]))**2
          
        }
      }
    }
    ASE_1 = 2*sqrt(ASE_1)/(n*(hx+hy)**2)
    
    ASE_0 = 2*sqrt(p - n*num**2)/(n*(hx+hy))
    
  }
  
  statistic = U/ASE_0
  
  results = data.frame(U, ASE_0, ASE_1, statistic)
  
  return(results)

}