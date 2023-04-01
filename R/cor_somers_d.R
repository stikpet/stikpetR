#' Somers' d
#' @param var1 the numeric scores of the first variable
#' @param var2 the numeric scores of the second variable
#' @param dir direction of the measure
#' @returns 
#' A dataframe with:
#' \item{d}{the Sommers d value}
#' \item{ASE_1}{the asymptotic standard error}
#' \item{ASE_0}{the asymptotic standard error, under hypothesis of independence}
#' \item{statistic}{the statistic from the test}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details
#' 
#' **Asymmetric versions**
#' The formula used is given by (Somers, 1962, p. 804):
#' \deqn{d_{y|x} = \frac{P - Q}{D_r}, d_{x|y} = \frac{P - Q}{D_c}}
#' With:
#' \deqn{P = \sum_{i=1}^r \sum_{j=1}^c P_{i,j}}
#' \deqn{Q = \sum_{i=1}^r \sum_{j=1}^c Q_{i,j}}
#' \deqn{P_{i,j} = F_{i,j}\times C_{i,j}}
#' \deqn{Q_{i,j} = F_{i,j}\times D_{i,j}}
#' \deqn{C_{i,j} = \sum_{h<i}\sum_{k<j} F_{h,k} + \sum_{h>i}\sum_{k>j} F_{h,k}}
#' \deqn{D_{i,j} = \sum_{h<i}\sum_{k>j} F_{h,k} + \sum_{h>i}\sum_{k<j} F_{h,k}}
#' \deqn{D_r = n^2 - \sum_{i=1}^r RS_i^2}
#' \deqn{D_c = n^2 - \sum_{j=1}^c CS_i^2}
#' \deqn{RS_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{CS_i = \sum_{i=1}^r F_{i,j}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} is the number of pairs
#' \item \eqn{r} the number of categories in the first variable (i.e. number of rows)
#' \item \eqn{c} the number of categories in the second variable (i.e. number of columns)
#' \item \eqn{F_{h,k}} is the number of cases that scored h for the first variable, and k for the second
#' \item \eqn{P} is double the number of concordant pairs
#' \item \eqn{Q} is double the number of discordant pairs
#' }
#' 
#' Note that \eqn{RS_i} are the frequencies of the scores in the first variable and
#' \eqn{CS_i} are the frequencies of the scores in the second variable.
#' 
#' For testing (SPSS, 2006, p. 121):
#' \deqn{z_{y|x} = \frac{d_{y|x}}{ASE_{d_{y|x}, 0}}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{y|x}\right|\right)\right)}
#' With:
#' \deqn{ASE_{d_{y|x}, 0} = \frac{2}{D_r}\times\sqrt{s}}
#' \deqn{s = \sum_i=1^r \sum_{j=1}^c F_{i,j}\times\left(C_{i,j} - D_{i,j}\right)^2 - \frac{\left(P - Q\right)^2}{n}}
#' 
#' and similar:
#' \deqn{z_{x|y} = \frac{d_{x|y}}{ASE_{d_{x|y}, 0}}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{x|y}\right|\right)\right)}
#' With:
#' \deqn{ASE_{d_{x|y}, 0} = \frac{2}{D_c}\times\sqrt{s}}
#' 
#' **Symmetric**
#' 
#' The symmetric version is the same result as Kendall tau b and calculated using:
#' \deqn{d = \frac{2\times\left(P - Q\right)}{D_r + D_c}}
#' and tested using:
#' \deqn{z_{d} = \frac{d}{ASE_{d, 0}}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{y|x}\right|\right)\right)}
#' With:
#' \deqn{ASE_{d,0} = \frac{4}{D_r + D_c}\times\sqrt{s}}
#' 
#' The function will also calculate the \eqn{ASE_1} which are defined as:
#' \deqn{ASE_{d_{y|x},1} = \frac{2\times\sqrt{\sum_{i=1}^3\sum_{j=1}^c F_{i,j}\times\left(D_r\times\left(C_{i,j} - D_{i,j}\right) - \left(P - Q\right)\times\left(n - RS_i\right)\right)^2}}{D_r^2}}
#' \deqn{ASE_{d_{x|y},1} = \frac{2\times\sqrt{\sum_{i=1}^3\sum_{j=1}^c F_{i,j}\times\left(D_c\times\left(C_{i,j} - D_{i,j}\right) - \left(P - Q\right)\times\left(n - CS_j\right)\right)^2}}{D_c^2}}
#' \deqn{ASE_{d,1} = \frac{2\times ASE_{\tau_b, 1}}{D_r + D_c}\times\sqrt{D_c\times D_c}}
#' With:
#' \deqn{ASE_{\tau_b, 1} = \frac{\sqrt{\sum_{i=1}^3\sum_{j=1}^c F_{i,j}\times\left(2\times\sqrt{D_r\times D_c}\times\left(C_{i,j} - D_{i,j}\right) + \tau_b\times v_{i,j}\right)^2 - n^3\times\tau_b^2\times\left(D_r + D_c\right)^2}}{D_r\times D_c}}
#' \deqn{v_{i,j} = RS_i\times D_c + CS_j\times D_r}
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' SomersDelta(ord1, ord2, direction = "row")
#' 
#' SomersDelta(ord1, ord2, direction = "column")
#' 
#' *library(ryouready)*
#' 
#' ord.somers.d(table(ord1,ord2))
#' 
#' @examples 
#' ord1 = c(5, 3, 3, 4, 3, 4, 3, 4, 4, 4, 5, 3, 1, 3, 2)
#' ord2 = c(5, 3, 3, 3, 3, 3, 5, 4, 3, 4, 5, 3, 2, 5, 2)
#' 
#' r_somers_d(ord1, ord2, dir="both")
#' r_somers_d(ord1, ord2, dir="rows")
#' r_somers_d(ord1, ord2, dir="columns")
#' 
#' @references 
#' Somers, R. H. (1962). A new asymmetric measure of association for ordinal variables. *American Sociological Review, 27*(6), 799–811. https://doi.org/10.2307/2090408
#' 
#' SPSS. (2006). SPSS 15.0 algorithms.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' 
#' @export

r_somers_d <- function(var1, var2, dir=c("both", "rows", "columns")){
  if (length(dir)>1) {
    dir="rows"
  }
  datFrame = na.omit(data.frame(var1, var2))
  ct = table(datFrame$var1, datFrame$var2)
  
  nr = nrow(ct)
  nc = ncol(ct)
  
  topLeft = matrix(0, nr, nc)
  for (iCell in 2:nr) {
    for (jCell in 2:nc){
      obs = 0
      for (i in 1:(iCell - 1)){
        for (j in 1:(jCell - 1)){
          obs = obs + ct[i,j]
        }
      }
      
      topLeft[iCell, jCell] = obs
      
    }
  }
  
  botRight = matrix(0, nr, nc)
  for (iCell in 1:(nr-1)) {
    for (jCell in 1:(nc-1)){
      obs = 0
      for (i in (iCell+1):nr){
        for (j in (jCell+1):nc){
          
          obs = obs + ct[i,j]
          
        }
      }
      
      botRight[iCell, jCell] = obs
      
    }
  }
  
  conc = topLeft + botRight
  
  topRight = matrix(0, nr, nc)
  for (iCell in 2:nr) {
    for (jCell in 1:(nc-1)){
      obs = 0
      for (i in 1:(iCell - 1)){
        for (j in (jCell + 1):nc){
          obs = obs + ct[i,j]
        }
      }
      
      topRight[iCell, jCell] = obs
      
    }
  }
  
  botLeft = matrix(0, nr, nc)
  for (iCell in 1:(nr-1)) {
    for (jCell in 2:nc){
      obs = 0
      for (i in (iCell+1):nr){
        for (j in 1:(jCell - 1)){
          
          obs = obs + ct[i,j]
          
        }
      }
      
      botLeft[iCell, jCell] = obs
      
    }
  }
  
  disc = topRight + botLeft
  Ps = ct*conc
  P = sum(Ps)
  
  Qs = ct*disc
  Q = sum(Qs)
  
  n = sum(ct)
  
  Rs = unname(rowSums(ct))
  Cs = unname(colSums(ct))
  Dr = n**2 - sum(Rs**2)
  Dc = n**2 - sum(Cs**2)
  
  s = sqrt(sum(ct*(conc - disc)**2) - (P - Q)**2/n)
  
  if (dir=="columns") {
    dyx = (P - Q)/Dr
    d = dyx
    
    ASE_yx1 = 0
    for (i in 1:nr) {
      for (j in 1:nc) {
        ASE_yx1 = ASE_yx1 + ct[i,j]*(Dr*(conc[i,j] - disc[i,j]) - (P - Q)*(n - Rs[i]))**2
        
      }
    }
    ASE_yx1 = 2*sqrt(ASE_yx1)/Dr**2
    ASE_1 = ASE_yx1
    
    ASE_0 = 2*s/Dr
    
  }
  
  else if (dir=="rows") {
    dxy = (P - Q)/Dc
    d = dxy
    
    ASE_xy1 = 0
    for (i in 1:nr) {
      for (j in 1:nc) {
        ASE_xy1 = ASE_xy1 + ct[i,j]*(Dc*(conc[i,j] - disc[i,j]) - (P - Q)*(n - Cs[j]))**2
        
      }
    }
    ASE_xy1 = 2*sqrt(ASE_xy1)/Dc**2
    ASE_1 = ASE_xy1
    
    ASE_0 = 2*s/Dc
    
  }
  
  else {
    d = (P - Q)/(0.5*(Dr + Dc))
    
    tau_b = (P - Q)/sqrt(Dr*Dc)
    ASE_1_tau_b = 0
    for (i in 1:nr) {
      for (j in 1:nc) {
        v = Rs[i]*Dc + Cs[j]*Dr
        ASE_1_tau_b = ASE_1_tau_b + ct[i,j]*(2*sqrt(Dr*Dc)*(conc[i,j] - disc[i,j]) + tau_b*v)**2
      }
    }
    
    ASE_1_tau_b = sqrt(ASE_1_tau_b - n**3*tau_b**2*(Dr+Dc)**2)/(Dr*Dc)
    ASE_1 = 2*ASE_1_tau_b/(Dr + Dc)*sqrt(Dr*Dc)
    
    ASE_0 = 4*s/(Dc + Dr)
  }
  
  z = d/ASE_0
  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic = z
  results = data.frame(d, ASE_1, ASE_0, statistic, pValue)
  
  return(results)
  
}