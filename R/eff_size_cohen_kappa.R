#' Cohen Kappa
#' 
#' @param nom1 the scores on the first variable
#' @param nom2 the scores on the second variable
#' @param ase c("exact", "approximate") optional to indicate which method to use to calculate asymptotic standard errors
#' @return dataframe with the effect size value, the asymptotic standard errors (assuming null and alternative), test statistic, degrees of freedom, and p-value (sig.)
#' 
#' @details 
#' The formula used is (Cohen, 1960, p. 40):
#' \deqn{\kappa = \frac{p_0 - p_c}{1 - p_c}}
#' 
#' With:
#' \deqn{p_0 = \frac{P}{n}}
#' \deqn{p_c = \frac{Q}{n^2}}
#' \deqn{P = \sum_{i=1}^r F_{i,i}}
#' \deqn{Q = \sum_{i=1}^r R_i \times C_i}
#' \deqn{R_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^r F_{i,j}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{r} is the number of rows (categories in the first variable)
#' \item \eqn{c} is the number of columns (categories in the second variable)
#' \item \eqn{n} is the total number of scores
#' \item \eqn{F_{i,j}} is the frequency (count) of scores equal to the i-th category in the first variable, and the j-th category in the second.
#' }
#' 
#' The approximate asymptotic standard errors (ase="approximate") use (Cohen, 1960, pp. 43-44):
#' \deqn{ASE_0 \approx \sqrt{\frac{p_c}{n\times \left(1 - p_c\right)}}}
#' \deqn{ASE_1 \approx \sqrt{\frac{p_0 \times \left(1 - p_0\right)}{n\times \left(1 - p_c\right)^2}}}
#' 
#' The exact asymptotic standard errors (ase="exact") use (Fleiss et al., 1969, p. 325):
#' \deqn{ASE_0 = \sqrt{\frac{SS_0}{n\times\left(1 - p_c\right)^2}}}
#' \deqn{ASE_1 = \sqrt{\frac{SS_1}{n\times\left(1 - p_c\right)^4}}}
#' With:
#' \deqn{SS_0 = \left(\sum_{i=1}^r p_{i,.} \times p_{.,i} \times \left(1 - \left(p_{.,i} + p_{i,.}\right)\right)^2\right) - p_c^2 + \left(1 - p_0\right)^2 \times \sum_{i=1}^r \underset{j\neq i}{\sum_{j=1}^c} p_{i,.} \times p_{.,j} \times \left(p_{.,i} + p_{j,.}\right)^2}
#' \deqn{SS_1 = \left(\sum_{i=1}^r p_{i,i} \times \left(\left(1 - p_c\right) - \left(p_{.,i} + p_{i,.}\right)\times\left(1 - p_0\right)\right)^2\right) - \left(p_0\times p_c - 2\times p_c + p_0\right)^2 + \left(1 - p_0\right)^2 \times \sum_{i=1}^r \underset{j\neq i}{\sum_{j=1}^c} p_{i,j}\times\left(p_{.,i} + p_{j,.}\right)^2}
#' \deqn{p_{i,j} = \frac{F_{i,j}}{n}}
#' \deqn{p_{i,.} = \frac{R_i}{n}}
#' \deqn{p_{.,j} = \frac{C_j}{n}}
#' 
#' The test is then performed using (Cohen, 1960, p. 44):
#' \deqn{z_{\kappa} = \frac{\kappa}{ASE_0}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{\kappa}\right|\right)\right)}
#' 
#' Where \eqn{\Phi\left(\dots\right)} is the cumulative density function of the standard normal distribution
#' 
#' @references 
#' Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement, 20*(1), 37–46. https://doi.org/10.1177/001316446002000104
#' 
#' Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample standard errors of kappa and weighted kappa. *Psychological Bulletin, 72*(5), 323–327. https://doi.org/10.1037/h0028106
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' nom1 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
#' nom2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
#' es_cohen_kappa(nom1, nom2)
#' es_cohen_kappa(nom1, nom2, ase="approximate")
#' 
#' @export
es_cohen_kappa <- function(nom1, nom2, ase="exact"){
  
  dFra = na.omit(data.frame(nom1, nom2))
  ct = table(dFra)
  cNames = colnames(ct)
  rNames = rownames(ct)
  #The row and column names might not match, so first fix that
  #only use the same row and column names
  ct = ct[rNames %in% cNames, ]
  ct = ct[, cNames %in% rNames]
  #now order the rows and columns
  ct = ct[order(rownames(ct)), ]
  ct = ct[ ,order(colnames(ct))]
  
  r = nrow(ct)
  c = ncol(ct)
  Rs = unname(rowSums(ct))
  Cs = unname(colSums(ct))
  n = sum(Rs)
  
  P = 0
  Q = 0
  for (i in 1:r) {
    P = P + ct[i,i]
    Q = Q + Rs[i]*Cs[i]
  }
  
  p0 = P/n
  pc = Q/n**2
  
  k = (p0 - pc)/(1 - pc)
  
  if (ase == "approximate") {
    #Approximate ASE_1
    ASE_1 = sqrt(p0*(1 - p0) / (n*(1 - pc)^2))
    
    #Approximate ASE_0
    ASE_0 = sqrt(pc/(n*(1 - pc)))
  }
  else if (ase == "exact"){
    #For exact ASE
    pij = ct/n
    pid = Rs/n
    pdj = Cs/n
    SS1t1 = 0
    SS1t2 = 0
    SS0t1 = 0
    SS0t2 = 0
    for (i in 1:r) {
      SS1t1 = SS1t1 + pij[i,i]*((1 - pc) - (pdj[i] + pid[i])*(1 - p0))**2
      SS0t1 = SS0t1 + pid[i] * pdj[i] * (1 - (pdj[i] + pid[i]))**2
      for (j in 1:c) {
        if (i != j) {
          SS1t2 = SS1t2 + pij[i,j]*(pdj[i] + pid[j])**2
          SS0t2 = SS0t2 + pid[i]*pdj[j]*(pdj[i] + pid[j])**2
        }
      }
      
    }
    SS1 = SS1t1 + (1 - p0)**2*SS1t2 - (p0*pc - 2*pc + p0)**2
    SS0 = SS0t1 + SS0t2 - pc**2
    ASE_1 = sqrt(SS1/(n*(1 - pc)**4))
    ASE_0 = sqrt(SS0/(n*(1 - pc)**2))
  }
  
  #Testing
  z = k/ASE_0
  pValue = 2*(1 - pnorm(abs(z)))
  
  results = data.frame(k, ASE_1, ASE_0, z, pValue)
  
  return(results)
  
}