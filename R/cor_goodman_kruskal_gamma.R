#' Goodman-Kruskal Gamma
#' @description 
#' A rank correlation coefficient. It ranges from -1 (perfect negative association) to 1 (perfect positive association). A zero would indicate no correlation at all.
#' 
#' A positive correlation indicates that if someone scored high on the first field, they also likely score high on the second, while a negative correlation would indicate a high score on the first would give a low score on the second.
#' 
#' Alternatives for Gamma are Kendall Tau, Stuart-Kendall Tau and Somers D, but also Spearman rho could be considered.
#' 
#' Gamma looks at so-called discordant and concordant pairs, and ignores tied pairs. Kendall Tau b does the same, but applies a correction for ties. Stuart-Kendall Tau c also, but also takes the size of the table into consideration. Somers d only makes a correction for tied pairs in one of the two directions. Spearman rho is more of a variation on Pearson correlation, but applied to ranks. See Goktas and isci. (2011) for more information on the comparisons.
#' 
#' @param ordField1 the numeric scores of the first variable
#' @param ordField2 the numeric scores of the second variable
#' @param levels1 vector, optional. the categories to use from ordField1
#' @param levels2 vector, optional. the categories to use from ordField2
#' @param ase  optional. Which asymptotic standard error to use. Either "appr" (default), 0, 1
#' @param useRanks boolean, optional. rank the data first or not. Default is False
#' 
#' @returns 
#' A dataframe with:
#' \item{g}{the Goodman-Kruskal Gamma value}
#' \item{statistic}{the z-value used for the test}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Goodman & Kruskal, 1954, p. 749):
#' \deqn{\gamma = \frac{P - Q}{P + Q}}
#' With:
#' \deqn{P = \sum_{i=1}^r \sum_{j=1}^c P_{i,j}}
#' \deqn{Q = \sum_{i=1}^r \sum_{j=1}^c Q_{i,j}}
#' \deqn{P_{i,j} = F_{i,j}\times C_{i,j}}
#' \deqn{Q_{i,j} = F_{i,j}\times D_{i,j}}
#' \deqn{C_{i,j} = \sum_{h<i}\sum_{k<j} F_{h,k} + \sum_{h>i}\sum_{k>j} F_{h,k}}
#' \deqn{D_{i,j} = \sum_{h<i}\sum_{k>j} F_{h,k} + \sum_{h>i}\sum_{k<j} F_{h,k}}
#' 
#' The test can be done with a generic approximation:
#' \deqn{z_{\gamma} = \gamma\times\sqrt{\frac{P+Q}{n\times\left(1-\gamma^2\right)}}}
#' 
#' If we assume the alternative hypothesis we can obtain (Goodman & Kruskal, 1963, p. 324; Goodman & Kruskal, 1972, p. 416; Brown & Benedetti, 1977, p. 310):
#' \deqn{z_{\gamma} = \frac{\gamma}{ASE_1}}
#' \deqn{ASE_1 = \frac{4}{\left(P+Q\right)^2}\times\sqrt{\sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(Q\times C_{i,j}-P\times D_{i,j}\right)^2}}
#' 
#' While if we assume the null hypothesis we can obtain (Brown & Benedetti, 1977, p. 311):
#' \deqn{z_{\gamma} = \frac{\gamma}{ASE_0}}
#' \deqn{ASE_0 = \frac{2}{P+Q}\times\sqrt{\sum_{i=1}^r\sum_{j=1}^c F_{i,j}\times\left(C_{i,j}- D_{i,j}\right)^2 - \frac{\left(P-Q\right)^2}{n}}}
#' 
#' The significance (p-value) in each case is then determined using:
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{gamma}\right|\right)\right)}
#' 
#' *Symbols*
#' 
#' \itemize{
#' \item \eqn{F_{i,j}} the count of scores equal to i in the first variable and j in the second
#' \item \eqn{n} the total sample size
#' \item \eqn{r} the number of unique categories in the first variable (number of rows)
#' \item \eqn{c} the number of unique categories in the second variable (number of columns)
#' \item \eqn{P} double the number of concordant pairs
#' \item \eqn{Q} double the number of discordant pairs
#' }
#' 
#' Note that Kendall \eqn{\tau_a} is the same as Goodman-Kruskall gamma.
#' 
#' **Alternatives**
#' 
#' *library(DescTools)*
#' 
#' GoodmanKruskalGamma(table(ord1, ord2), conf.level=0.95)
#' 
#' *library(MESS)*
#' 
#' gkgamma(table(ord1, ord2))
#' 
#' *library(ryouready)*
#' 
#' ord.gamma(table(ord1, ord2))
#' 
#' @references 
#' Brown, M. B., & Benedetti, J. K. (1977). Sampling behavior of test for correlation in two-way contingency tables. *Journal of the American Statistical Association, 72*(358), 309-315. doi:10.2307/2286793
#' 
#' Goktas, A., & isci, O. (2011). A comparison of the most commonly used measures of association for doubly ordered square contingency tables via simulation. *Advances in Methodology and Statistics, 8*(1). doi:10.51936/milh5641
#' 
#' Goodman, L. A., & Kruskal, W. H. (1954). Measures of association for cross classifications. J*ournal of the American Statistical Association, 49*(268), 732-764. doi:10.2307/2281536
#' 
#' Goodman, L. A., & Kruskal, W. H. (1963). Measures of Association for Cross Classifications III: Approximate Sampling Theory. *Journal of the American Statistical Association, 58*(302), 310-364. doi:10.1080/01621459.1963.10500850
#' 
#' Goodman, L. A., & Kruskal, W. H. (1972). Measures of Association for Cross Classifications IV: Simplification of Asymptotic Variances. *Journal of the American Statistical Association, 67*(338), 415-421. doi:10.1080/01621459.1972.10482401
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
r_goodman_kruskal_gamma <- function(ordField1, ordField2, levels1=NULL, levels2=NULL, ase="appr", useRanks=FALSE){
  
  if (!is.null(levels1)){
    myFieldOrd = factor(ordField1, ordered = TRUE, levels = levels1)
    ordField1 = as.numeric(myFieldOrd)
  }
  
  if (!is.null(levels2)){
    myFieldOrd = factor(ordField2, ordered = TRUE, levels = levels2)
    ordField2 = as.numeric(myFieldOrd)
  }
  
  ct = table(ordField1, ordField2)
  
  nr = nrow(ct)
  nc = ncol(ct)
  
  if (useRanks){
    colnames(ct)=1:nc
    rownames(ct)=1:nr
  }
  
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
  C = topLeft + botRight
  
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
  D = topRight + botLeft
  
  P = sum(ct*C)
  Q = sum(ct*D)
  
  g = (P - Q)/(P + Q)
  
  n = sum(ct)
  
  if (ase=="appr"){
    z = g * ((P + Q) / (n * (1 - g**2)))**0.5 
  }
  else{
    if (ase==0){
      ASE = 4*sqrt(sum(ct*(Q*C - P*D)**2))/(P + Q)**2}
    else{
      ASE = 2*sqrt(sum(ct*(C - D)**2) - (P - Q)**2/n)/(P + Q)}
    z = g/ASE        
    
  }
  
  
  pValue = 2*(1 - pnorm(abs(z)))
  
  statistic = z
  results = data.frame(g, z, pValue)
  colnames(results)=c("gamma", "statistic", "p-value")
  return(results)
  
}




