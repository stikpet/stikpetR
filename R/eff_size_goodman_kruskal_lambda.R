#' Goodman-Kruskal Lambda
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param categories1 optional, categories to use for field1
#' @param categories2 optional, categories to use for field2
#' @param ties c("first", "random", "average") optional to indicate what to do in case of multimodal situations.
#' @return dataframe with the effect size value, asymptotic standard error (both assuming null and alternative), the test statistic, and p-value
#' 
#' @details 
#' The formula used is (Goodman & Kruskal, 1954, p. 743):
#' \deqn{\lambda_{Y|X} = \frac{\left(\sum_{i=1}^r f_{i,max}\right) - C_{max}}{n - C_{max}}}
#' \deqn{\lambda_{X|Y} = \frac{\left(\sum_{j=1}^c f_{max,j}\right) - R_{max}}{n - R_{max}}}
#' \deqn{\lambda = \frac{\left(\sum_{i=1}^r f_{i,max}\right) + \left(\sum_{j=1}^c f_{max,j}\right)- C_{max} - R_{max}}{2\times n - C_{max} - R_{max}}} 
#' 
#' The asymptotic standard errors are calculated using (Hartwig, 1973, p. 308):
#' \deqn{ASE\left(\lambda_{Y|X}\right)_0 = \sqrt{\left(\sum_{i,j}\left(F_{i,j}\times\left(\delta_{i,j}^c - \delta_j^c\right)^2\right)\right) - \frac{\left(\left(\sum_{i=1}^r f_{i, max}\right) - C_{max}\right)^2}{n}}{n - C_{max}}}
#' \deqn{ASE\left(\lambda_{X|Y}\right)_0 = \sqrt{\left(\sum_{i,j}\left(F_{i,j}\times\left(\delta_{i,j}^r - \delta_i^r\right)^2\right)\right) - \frac{\left(\left(\sum_{j=1}^c f_{max, j}\right) - R_{max}\right)^2}{n}}{n - R_{max}}}
#' \deqn{ASE\left(\lambda\right)_0 = \sqrt{\left(\sum_{i,j}\left(F_{i,j}\times\left(\delta_{i,j}^c - \delta_j^c + \delta_{i,j}^r - \delta_i^r\right)^2\right)\right) - \frac{\left(\left(\sum_{i=1}^r f_{i, max}\right) - C_{max} + \left(\sum_{j=1}^c f_{max, j}\right) - R_{max}\right)^2}{n}}{2\times n - C_{max} - R_{max}}}
#' 
#' \deqn{ASE\left(\lambda_{Y|X}\right)_1 = \sqrt{\frac{\left(n - \sum_{i=1}^r f_{i, max}\right) \times \left(\left(\sum_{i=1}^r f_{i, max}\right) + C_{max} - 2\times\sum_{i,j}\left(F_{i,j}\times\delta_{i,j}^c\times\delta_j^c\right)\right)}{\left(n - C_{max}\right)^3}}}
#' \deqn{ASE\left(\lambda_{X|Y}\right)_1 = \sqrt{\frac{\left(n - \sum_{j=1}^c f_{max, j}\right) \times \left(\left(\sum_{j=1}^c f_{max, j}\right) + R_{max} - 2\times\sum_{i,j}\left(F_{i,j}\times\delta_{i,j}^r\times\delta_i^r\right)\right)}{\left(n - R_{max}\right)^3}}}
#' \deqn{ASE\left(\lambda\right)_1 = \frac{\left(\sum_{i,j}\left(F_{i,j}\times\left(\delta_{i,j}^c + \delta_{i,j}^r - \delta_{j}^c - \delta_{i}^r + \lambda\times\left(\delta_{j}^c + \delta_{i}^r\right)\right)^2\right)\right) - 4\times n\times\lambda^2}{2\times n - C_{max} - R_{max}}}
#' 
#' With:
#' \deqn{\delta_{i,j}^c \begin{cases} 1 & \text{ if } j = \text{column index} f_{i,max} \\ 0 & \text{ else } \end{cases}}
#' \deqn{\delta_{j}^c \begin{cases} 1 & \text{ if } j = \text{column index} C_{max} \\ 0 & \text{ else } \end{cases}}
#' \deqn{\delta_{i,j}^r \begin{cases} 1 & \text{ if } i = \text{row index} f_{max,j} \\ 0 & \text{ else } \end{cases}}
#' \deqn{\delta_{i}^r \begin{cases} 1 & \text{ if } i = \text{column index} R_{max} \\ 0 & \text{ else } \end{cases}}
#' 
#' The test is performed using:
#' \deqn{z_{i} = \frac{\lambda_{i}}{ASE\left(i\right)_0}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_i\right|\right)\right)}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{F_{i,j}} the absolute frequency (observed count) from row i and column j.
#' \item \eqn{c} the number of columns
#' \item \eqn{r} the number of rows
#' \item \eqn{R_i} row total of row i, it can be calculated using \eqn{R_i=\sum_{j=1}^c F_{i,j}}
#' \item \eqn{C_j} column total of column j, it can be calculated using \eqn{C_j=\sum_{i=1}^r F_{i,j}}
#' \item \eqn{n} the total number of cases, it can be calculated in various ways, \eqn{n = \sum_{j=1}^c C_j =\sum_{i=1}^r R_i = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \item \eqn{f_{i,max}} is the maximum count of row i, i.e. \eqn{f_{i,max}=\text{max}\left\{ F_{i,1},F_{i,2},…,F_{i,c}\right\}}
#' \item \eqn{f_{max,j}} is the maximum count of column j, i.e. \eqn{f_{max,j}=\text{max}\left\{ F_{1,j},F_{2,j},…,F_{r,j}\right\}}
#' \item \eqn{R_{max}} is the maximum of the row totals, i.e. \eqn{R_{max}=\text{max}\left\{ R_1,R_2,…,R_r\right\}}
#' \item \eqn{C_{max}} is the maximum of the column totals, i.e. \eqn{C_{max}=\text{max}\left\{ C_1,C_2,…,C_c\right\}}
#' \item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' Unfortunately not much is written about how to deal with situations if more than one row / column / cell
#' has the highest (i.e. a multimodal situation). 
#' Hartwig proposed three options in case multi-modal situation occurs: choose random, 
#' choose the largest ASE, or average them.
#' This function can allow you to simply choose the first only (I think SPSS uses this),
#' average them, or simply choose one at random.
#' 
#' @references 
#' Goodman, L. A., & Kruskal, W. H. (1954). Measures of Association for Cross Classifications. *Journal of the American Statistical Association, 49*(268), 732–764. https://doi.org/10.2307/2281536
#' 
#' Gray, L. N., & Campbell, R. (1975). Statistical significance of the Lambda coefficients: A comment. *Behavioral Science, 20*(4), 258–259. https://doi.org/10.1002/bs.3830200407
#' 
#' Hartwig, F. (1973). Statistical significance of the lambda coefficients. *Behavioral Science, 18*(4), 307–310. https://doi.org/10.1002/bs.3830180409
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
es_goodman_kruskal_lambda <- function(field1, field2, categories1=NULL, categories2=NULL, ties="first"){
  
  #The average method only averages column and row maximums (rm and cm),
  #it uses "first" for ties in fim and fmj.
  
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
  
  Rs = rowSums(ct)
  Cs = colSums(ct)
  n = sum(Rs)
  
  rm = max(Rs)
  cm = max(Cs)
  
  ASE_0s = c(0)
  ASE_0sXY = c(0)
  ASE_0sYX = c(0)
  
  ASE_1s = c(0)
  ASE_1sXY = c(0)
  ASE_1sYX = c(0)  
  nRuns = 0
  RsLoop = TRUE
  RsLoopIndex = 1
  while (RsLoop) {
    RsIndex = seq(along=Rs)[Rs == rm]
    tNrm = length(RsIndex)
    if (ties=="first") {
      RsUniqueIndex = RsIndex[1]
      RsLoop = FALSE
    }
    else if (ties=="random") {
      RsUniqueIndex = RsIndex[sample(1:tNrm,1)]
      RsLoop = FALSE
    }
    
    else if(ties=="average"){
      RsUniqueIndex = RsIndex[RsLoopIndex]
      if (tNrm==RsLoopIndex) {
        RsLoop = FALSE
      }
      else{
        RsLoopIndex = RsLoopIndex + 1
      }
    }
    
    CsLoop = TRUE
    CsLoopIndex = 1
    
    while (CsLoop){
      
      nRuns = nRuns + 1
      
      CsIndex = seq(along=Cs)[Cs == cm]
      tNcm = length(CsIndex)
      if (ties=="first") {
        CsUniqueIndex = CsIndex[1]
        CsLoop = FALSE
      }
      else if (ties=="random") {
        CsUniqueIndex = CsIndex[sample(1:tNcm,1)]
        CsLoop = FALSE
      }
      else if(ties=="average"){
        CsUniqueIndex = CsIndex[CsLoopIndex]
        if (tNcm==CsLoopIndex) {
          CsLoop = FALSE
        }
        else{
          CsLoopIndex = CsLoopIndex + 1
        }
      }
      
      #fim is the maximum count of row i
      fim = rep(0, r)
      fimIndex = list()
      fimUniqueIndex = rep(0, r)
      for (i in 1:r) {
        fim[i] = max(ct[i,])
        fimIndex[[i]] = seq(along=ct[i,])[ct[i,]==fim[i]]
        if (ties=="first" || ties=="average") {
          fimUniqueIndex[i] = fimIndex[[i]][1]
        }
        else if (ties=="random") {
          fimUniqueIndex[i] = fimIndex[[i]][sample(1:length(fimIndex[[i]]),1)]
        }
      }
      
      dijc = matrix(0, nrow=r, ncol=c)
      djc = matrix(0, nrow=r, ncol=c)
      for (i in 1:r) {
        for (j in 1:c) {
          if (j==fimUniqueIndex[i]) {
            dijc[i,j] = 1
          }
          if (j == CsUniqueIndex) {
            djc[i,j] = 1
          }
        }
      }
      
      #fmj is the maximum count of column j
      fmj = rep(0, c)
      fmjIndex = list()
      fmjUniqueIndex = rep(0,c)
      for (i in 1:c) {
        fmj[i] = max(ct[,i])
        fmjIndex[[i]] = seq(along=ct[,i])[ct[,i]==max(ct[,i])]
        
        if (ties=="first" || ties=="average") {
          fmjUniqueIndex[i] = fmjIndex[[i]][1]
        }
        else if (ties=="random") {
          fmjUniqueIndex[i] = fmjIndex[[i]][sample(1:length(fmjIndex[[i]]),1)]
        }
      }
      
      dijr = matrix(0, nrow=r, ncol = c)
      djr = matrix(0, nrow=r, ncol = c)
      for (i in 1:r) {
        for (j in 1:c) {
          if (i==fmjUniqueIndex[j]) {
            dijr[i,j] = 1
          }
          if (i == RsUniqueIndex) {
            djr[i,j] = 1
          }
        }
      }
      
      
      
      Lyx = (sum(fim) - cm)/(n - cm)
      L = Lyx
      
      ASE_Lyx_0 = sum(ct*(dijc - djc)**2)
      ASE_Lyx_1 = sum(ct*dijc*djc)
      
      ASE_Lyx_0 = sqrt(ASE_Lyx_0 - (sum(fim)-cm)**2/n)/(n - cm)
      ASE_Lyx_1 = sqrt((n - sum(fim))*(sum(fim)+cm - 2*ASE_Lyx_1)/(n - cm)**3)
      
      ASE_0sYX[nRuns] = ASE_Lyx_0
      ASE_1sYX[nRuns] = ASE_Lyx_1
      
      
      
      Lxy = (sum(fmj) - rm)/(n - rm)
      L = Lxy
      
      ASE_Lxy_0 = sum(ct*(dijr - djr)**2)
      ASE_Lxy_1 = sum(ct*dijr*djr)
      
      ASE_Lxy_0 = sqrt(ASE_Lxy_0 - (sum(fmj)-rm)**2/n)/(n - rm)
      ASE_Lxy_1 = sqrt((n - sum(fmj))*(sum(fmj)+rm - 2*ASE_Lxy_1)/(n - rm)**3)
      
      ASE_0sXY[nRuns] = ASE_Lxy_0
      ASE_1sXY[nRuns] = ASE_Lxy_1
      
      
      
      L = (sum(fim) + sum(fmj) - cm - rm)/(2*n - rm - cm)
      
      ASE_0b = sum(ct*(dijr + dijc - djr - djc)**2)
      ASE_1b = sum(ct*(dijr+dijc-djr-djc+L*(djr+djc))**2)
      
      ASE_0s[nRuns] = sqrt(ASE_0b - (sum(fim)+sum(fmj)-rm-cm)**2/n)/(2*n - rm - cm)
      ASE_1s[nRuns] = sqrt(ASE_1b - 4*n*L**2)/(2*n - rm - cm)
      
    }
  }
  
  if (ties=="max") {
    ASE_0 = c(max(ASE_0s), max(ASE_0sXY), max(ASE_0sYX))
    ASE_1 = c(max(ASE_1s), max(ASE_1sXY), max(ASE_1sYX))
  }
  else if(ties=="average"){
    ASE_0 = c(mean(ASE_0s), mean(ASE_0sXY), mean(ASE_0sYX))
    ASE_1 = c(mean(ASE_1s), mean(ASE_1sXY), mean(ASE_1sYX))
  }
  else {
    ASE_0 = c(ASE_0s[1], ASE_0sXY[1], ASE_0sYX[1])
    ASE_1 = c(ASE_1s[1], ASE_1sXY[1], ASE_1sYX[1])
  }
  
  dependent = c("symmetric", "field1", "field2")
  Ls = c(L, Lxy, Lyx)    
  statistic = Ls/ASE_0
  p = 2*(1-pnorm(abs(statistic)))
  
  results = data.frame(dependent, Ls, n, ASE_0, ASE_1, statistic, p)
  
  return(results)
  
}