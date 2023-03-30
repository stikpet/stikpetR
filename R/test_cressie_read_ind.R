#' Cressie-Read Test of Independence / Power Divergence Test
#' 
#' @param nom1 A vector with the data of the first variable
#' @param nom2 A vector with the data of the second variable
#' @param cc c(NULL, "yates", "pearson", or "williams") Optional continuity correction (default is NULL)
#' @param lambda Optional value for lambda (default is 2/3)
#' @return dataframe with test statistic, degrees of freedom, p-value, minimum expected count, proportion of expected counts below 5, and test used
#' 
#' 
#' @details 
#' The formula used is (Cressie & Read, 1984, p. 442):
#' \deqn{\chi_{C}^{2} = \begin{cases} 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(F_{i,j}\times ln\left(\frac{F_{i,j}}{E_{i,j}}\right)\right) & \text{ if } \lambda=0 \\ 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(E_{i,j}\times ln\left(\frac{E_{i,j}}{F_{i,j}}\right)\right) & \text{ if } \lambda=-1 \\ \frac{2}{\lambda\times\left(\lambda + 1\right)} \times \sum_{i=1}^{r}\sum_{j=1}^{c} F_{i,j}\times\left(\left(\frac{F_{i,j}}{E_{i,j}}\right)^{\lambda} - 1\right) & \text{ else } \end{cases}}
#' \deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{C}^{2},df\right)}
#' 
#' With:
#' \deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
#' \deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^c F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^r F_{i,j}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{r} the number of categories in the first variable (the number of rows)
#' \item \eqn{c} the number of categories in the second variable (the number of columns)
#' \item \eqn{F_{i,j}} the observed count in row i and column j
#' \item \eqn{E_{i,j}} the expected count in row i and column j
#' \item \eqn{R_i} the i-th row total
#' \item \eqn{C_j} the j-th column total
#' \item \eqn{n} the sum of all counts
#' \item \eqn{\chi^2\left(\dots\right)}	the chi-square cumulative density function
#' }
#' 
#' Cressie and Read (1984, p. 463) suggest to use \eqn{\lambda = \frac{2}{3}},  which
#' is therefor the default in this function.
#' 
#' The Yates correction (yates) is calculated using (Yates, 1934, p. 222):
#' 
#' Use instead of \eqn{F_{i,j}} the adjusted version defined by:
#' \deqn{F_{i,j}^\ast = \begin{cases} F_{i,j} - 0.5 & \text{ if } F_{i,j}>E_{i,j}  \\ F_{i,j} & \text{ if } F_{i,j}= E_{i,j}\\ F_{i,j} + 0.5 & \text{ if } F_{i,j}<E_{i,j} \end{cases}}
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right)-1\right) \times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_j}\right)-1\right)}{6\times n\times df}}
#' 	
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. *Journal of the Royal Statistical Society: Series B (Methodological), 46*(3), 440–464. https://doi.org/10.1111/j.2517-6161.1984.tb01318.x
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @examples  
#' nom1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male","male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
#' nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
#' ts_cressie_read_ind(nom1, nom2)
#' ts_cressie_read_ind(nom1, nom2, cc="yates")
#' ts_cressie_read_ind(nom1, nom2, cc="pearson")
#' ts_cressie_read_ind(nom1, nom2, cc="williams")
#'  
#' @export
ts_cressie_read_ind <- function(nom1, nom2, cc=NULL, lambda=2/3){
  
  datFr = na.omit(data.frame(nom1, nom2))
  
  obs = table(datFr)
  
  R = rowSums(obs)
  nr = length(R)
  C = colSums(obs)
  nc = length(C)
  n = sum(R)
  
  expCount = matrix(1, nrow=nr, ncol=nc)
  for (j in 1:nc) {
    for (i in 1:nr){
      expCount[i,j] = R[i]*C[j]/n
    }
  }
  
  if (is.null(cc)==FALSE && cc=="yates") {
    #Yates correction
    for (j in 1:nc) {
      for (i in 1:nr){
        if (obs[i,j] > expCount[i,j]) {
          obs[i,j] = obs[i,j] - 0.5
        }
        else if (obs[i,j] < expCount[i,j]) {
          obs[i,j] = obs[i,j] + 0.5
        }
      }
    }
  }
  
  if (lambda==-1){
    chiVal = 2*sum(expCount*(log(expCount/obs))) 
  }
  else if (lambda!=0) {
    chiVal = 2/(lambda*(lambda + 1))*sum(obs*((obs/expCount)**lambda - 1))
  }
  else{
    chiVal = 2*sum(obs*(log(obs/expCount)))  
  }
  
  if (is.null(cc)==FALSE){
    if (cc=="pearson") {
      #Pearson Correction
      chiVal = (n - 1)/n * chiVal
    }
    else if (cc=="williams"){
      #Williams Correction
      q = 1 + (n*sum(1/R)-1) * (n*sum(1/C)-1) / (6*n*(nr - 1)*(nc - 1))
      chiVal = chiVal/q
    }
  }
  
  df = (nr - 1)*(nc - 1)
  
  pValue = 1 - pchisq(chiVal, df)
  
  minExp = min(expCount)
  propBelow5 = sum(expCount < 5)/(nr*nc)
  
  statistic = chiVal
  testResults <- data.frame(statistic, df, pValue, minExp, propBelow5)
  
  return(testResults)
  
}