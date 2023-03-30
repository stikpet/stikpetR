#' James One-Way Test
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @param order the order of the James test to perform (see details)
#' @param alpha the sig. level to use if order=1
#' @param iters boolean to indicate to use iterations to find p-value
#' @param alt boolean to indicate the use of an alternative calculation of degrees of freedom (see details)
#' 
#' @returns 
#' A dataframe with:
#' \item{statistic}{the J-statistic from the test}
#' \item{df}{the degrees of freedom}
#' \item{Jcrit}{critical J value}
#' \item{"reject H0"}{indication to reject H0 or not}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula use for the test statistic (James, 1951, p. 324):
#' \deqn{J = \sum_{j=1}^k w_j\times\bar{x}_j^2 - \frac{\left(\sum_{s=1}^k w_s\times\bar{x}_s\right)^2}{w}=\chi_{Cochran}^2}
#' With:
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' #' *Symbols:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{x_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{w_j} the weight for category j
#' \item \eqn{h_j} the adjusted weight for category j
#' }
#' 
#' For large group size (order=0) the same result as the Cochran test (James, 1951, p. 324):
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(J, df\right)}
#' 
#' The first order James test (order=1) is done using (James, 1951, p. 324):
#' \deqn{J_{crit} = \chi_{crit}^2\times\left(1 + \frac{3\times\chi_{crit}^2 + k + 1}{2\times\left(k^2 - 1\right)}\times\lambda\right)}
#' With:
#' \deqn{\lambda = \sum_{j=1}^k \frac{\left(1 - h_j\right)^2}{v_j}}
#' \deqn{\chi_{crit}^2 = Q\left(\chi^2\left(1-\alpha, df\right)\right)}
#' \deqn{v_j = n_j - 1}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\chi_{crit}^2} the critical chi-square value at alpha level
#' }
#' 
#' 
#' The second order James test (order=2) is done using (James)
#' \deqn{J_{crit}=C + \frac{1}{2}\times\left(3\times\chi_4+\chi_2\right)\times\lambda+\frac{1}{16}\times\left(3\times\chi_4+\chi_2\right)^2\times\left(1-\frac{k-3}{C})\times\lambda^2 + \frac{1}{2}\times\left(3\times\chi_4+\chi_2\right) \times\left(\left(8\times R_{23}-10\times R_{22}+4\times\ R_{21}-6\times R_{12}^2+8\times R_{12}\times R_{11}-4\times\ R_{11}^2\right)+\left(2\times R_{23}-4\times R_{22}+2\times R_{21}-2\times R_{12}^2+4\times R_{12}\times R_{11}-2\times R_{11}^2\right)\times\left(\chi_2-1\right)+\frac{1}{4}\times\left(-R_{12}^2+4\times R_{12}\times R_{11}-2\times R_{12}\times R_{10}-4\times R_{11}^2+4\times R_{11}\times R_{10}-R_{10}^2\right) \times\left(3\times\chi_4-2\times\chi_2-1)\right) +\right(R_{23}-3\times R_{22}+3\times R_{21}-R_{20}\right)\times\left(5\times\chi_6+2\times\chi_4+\chi_2\right) + \frac{3}{16}\times\left(R_{12}^2-4\times R_{23}+6\times R_{22}-4\times R_{21}+R_{20}\right)\times \left(35\times\chi_8+15\times\chi_6+9\times\chi_4+5\times\chi_2\right) + \frac{1}{16}\times\left(-2\times R_{22}^2+4\times R_{21}-R_{20}+2\times R_{12}\times R_{10}-4\times R_{11}\times R_{10}+R_{10}^2\right)\times\left(9\times\chi_8-3\times\chi_6-5\times\chi_4-\chi_2\right)+\frac{1}{4}\times\left(-R_{22}+R_{11}^2\right)\times\left(27\times\chi_8+3\times\chi_6+\chi_4+\chi_2\right)+\frac{1}{4}\times\left(R_{23}-R_{12}\times R_{11}\right)\times\left(45\times\chi_8+9\times\chi_6+7\times\chi_4+3\times\chi_2\right)}
#' With:
#' \deqn{\lambda_2 = \sum_{j=1}^k \frac{\left(1 - h_j\right)^2}{v_j^*}}
#' \deqn{v_j^* = n_j - 2}
#' \deqn{\chi_{2\times s} = \frac{\left(\chi_{crit}^2\right)^s}{\prod_{i=1}^s\left(k + 2\times i - 3\right)}}
#' \deqn{R_{xy} = \sum_{j=1}^k \frac{h_j^y}{\left(v_j^*\right)^x}}
#' 
#' If setting 'iters=TRUE' a binary search will be done estimating the p-value to find
#' \eqn{J_{crit} = J} by changing the \eqn{\alpha}
#' 
#' The use of \eqn{v_j^* = n_j - 2} for the James order 2 test is based on James (1951, p. 328) 
#' which can also be found in Deshon and Alexander (1994, p. 331).
#' 
#' However, others use \eqn{v_j^* = n_j - 1} for example Myers (1998, p. 209) and Cribbie et al. (2012, p. 62)
#' By setting 'alt=TRUE' this alternative version will be used.
#' 
#' @references 
#' Cribbie, R. A., Fiksenbaum, L., Keselman, H. J., & Wilcox, R. R. (2012). Effect of non-normality on test statistics for one-way independent groups designs: Effects of non-normality on test statistics. *British Journal of Mathematical and Statistical Psychology, 65*(1), 56–73. https://doi.org/10.1111/j.2044-8317.2011.02014.x
#' 
#' James, G. S. (1951). The comparison of several groups of observations when the ratios of the population variances are unknown. *Biometrika, 38*(3–4), 324–329. https://doi.org/10.1093/biomet/38.3-4.324
#' 
#' Myers, L. (1998). Comparability of the james’ second-order approximation test and the alexander and govern  A  statistic for non-normal heteroscedastic data. *Journal of Statistical Computation and Simulation, 60*(3), 207–222. https://doi.org/10.1080/00949659808811888
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' scores = c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, NA, 90, 25, 40, 70, 65, NA, 70, 98, 40, 65, 60, 35, NA, 50, 40, 75, NA, 65, 70, NA, 20, 80, 35, NA, 68, 70, 60, 70, NA, 80, 98, 10, 40, 63, 75, 80, 40, 90, 100, 33, 36, 65, 78, 50)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' ts_james(scores, groups, order=0)
#' ts_james(scores, groups, order=1, iters=TRUE)
#' ts_james(scores, groups, order=1, iters=FALSE)
#' ts_james(scores, groups, order=2, iters=TRUE)
#' ts_james(scores, groups, order=2, iters=FALSE)
#' ts_james(scores, groups, order=2, iters=TRUE, alt=TRUE)
#' ts_james(scores, groups, order=2, iters=FALSE, alt=TRUE)
#' 
#' @export
ts_james <- function(scores, groups, order=c(0, 1, 2), alpha=0.05, iters=FALSE, alt=FALSE){
  
  if (length(order)>1) {
    order=2
  }
  
  counts <- setNames(aggregate(scores~groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(scores~groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(scores~groups, FUN=var), c("category", "var"))
  res <- merge(counts, means, by = 'category')
  res <- merge(res, vars, by = 'category')
  
  #number of categories
  k <- dim(res)[1]
  
  #(first) degrees of freedom
  df1 = k - 1
  
  #weights (w), adjusted weights (h), weighted mean, and cochran test statistic
    res$w = res$n/res$var
    
    w = sum(res$w)
    res$h = res$w/w
    yw = sum(res$h*res$mean)
    chi2Cochran = sum(res$w*(res$mean - yw)**2)

  
  #lambda
    lamb = sum((1 - res$h)**2/(res$n - 1))
  
  #james test
    J = chi2Cochran
    
    if (order==0){
      chi2Stat = J
      pVal = pchisq(J, df1, lower.tail=FALSE) 
      reject = pVal < alpha
      comment = "for large category sizes"
      
      out <- data.frame(J, df1, pVal, reject)
      colnames(out)<-c("statistic", "df", "pValue", "reject H0")}
    
    else{
      cCrit = qchisq(1-alpha, df1)
      
      if (order==1){
        Jcrit = cCrit*(1 + (3*cCrit  + k + 1)/(2*(k**2 - 1))*lamb)
        
        if (iters){
          comment = "first-order with iterations for p-value approximation"
          pLow = 0
          pHigh = 1
          pVal = 0.05
          nIter = 1
          whileDo = TRUE
          
          while (whileDo){
            cCrit = qchisq(1-pVal, df1)
            Jcrit = cCrit*(1 + (3*cCrit  + k + 1)/(2*(k**2 - 1))*lamb)
            
            if (Jcrit < J){
              pHigh = pVal
              pVal = (pLow + pVal)/2}
            else if (Jcrit > J){
              pLow = pVal
              pVal = (pHigh + pVal)/2}
            
            nIter = nIter + 1
            
            if (Jcrit == J | nIter >= 800){
              whileDo = FALSE}
          }
          
          reject = pVal < alpha
          
          out <- data.frame(J, df1, pVal, reject)
          colnames(out)<-c("statistic", "df", "pValue", "reject H0")}
        
        else{
          comment = "first-order"
          reject = J > Jcrit
          
          out <- data.frame(J, df1, Jcrit, reject)
          colnames(out)<-c("statistic", "df", "J-critical", "reject H0")}
      }
      
      else{
        if (!alt){
          comment = "second order"
          res$v = res$n - 2
          lamb = sum((1 - res$h)**2/res$v)}
        
        else{
          comment = "second order with alternative v (v = n -1)"
          res$v = res$n - 1}
        
        R10 = sum(res$h**0 / res$v**1)
        R11 = sum(res$h**1 / res$v**1)
        R12 = sum(res$h**2 / res$v**1)
        R20 = sum(res$h**0 / res$v**2)
        R21 = sum(res$h**1 / res$v**2)
        R22 = sum(res$h**2 / res$v**2)
        R23 = sum(res$h**3 / res$v**2)
        
        c2 = cCrit**1/(k + 2*1 - 3)
        c4 = c2 * cCrit/(k + 2*2 - 3)
        c6 = c4 * cCrit/(k + 2*3 - 3)
        c8 = c6 * cCrit/(k + 2*4 - 3)
        
        Jcrit = cCrit + 1/2*(3*c4+c2)*lamb + 
          1/16*(3*c4 + c2 )**2*(1-(k-3)/cCrit)*lamb**2 + 
          1/2*(3*c4 + c2 )*
          ((8*R23 - 10*R22 + 4*R21 - 6*R12**2 + 8*R12*R11 - 4*R11**2) + 
             (2*R23 - 4*R22 + 2*R21 - 2*R12**2 + 4*R12*R11 - 2*R11**2)*(c2 - 1) + 
             1/4*(-R12**2 + 4*R12*R11 - 2*R12*R10 - 4*R11**2 + 4*R11*R10 - R10**2 )*(3*c4 - 2*c2 - 1)) + 
          (R23 - 3*R22 + 3*R21 - R20)*(5*c6 + 2*c4 + c2) + 
          3/16*(R12**2 - 4*R23 + 6*R22 - 4*R21 + R20)*(35*c8 + 15*c6 + 9*c4 + 5*c2) + 
          1/16*(-2*R22**2 + 4*R21 - R20 + 2*R12*R10 - 4*R11*R10 + R10**2)*(9*c8 - 3*c6 - 5*c4 - c2) + 
          1/4*(-R22 + R11**2 )*(27*c8 + 3*c6 + c4 + c2) + 
          1/4*(R23 - R12*R11)*(45*c8 + 9*c6 + 7*c4 + 3*c2)
        
        if (iters){
          comment = paste(comment, ", using iterations for p-value approximation")
          pLow = 0
          pHigh = 1
          pVal = 0.05
          nIter = 1
          whileDo = TRUE
          
          while (whileDo){
            cCrit = qchisq(1-pVal, df1)
            
            #(re)calculate chi values
            c2 = cCrit**1/(k + 2*1 - 3)
            c4 = c2 * cCrit/(k + 2*2 - 3)
            c6 = c4 * cCrit/(k + 2*3 - 3)
            c8 = c6 * cCrit/(k + 2*4 - 3)
            
            #calculate Jcrit
            Jcrit = cCrit + 1/2*(3*c4+c2)*lamb + 
              1/16*(3*c4 + c2 )**2*(1-(k-3)/cCrit)*lamb**2 + 
              1/2*(3*c4 + c2 )*
              ((8*R23 - 10*R22 + 4*R21 - 6*R12**2 + 8*R12*R11 - 4*R11**2) + 
                 (2*R23 - 4*R22 + 2*R21 - 2*R12**2 + 4*R12*R11 - 2*R11**2)*(c2 - 1) + 
                 1/4*(-R12**2 + 4*R12*R11 - 2*R12*R10 - 4*R11**2 + 4*R11*R10 - R10**2 )*(3*c4 - 2*c2 - 1)) + 
              (R23 - 3*R22 + 3*R21 - R20)*(5*c6 + 2*c4 + c2) + 
              3/16*(R12**2 - 4*R23 + 6*R22 - 4*R21 + R20)*(35*c8 + 15*c6 + 9*c4 + 5*c2) + 
              1/16*(-2*R22**2 + 4*R21 - R20 + 2*R12*R10 - 4*R11*R10 + R10**2)*(9*c8 - 3*c6 - 5*c4 - c2) + 
              1/4*(-R22 + R11**2 )*(27*c8 + 3*c6 + c4 + c2) + 
              1/4*(R23 - R12*R11)*(45*c8 + 9*c6 + 7*c4 + 3*c2)
            
            if (Jcrit < J){
              pHigh = pVal
              pVal = (pLow + pVal)/2}
            else if (Jcrit > J){
              pLow = pVal
              pVal = (pHigh + pVal)/2}
            
            nIter = nIter + 1
            
            if (Jcrit == J | nIter >= 500){
              whileDo = FALSE}
          }
          
          reject = pVal < alpha
          
          out <- data.frame(J, df1, pVal, reject)
          colnames(out)<-c("statistic", "df", "pValue", "reject H0")}
        
        else{
          reject = J > Jcrit
          
          out <- data.frame(J, df1, Jcrit, reject)
          colnames(out)<-c("statistic", "df", "J-critical", "reject H0")}
      }
    }
  return (out)
}