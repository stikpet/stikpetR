#' Spearman Rho / Rank Correlation Coefficient
#' 
#' @description 
#' The Spearman Rank Correlation Coefficient is the Pearson Correlation Coefficient, 
#' after the scores first have been converted to ranks.
#' 
#' This function makes use of *di_spearman()* for the test of this correlation,
#' which requires the *pspearman* library for exact computations.
#' 
#' @param ordField1 the numeric scores of the first variable
#' @param ordField2 the numeric scores of the second variable
#' @param levels1 vector, optional. the categories to use from ordField1
#' @param levels2 vector, optional. the categories to use from ordField2
#' @param test the test to be used. Either "t" (default), "as89", "exact", "iman-conover", "z-fieller", "z-olds", "none"
#' @param cc boolean to indicate the use of a continuity correction
#' 
#' @returns 
#' A dataframe with:
#' \item{rs}{the correlation coefficient}
#' \item{pValue}{the significance (p-value)}
#' \item{statistic}{the statistic from the test (only if applicable)}
#' \item{df}{the degrees of freedom (only if applicable)}
#' 
#' 
#' @details 
#' The formula used is (Spearman, 1904, p. 77):
#' \deqn{r_s = \frac{SS_{r_x, r_y}}{SS_{r_x}\times SS_{r_y}}}
#' With:
#' \deqn{SS_{r_x} = \sum_{i=1}^n \left(r_{x_i} - \bar{r}_x\right)^2}
#' \deqn{SS_{r_y} = \sum_{i=1}^n \left(r_{y_i} - \bar{r}_y\right)^2}
#' \deqn{SS_{r_x, r_y} = \sum_{i=1}^n \left(r_{x_i} - \bar{r}_x\right) \times \left(r_{y_i} - \bar{r}_y\right)}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{r_{x_i}} the i-th rank of the scores of the first variable
#' \item \eqn{r_{y_i}} the i-th rank of the scores of the second variable
#' \item \eqn{n} the total sample size (number of ranks)
#' }
#' 
#' If all the ranks are distinct (i.e. no ties) the formula can also be written as:
#' 
#' \deqn{r_s = 1 - \frac{6}{n\times\left(n^2 - 1\right)}\times S}
#' With:
#' \deqn{S = \sum_{i=1}^n d_i^2}
#' \deqn{d_i^2 = \left(r_{x_i} - r_{y_i}\right)^2}
#' 
#' The test can be performed in different ways. Options to choose from are:
#' \itemize{
#' \item "t" uses a Student t distribution approximation
#' \item "z-fieller" uses a standard normal approximation from Fieller
#' \item "z-old" uses standard normal approximation from Old
#' \item "iman-conover" a combination of z and t distribution from Iman and Conover
#' \item "AS89" uses the AS 89 algorithm
#' \item "exact" uses an exact distribution
#' }
#' 
#' See for the details of each the *di_spearman()* function
#' 
#' A continuity correction can be applied (Zar, 1972, p. 579):
#' \deqn{r_s^{cc} = \left|r_s\right| - \frac{6}{n^3 - n}}
#' 
#' **Alternatives**
#' 
#' *R's stats*
#' 
#' Using the t-approximation:
#' 
#' cor.test(ord1, ord2, method="spearman")
#' 
#' Using AS89
#' 
#' cor.test(ord1, ord2, method="spearman", exact=TRUE)
#' 
#' *library(pspearman)*
#' 
#' spearman.test(ord1, ord2, approximation="t-distribution")
#' 
#' spearman.test(ord1, ord2, approximation="AS89")
#' 
#' spearman.test(ord1, ord2, approximation="exact")
#' 
#' @references 
#' Göktaş, A., & İşçi, Ö. (2011). A comparison of the most commonly used measures of association for doubly ordered square contingency tables via simulation. *Advances in Methodology and Statistics, 8*(1). doi:10.51936/milh5641
#' 
#' Spearman, C. (1904). The proof and measurement of association between two things. *The American Journal of Psychology, 15*(1), 72–101.
#' 
#' Zar, J. H. (1972). Significance testing of the Spearman rank correlation coefficient. *Journal of the American Statistical Association, 67*(339), 578–580. doi:10.1080/01621459.1972.10481251
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
r_spearman_rho <- function(ordField1, ordField2, levels1=NULL, levels2=NULL, 
                           test=c("t", "z-fieller", "z-olds", "iman-conover", "as89", "exact"), 
                           cc=FALSE){
  
  if (length(test)>1) {
    test = "t"
  }
  
  if (!is.null(levels1)){
    myFieldOrd = factor(ordField1, ordered = TRUE, levels = levels1)
    ordField1 = as.numeric(myFieldOrd)
  }
  
  if (!is.null(levels2)){
    myFieldOrd = factor(ordField2, ordered = TRUE, levels = levels2)
    ordField2 = as.numeric(myFieldOrd)
  }  
  
  dFr = na.omit(data.frame(ordField1, ordField2))
  
  rx = rank(dFr$ordField1)
  ry = rank(dFr$ordField2)
  
  mrx = mean(rx)
  mry = mean(ry)
  
  n = length(rx)
  
  SSrx = var(rx)*(n - 1)
  SSry = var(ry)*(n - 1)
  
  SSxy = sum((rx - mrx)*(ry - mry))
  
  rs = SSxy / sqrt(SSrx*SSry)  
  
  if (cc && test!="exact") {
    #(Zar, 1972, p. 579)
    rs = abs(rs) - 6/(n**3 - n)
  }
  
  if (test=="none"){
    results = rs}
  else{
    df = n - 2
    
    if (test=="as89") {
      distRes = di_spearman(n, rs, method="AS89")
    }
    else{
      distRes = di_spearman(n, rs, method=test)
    }
    
    pValue = distRes$pValue
    
    if (test=="t" || test=="iman-conover" || test=="as89") {
      statistic = distRes$statistic
      results = data.frame(rs, pValue, statistic, df)
    }
    else if (test=="z-fieller" || test=="z-olds"){
      statistic = distRes$statistic
      results = data.frame(rs, pValue, statistic)
    }
    
    else if (test=="exact") {
      results = data.frame(rs, pValue)
    }
  }
  return(results)
}