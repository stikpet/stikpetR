#' Özdemir-Kurt Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
#' 
#' There are quite some alternatives for this, the stikpet library has Fisher, Welch, James, Box, Scott-Smith, Brown-Forsythe, Alexander-Govern, Mehrotra modified Brown-Forsythe, Hartung-Agac-Makabi, Özdemir-Kurt and Wilcox as options. See the notes from ts_fisher_owa() for some discussion on the differences.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{the test statistic (B2 value)}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Özdemir & Kurt, 2006, pp. 85-86):
#' \deqn{ B^2 = \sum_{j=1}^k\left(c_j\times\sqrt{\ln\left(1+\frac{t_j^2}{v_i}\right)}\right)^2 }
#' \deqn{ df = k - 1 }
#' \deqn{sig. = 1 - \chi^2\left(B^2, df\right)}
#' With:
#' \deqn{\chi_{crit}^2 = Q\left(chi_{crit}^2\left(1 - \alpha, df\right)\right)}
#' \deqn{t_j = \frac{\bar{x}_j - \bar{x}_w}{\sqrt{\frac{s_j^2}{n_j}}} }
#' \deqn{c_j = \frac{4\times v_j^2 + \frac{5\times\left(2\times z_{crit}^2+3\right)}{24}}{4\times v_j^2+v_j+\frac{4\times z_{crit}^2+9}{12}}\times\sqrt{v_j} }	
#' \deqn{v_j = n_j - 1}
#' \deqn{\bar{x}_w = \sum_{j=1}^k h_j\times \bar{x}_j}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{z_{crit} = Q\left(\Phi\left(1 - \frac{\alpha}{2}\right)\right)}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{w_j^{\ast}} the modified weight for category j
#' \item \eqn{h_j^{\ast}} the adjusted modified weight for category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\alpha} the significance level (usually 0.05)
#' \item \eqn{Q\left(\dots\right)} the quantile (inverse) distribution function
#' \item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' \item \eqn{\chi^2\left(\dots\right)} the cumulative density function of the chi-square distribution
#' }
#' 
#' A binary search for a p-value is done such that \eqn{B^2 = \chi_{crit}^2}.
#' 
#' @references 
#' Özdemir, A. F., & Kurt, S. (2006). One way fixed effect analysis of variance under variance heterogeneity and a solution proposal. *Selçuk Journal of Applied Mathematics, 7*(2), 81–90.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_ozdemir_kurt_owa <- function(nomField, scaleField, categories=NULL){
  
  datFrame = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$nomField %in% categories),]}
  
  counts <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  myRes$se <- sqrt(myRes$var/myRes$n)
  myRes$w <- (1/myRes$se^2)/sum(1/myRes$se^2)
  wMean <- sum(myRes$w * myRes$mean)
  myRes$t <- (myRes$mean - wMean)/myRes$se
  myRes$v <- myRes$n - 1
  k <- dim(myRes)[1]
  df <- k - 1
  
  
  pLow <- 0
  pHigh <- 1
  pVal <- 0.05
  nIter = 1
  repeat {
    zCrit <- qnorm(1-pVal/2)
    myRes$c <- (4*myRes$v^2 + 5*(2*zCrit^2 + 3)/24)/(4*myRes$v^2 + myRes$v + (4*zCrit^2 + 9)/12) * sqrt(myRes$v)
    B2 <- sum((myRes$c*sqrt(log(1 + myRes$t^2/myRes$v)))^2)
    chiCrit <- qchisq(1-pVal, df, lower.tail=TRUE)
    if (chiCrit == B2 || nIter > 100) {
      break
    }
    
    if (chiCrit < B2) {
      pHigh <- pVal
      pVal <- (pLow + pVal)/2
    } 
    else if (chiCrit > B2) {
      pLow <- pVal
      pVal <- (pHigh + pVal)/2
    }
    
    nIter = nIter + 1
    
  }
  
  n <- sum(myRes$n)
  testResults <- data.frame(n, B2, df, pVal)
  colnames(testResults)<-c("n", "statistic", "df", "p-value")
  
  return (testResults)
  
}