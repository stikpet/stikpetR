#' Alexander-Govern Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
#' 
#' Schneider and Penfield (1997) looked at the Welch, Alexander-Govern and the James test (they ignored the Brown-Forsythe since they found it to perform worse than Welch or James), and concluded: "Under variance heterogeneity, Alexander-Govern's approximation was not only comparable to the Welch test and the James second-order test but was superior, in certain instances, when coupled with the power results for those tests" (p. 285).
#' 
#' There are quite some alternatives for this, the stikpet library has Fisher, Welch, James, Box, Scott-Smith, Brown-Forsythe, Alexander-Govern, Mehrotra modified Brown-Forsythe, Hartung-Agac-Makabi, Ozdemir-Kurt and Wilcox as options. See the notes from ts_fisher_owa() for some discussion on the differences.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{k}{the number of categories}
#' \item{statistic}{the test statistic (chi-square value)}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Alexander & Govern, 1994, pp. 92-94):
#' \deqn{A = \sum_{j=1}^k z_j^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(A, df\right)}
#' 
#' With:
#' \deqn{z_j = c_j + \frac{c_j^3 + 3\times c_j}{b_j} - \frac{4\times c_j^7 + 33\times c_j^5 + 240\times c_j^3 + 855\times c_j}{10\times b_j^2 + 8\times b_j\times c_j^4 + 1000\times b_j}}
#' \deqn{c_j = \sqrt{a_j\times\ln\left(1 + \frac{t_j^2}{n_j - 1}\right)}}
#' \deqn{b_j = 48\times a_j^2}
#' \deqn{a_j = n_j - 1.5}
#' \deqn{t_j = \frac{\bar{x}_j - \bar{y}_w}{\sqrt{\frac{s_j^2}{n_j}}}}
#' \deqn{\bar{y}_w = \sum_{j=1}^k h_j\times \bar{x}_j}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{n} the total sample size
#' \item \eqn{k} the number of categories
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\chi^2\left(\dots, \dots\right)} the cumulative distribution function of the chi-square distribution.
#' }
#' 
#' @references 
#' Alexander, R. A., & Govern, D. M. (1994). A new and simpler approximation for ANOVA under variance heterogeneity. *Journal of Educational Statistics, 19*(2), 91-101. doi:10.2307/1165140
#' 
#' Schneider, P. J., & Penfield, D. A. (1997). Alexander and Govern's approximation: Providing an alternative to ANOVA under variance heterogeneity. *The Journal of Experimental Education, 65*(3), 271-286. doi:10.1080/00220973.1997.9943459
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_alexander_govern_owa <- function(nomField, scaleField, categories=NULL){
  
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
  myRes$recSe <- 1/(myRes$se^2)
  myRes$w <- myRes$recSe/sum(myRes$recSe)
  
  muHat <- sum(myRes$w*myRes$mean)
  t <- (myRes$mean - muHat)/myRes$se
  v <- myRes$n - 1
  a <- v - 0.5
  b <- 48*a^2
  c <- sqrt(a*log(1 + t^2/v))
  z = c + (c^3 + 3*c)/b - (4*c^7 + 33*c^5 + 240*c^3 + 855*c)/(10*b^2 + 8*b*c^4 + 1000*b)
  A <- sum(z^2)
  k <- dim(myRes)[1]
  df <- k - 1
  pVal <- pchisq(A, df, lower.tail=FALSE)
  
  n = sum(myRes$n)   
  testResults <- data.frame(n, A, df, pVal)
  colnames(testResults)<-c("n", "statistic", "df", "p-value")
  
  return (testResults)
  
}



