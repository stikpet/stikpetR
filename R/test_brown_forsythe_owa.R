#' Brown-Forsythe Means Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
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
#' \item{statistic}{the test statistic (F value)}
#' \item{df1}{the degrees of freedom 1}
#' \item{df2}{the degrees of freedom 2}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Brown & Forsythe, 1974, p. 130):
#' \deqn{F_{BF} = \frac{\sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}{\sum_{j=1}^k\left(1 - \frac{n_j}{n}\right)\times s_j^2}}
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 =\frac{\left(\sum_{j=1}^k\left(1 - \frac{n_j}{n}\right)\times s_j^2\right)^2}{\sum_{j=1}^k\frac{\left(1 - \frac{n_j}{n}\right)\times s_j^4}{n_j - 1}}} 
#' \deqn{sig. = 1 - F\left(F_{BF}, df_1, df_2\right)}
#' With:
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^k n_j\times\bar{x}_j}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{F\left(\dots,\dots,\dots\right)} the cumulative distribution function of the F distribution.
#' }
#' 
#' This appears to give the same results as the Box correction, except for \eqn{df_2}. 
#' 
#' @references 
#' Brown, M. B., & Forsythe, A. B. (1974). The small sample behavior of some statistics which test the equality of several means. *Technometrics, 16*(1), 129-132. https://doi.org/10.1080/00401706.1974.10489158
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_brown_forsythe_owa <- function(nomField, scaleField, categories=NULL){
  
  datFrame = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$nomField %in% categories),]}
  
  counts <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  n <- sum(myRes$n)
  xbar <- sum(myRes$n*myRes$mean)/n
  Fstat <- sum(myRes$n*(myRes$mean - xbar)^2) / sum((1 - myRes$n/n)*myRes$var) 
  k <- dim(myRes)[1]
  df1 <- k - 1
  df2 <- 1/(sum((((1 - myRes$n/n)*myRes$var/sum((1 - myRes$n/n)*myRes$var))^2)/(myRes$n - 1)))
  pVal <- pf(Fstat, df1, df2, lower.tail = FALSE)
  
  testResults <- data.frame(n, k, Fstat, df1, df2, pVal)
  colnames(testResults)<-c("n", "k", "statistic", "df1", "df2", "p-value")
  
  return (testResults)
  
}



