#' Box F-Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' Box proposed a correction to the original Fisher one-way ANOVA, on both the test-statistic and the degrees of freedom.
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
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Box, 1954, p. 299):
#' \deqn{F_{B} = \frac{F_F}{c}}
#' \deqn{df_1 = \frac{\left(\sum_{j=1}^k\left(n - n_j\right)\times s_j^2\right)^2}{\left(\sum_{j=1} n_j\times s_j^2\right)^2 + n\times\sum_{j=1}^k\left(n - 2\times n_j\right)\times s_j^4}}
#' \deqn{df_2 =\frac{\left(\sum_{j=1}^k\left(n_j-1\right)\times s_j^2\right)^2}{\sum_{j=1}^k\left(n_j-1\right)\times s_j^4}} 
#' \deqn{sig. = 1 - F\left(F_B, df_1, df_2\right)}
#' 
#' With:
#' \deqn{c = \frac{n - k}{n\times\left(k - 1\right)}\times \frac{\sum_{j=1}^k\left(n-n_j\right)\times s_j^2}{\sum_{j=1}^k\left(n_j-1\right)\times s_j^2}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{F-F} the F statistic of the classic/Fisher one-way ANOVA. See *ts_fisher_owa()* for details.
#' \item \eqn{n} the total sample size
#' \item \eqn{k} the number of categories
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{\bar{x}_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{F\left(\dots,\dots,\dots\right)} the cumulative distribution function of the F distribution.
#' }
#' 
#' This also appears to give the same results for the test statistic, df_1 as the Brown-Forsythe test for means but a different df_2
#' 
#' The *doex* and *onewaytests* libraries used to have a different method for calculating df_1 but after personal communication with the creators of those packages they mentioned to fix it in an update.
#' 
#' Asiribo and Gurland (1990) derive the same correction as Box, although their notation for \eqn{df_1^*} is different, but will give the same result. 
#' 
#' @references 
#' Asiribo, O., & Gurland, J. (1990). Coping with variance heterogeneity. *Communications in Statistics - Theory and Methods, 19*(11), 4029-4048. doi:10.1080/03610929008830427
#' 
#' Box, G. E. P. (1954). Some theorems on quadratic forms applied in the study of analysis of variance problems, I: Effect of inequality of variance in the one-way classification. *The Annals of Mathematical Statistics, 25*(2), 290-302. doi:10.1214/aoms/1177728786
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_box_owa <- function(nomField, scaleField, categories=NULL){
  
  datFrame = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$nomField %in% categories),]}
  
  counts <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  k <- dim(myRes)[1]
  n = sum(myRes$n)
  xbar = mean(datFrame$scaleField)
  
  b = (n - k)/(n*(k - 1)) * sum((n - myRes$n)*myRes$var)/sum((myRes$n - 1)*myRes$var)
  Fstat = ((n - k)/(k - 1))*((sum(myRes$n*(myRes$mean-xbar)^2))/(sum((myRes$n-1)*myRes$var)))
  Fadj = Fstat / b
  
  df1 = sum((n - myRes$n)*myRes$var)^2/(sum(myRes$n*myRes$var)^2 + n*sum((n - 2*myRes$n)*myRes$var^2))
  
  df2 = sum((myRes$n - 1)*myRes$var)^2/sum((myRes$n-1)*myRes$var^2)
  
  pVal <- pf(Fadj, df1, df2, lower.tail = FALSE)
  
  testResults <- data.frame(n, k, Fadj, df1, df2, pVal)
  colnames(testResults)<-c("n", "k", "statistic", "df1", "df2", "pValue")
  
  return (testResults)
  
}



