#' Wilcox Test
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
#' \item{statistic}{the test statistic (chi-square value value)}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' The formula used (Wilcox, 1988, pp. 110-111)
#' \deqn{H = \frac{\sum_{j=1}^k \left(W_j - \bar{W}\right)^2}{\hat{\theta}}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(H, df\right)}
#' With:
#' \deqn{W_j = b_j\times x_{n_j,j} + \frac{1 - b_j}{n_j}\times\sum_{i=1}^{n_j-1} x_{i,j}}
#' \deqn{\bar{W} = \frac{\sum_{j=1}^k W_j}{k}}
#' \deqn{b_j = \frac{1 + \sqrt{\frac{\left(n_j - 1\right)\times\left(n_j\times\hat{\theta}\right)}{s_j^2}}}{n_j}}
#' \deqn{\hat{\theta} = \text{max}\left\{\frac{s_1^2}{n_1}, \frac{s_2^2}{n_2}, \dots, \frac{s_k^2}{n_k}\right\}}
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
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\chi^2\left(\dots\right)} the cumulative density function of the chi-square distribution
#' }
#' 
#' The original article has an error in the formula for \eqn{b_j}. There are missing brackets. 
#' Using the population version in the article of \eqn{c_j} the formula used here was adapted.
#' 
#' @references 
#' Wilcox, R. R. (1988). A new alternative to the ANOVA F and new results on James’s second-order method. *British Journal of Mathematical and Statistical Psychology, 41*(1), 109–117. https://doi.org/10.1111/j.2044-8317.1988.tb00890.x
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_wilcox_owa <- function(nomField, scaleField, categories=NULL){
  
  datFrame = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$nomField %in% categories),]}
  
  counts <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=var), c("category", "var"))
  maxs <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=max), c("category", "max"))
  sums <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=sum), c("category", "sum"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  myRes <- merge(myRes, maxs, by = 'category')  
  myRes <- merge(myRes, sums, by = 'category')  
  
  k <- dim(myRes)[1]
  
  d = max(myRes$var/myRes$n)
  bj = (1 + sqrt((myRes$n - 1)*(myRes$n*d - myRes$var)/myRes$var))/myRes$n
  
  Wj = bj*myRes$max + (1-bj)/myRes$n * (myRes$sum - myRes$max)  
  Wm = sum(Wj)/k
  
  H = sum((Wj - Wm)**2)/d
  df = k - 1
  pValue = 1 - pchisq(H, df)
  
  n <- sum(myRes$n)  
  testResults <- data.frame(n, H, df, pValue)
  colnames(testResults)<-c("n", "statistic", "df", "p-value")
  
  return (testResults)
  
}