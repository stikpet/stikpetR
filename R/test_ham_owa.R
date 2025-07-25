#' Hartung-Argac-Makambi Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
#' 
#' This test is a modification of the Welch one-way ANOVA.
#' 
#' There are quite some alternatives for this, the stikpet library has Fisher, Welch, James, Box, Scott-Smith, Brown-Forsythe, Alexander-Govern, Mehrotra modified Brown-Forsythe, Hartung-Agac-Makabi, Ozdemir-Kurt and Wilcox as options. See the notes from ts_fisher_owa() for some discussion on the differences.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' @param version the phi method calculation to use (see details)
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
#' The formula used (Hartung et al., 2002, p. 206):
#' \deqn{W = \frac{\frac{1}{k-1}\times\sum_{j=1}^k w_j^{\ast}\times\left(\bar{x}_j - \bar{y}_w^{\ast}\right)^2}{1 + \frac{2\times\left(k-2\right)}{k^2-1}\times \lambda^*}}
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 = \frac{k^2-1}{3\times\lambda^{\ast}}}
#' \deqn{sig. = 1 - F\left(W, df_1, df_2\right)}
#' With:
#' \deqn{\bar{y}_w^{\ast} = \sum_{j=1}^k h_j^{\ast}\times \bar{x}_j}
#' \deqn{h_j^{\ast} = \frac{w_j^{\ast}}{w^{\ast}}}
#' \deqn{w_j^{\ast} = \frac{n_j}{s_j^2}\times\frac{1}{\phi_j}}
#' \deqn{w^{\ast} = \sum_{j=1}^k w_j^{\ast}}
#' \deqn{\phi_j = \frac{n_j + 2}{n_j + 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{j=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\lambda^{\ast} = \sum_{j=1}^k \frac{\left(1 - h_j^{\ast}\right)^2}{n_j - 1}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{x_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{w_j^{\ast}} the modified weight for category j
#' \item \eqn{h_j^{\ast}} the adjusted modified weight for category j
#' \item \eqn{df_i} the i-th degrees of freedom
#' }
#' 
#' Note that the numerator in \eqn{W} is the same as the Cochran test statistic.
#' 
#' Cavis and Yazici (2020, p. 6) uses \eqn{\phi_j=\frac{n_j-1}{n_j-3}}. However the original article states that these are unbalanced weights of the Welch test and in their experience, using these makes the test too conservative. In the original article they find from their simulation experience 
#' that using (n_j+2)/(n_j+1) gives reliable results for small sample sizes, and a large number of populations (Hartung et al. p. 207).
#' 
#' By setting 'version=2' the same version for \eqn{\phi} as in the Doex library will be used
#' 
#' @references 
#' Cavus, M., & Yazici, B. (2020). Testing the equality of normal distributed and independent groups' means under unequal variances by doex package. *The R Journal, 12*(2), 134. https://doi.org/10.32614/RJ-2021-008
#' 
#' Hartung, J., Argac, D., & Makambi, K. H. (2002). Small sample properties of tests on homogeneity in one-way anova and meta-analysis. *Statistical Papers, 43*(2), 197-235. https://doi.org/10.1007/s00362-002-0097-8
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_ham_owa <- function(nomField, scaleField, categories=NULL, version=c(1, 2)){
  
  if (length(version)>1) {
    version=1
  }
  
  datFrame = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$nomField %in% categories),]}
  
  counts <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scaleField~datFrame$nomField, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  if (version==1) {
    myRes$phi <- (myRes$n + 2)/(myRes$n + 1)
  } else{
    myRes$phi <- (myRes$n - 1)/(myRes$n - 3) 
  }
  myRes$w <- myRes$n/(myRes$phi*myRes$var)
  wSum <- sum(myRes$w)
  myRes$h <- myRes$w/wSum
  wMean <- sum(myRes$h * myRes$mean)
  k <- dim(myRes)[1]
  W <- sum(myRes$w*(myRes$mean - wMean)^2) / ((k - 1) + 2*(k - 2)*1/(k + 1)*sum((1 - myRes$h)^2/(myRes$n - 1)))
  df1 <- k - 1
  df2 <- (k^2 - 1)/(3*sum((1 - myRes$h)^2/(myRes$n - 1)))
  pVal <- pf(W, df1, df2, lower.tail = FALSE)
  
  n <- sum(myRes$n)
  testResults <- data.frame(n, k, W, df1, df2, pVal)
  colnames(testResults)<-c("n", "k", "statistic", "df1", "df2", "p-value")
  
  return (testResults)
  
}



