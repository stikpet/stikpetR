#' Hartung-Argaç-Makambi Test
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @param version the phi method calculation to use (see details)
#' @returns 
#' A dataframe with:
#' \item{statistic}{the F-statistic from the test}
#' \item{df1}{the first degrees of freedom}
#' \item{df2}{the second degrees of freedom}
#' \item{pValue}{the significance (p-value)}
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
#' Cavis and Yazici (2020, p. 6) uses \eqn{\phi_j=\frac{n_j-1}{n_j-3}}. 
#' However the original article states that these are unbalanced weights of the Welch test 
#' and in their experience, using these makes the test too conservative. 
#' In the original article they find from their simulation experience 
#' that using (n_j+2)/(n_j+1) gives reliable results for small sample sizes, 
#' and a large number of populations (Hartung et al. p. 207).
#' 
#' By setting 'version=2' the same version for \eqn{\phi} as in the Doex library will be used
#' 
#' @references 
#' Cavus, M., & Yazıcı, B. (2020). Testing the equality of normal distributed and independent groups’ means under unequal variances by doex package. *The R Journal, 12*(2), 134. https://doi.org/10.32614/RJ-2021-008
#' 
#' Hartung, J., Argaç, D., & Makambi, K. H. (2002). Small sample properties of tests on homogeneity in one-way anova and meta-analysis. *Statistical Papers, 43*(2), 197–235. https://doi.org/10.1007/s00362-002-0097-8
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
#' ts_hartung_agac_makabi(scores, groups)
#' ts_hartung_agac_makabi(scores, groups, version=2)
#' 
#' @export 
ts_hartung_agac_makabi <- function(scores, groups, version=c(1, 2)){
  
  if (length(version)>1) {
    version=1
  }
  
  datFrame = na.omit(data.frame(groups, scores))
  
  counts <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(datFrame$scores~datFrame$groups, FUN=var), c("category", "var"))
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
  
  testResults <- data.frame(W, df1, df2, pVal)
  colnames(testResults)<-c("statistic", "df1", "df2", "pValue")
  
  return (testResults)
  
}
