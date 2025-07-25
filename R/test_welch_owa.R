#' Welch One-Way ANOVA
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
#' 
#' Delacre et al. (2019) recommend to use the Welch ANOVA instead of the classic and Brown-Forsythe versions, but there are quite some alternatives for this, the stikpet library has Fisher, Welch, James, Box, Scott-Smith, Brown-Forsythe, Alexander-Govern, Mehrotra modified Brown-Forsythe, Hartung-Agac-Makabi, Ozdemir-Kurt and Wilcox as options. See the notes from ts_fisher_owa() for some discussion on the differences.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' 
#' @returns 
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{the test statistic (F value)}
#' \item{df1}{the degrees of freedom}
#' \item{df2}{the degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used is (Welch, 1951, pp. 334-335):
#' \deqn{F_w = \frac{\frac{1}{k - 1} \times \sum_{j=1}^k w_j\times\left(\bar{x}_j - \bar{y}_w\right)^2}{1 + 2\times\lambda\times\frac{k-2}{k^2 - 1}}}
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 = \frac{k^2 - 1}{3\times\lambda}}
#' \deqn{sig. = 1 - F\left(F_W, df_1, df_2\right)}
#' 
#' With:
#' \deqn{\lambda = \sum_{j=1}^k \frac{\left(1 - h_j\right)^2}{n_j - 1}}
#' \deqn{\bar{y}_w = \frac{\sum_{j=1}^k w_j\times\bar{x}_j}{\sum_{j=1}^k w_j} = \sum_{j=1}^k h_j\times\bar{x}_j}
#' \deqn{h_j = \frac{w_j}{w}}
#' \deqn{w = \sum_{j=1}^k w_j}
#' \deqn{w_j = \frac{n_j}{s_j^2}}
#' \deqn{s_j^2 = \frac{\sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}{n_j - 1}}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols:*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{n_j} the sample size of category j
#' \item \eqn{x_j} the sample mean of category j
#' \item \eqn{s_j^2} the sample variance of the scores in category j
#' \item \eqn{w_j} the weight for category j
#' \item \eqn{h_j} the adjusted weight for category j
#' \item \eqn{df_i} the i-th degrees of freedom
#' }
#' 
#' The formula can also be written as:
#' \deqn{F_W = \frac{\chi_{Cochran}^2}{k - 1 + 2\times\lambda\times\frac{k - 2}{k + 1}}}
#' Where \eqn{\chi_{Cochran}^2} is the test statistic of the Cochran one-way test
#' 
#' Cavus and Yazici (2020) make a difference between the Welch and the Welch-Aspin ANOVA. The only difference in the article is that with the Welch 2x(k-2) is used, 
#' while in the Welch-Aspin version 2xk-2. I think this is a mistake in their formula, since the article they refer to from Aspin is about two means.
#' 
#' Johansen F test (Johansen, 1980) will give the same results
#' 
#' @references 
#' Cavus, M., & Yazici, B. (2020). Testing the equality of normal distributed and independent groups' means under unequal variances by doex package. *The R Journal, 12*(2), 134. doi:10.32614/RJ-2021-008
#' 
#' Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology, 32*(1), 1-12. doi:10.5334/irsp.198
#' 
#' Johansen, S. (1980). The Welch-James approximation to the distribution of the residual sum of squares in a weighted linear regression. *Biometrika, 67*(1), 85-92. doi:10.1093/biomet/67.1.85
#' 
#' Welch, B. L. (1947). The generalization of `Student's' problem when several different population variances are involved. *Biometrika, 34*(1/2), 28-35. doi:10.2307/2332510
#' 
#' Welch, B. L. (1951). On the comparison of several mean values: An alternative approach. *Biometrika, 38*(3/4), 330-336. doi:10.2307/2332579
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_welch_owa <- function(nomField, scaleField, categories=NULL){
  
  dfr = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$nomField %in% categories),]}
  
  #overall mean and count
  xBar = mean(dfr$scaleField)
  n = nrow(dfr)
  
  #means and counts and variances per category
  xBars = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=mean)$x
  ns = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=length)$x
  vars = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=var)$x
  
  #number of categories
  k = length(ns)
  
  wj = ns/vars
  w = sum(wj)
  
  hj = wj/w
  yw = sum(hj*xBars)
  
  lambda = sum((1 - hj)**2/(ns - 1))
  
  Fvalue = (sum(wj*(xBars - yw)**2)/(k - 1)) / (1 + 2*(k - 2)/(k**2 - 1)*lambda)
  
  df1 = k - 1
  df2 = (k**2 - 1)/(3*lambda)
  
  pValue = 1 - pf(Fvalue, df1, df2)
  
  statistic = Fvalue
  res = data.frame(n, statistic, df1, df2, pValue)
  colnames(res) = c("n", "statistic", "df", "p-value")  
  
  
  return(res)
  
}



