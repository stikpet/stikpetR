#' Cochran One-Way ANOVA
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' Note that according to Hartung et al. (2002, p. 225) the Cochran test is the standard test in meta-analysis, but should not be used, since it is always too liberal.
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
#' \item{statistic}{the chi-square-statistic from the test}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used is (Cavus & Yazıcı, 2020, p. 5; Hartung et al., 2002, p. 202; Mezui-Mbeng, 2015, p. 787):
#' \deqn{\chi_C^2 = \sum_{j=1}^k w_j\times\left(\bar{x}_j - \bar{y}_w\right)^2}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_C^2, df\right)}
#' 
#' With:
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
#' \item \eqn{df} the degrees of freedom
#' \item \eqn{\chi^2\left(\dots,\dots\right)} the cumulative distribution function of the chi-square distribution.
#' }
#' 
#' Couldn’t really find the formula in the original article which is from Cochran (1937)
#' 
#' @references 
#' Cavus, M., & Yazıcı, B. (2020). Testing the equality of normal distributed and independent groups’ means under unequal variances by doex package. *The R Journal, 12*(2), 134. https://doi.org/10.32614/RJ-2021-008
#' 
#' Cochran, W. G. (1937). Problems arising in the analysis of a series of similar experiments. *Supplement to the Journal of the Royal Statistical Society, 4*(1), 102–118. https://doi.org/10.2307/2984123
#' 
#' Hartung, J., Argaç, D., & Makambi, K. H. (2002). Small sample properties of tests on homogeneity in one-way anova and meta-analysis. *Statistical Papers, 43*(2), 197–235. https://doi.org/10.1007/s00362-002-0097-8
#' 
#' Mezui-Mbeng, P. (2015). A note on Cochran test for homogeneity in two ways ANOVA and meta-analysis. *Open Journal of Statistics, 5*(7), 787–796. https://doi.org/10.4236/ojs.2015.57078
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_cochran_owa <- function(nomField, scaleField, categories=NULL){
  
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
  
  chi2Coch = sum(wj*(xBars - yw)**2)
  
  df = k - 1
  
  pValue = 1 - pchisq(chi2Coch, df)
  
  statistic = chi2Coch
  res = data.frame(n, statistic, df, pValue)
  colnames(res) = c("n", "statistic", "df", "p-value")  
  
  return(res)
  
}