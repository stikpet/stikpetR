#' Cochran One-Way ANOVA
#' 
#' @param scores the numeric scores variable
#' @param groups the groups variable
#' @returns 
#' A dataframe with:
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
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' scores = c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, NA, 90, 25, 40, 70, 65, NA, 70, 98, 40, 65, 60, 35, NA, 50, 40, 75, NA, 65, 70, NA, 20, 80, 35, NA, 68, 70, 60, 70, NA, 80, 98, 10, 40, 63, 75, 80, 40, 90, 100, 33, 36, 65, 78, 50)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' ts_cochran_owa(scores, groups)
#' 
#' @export
ts_cochran_owa <- function(scores, groups){
  
  dfr = na.omit(data.frame(scores, groups))
  
  #overall mean and count
  xBar = mean(dfr$scores)
  n = nrow(dfr)
  
  #means and counts and variances per category
  xBars = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=mean)$x
  ns = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=length)$x
  vars = aggregate(dfr$scores, by=list(group=dfr$groups), FUN=var)$x
  
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
  results = data.frame(statistic, df, pValue)
  
  return(results)
  
}