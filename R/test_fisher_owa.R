#' Fisher/Classic One-Way ANOVA / F-Test
#' @description 
#' Tests if the means (averages) of each category could be the same in the population.
#' 
#' If the p-value is below a pre-defined threshold (usually 0.05), the null hypothesis is rejected, and there are then at least two categories who will have a different mean on the scaleField score in the population.
#' 
#' There are quite some alternatives for this, the stikpet library has Welch, James, Box, Scott-Smith, Brown-Forsythe, Alexander-Govern, Mehrotra modified Brown-Forsythe, Hartung-Agac-Makabi, Özdemir-Kurt and Wilcox as options. See the notes for some discussion on the differences.
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField

#' @returns 
#' A dataframe with an ANOVA table showing:
#' \item{variance}{which variance is shown in that row}
#' \item{SS}{sum of squared deviations from the mean}
#' \item{df}{degrees of freedom}
#' \item{MS}{the mean square}
#' \item{F}{the F-statistic value}
#' \item{pValue}{the significance (p-value)}
#' 
#' @details 
#' The formula used:
#' \deqn{F_F = \frac{df_w\times SS_b}{df_b\times SSw}}
#' \deqn{df_b = k - 1}
#' \deqn{df_w = n - k}
#' \deqn{sig. = 1 - F\left(F_F, df_b, df_w\right)}
#' With:
#' \deqn{SS_b = \sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}
#' \deqn{SS_w = SS_t - SS_b}
#' \deqn{SS_t = \sum_{j=1}^k \sum_{i=1}^{n_j}\left(x_{i,j} - \bar{x}\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' \deqn{\bar{x} = \frac{\sum_{j=1}^k n_j\times\bar{x}_j }{n} = \frac{\sum_{j=1}^k \sum_{i=1}^{n_j} x_{i,j}}{n}}
#' \deqn{n = \sum_{j=1}^k n_j}
#' 
#' Alternative format of the F-statistic equation (but the same result):
#' \deqn{F_F = \frac{MS_b}{MSw}}
#' 
#' *Symbols*
#' \itemize{
#' \item \eqn{x_{i,j}} the i-th score in category j
#' \item \eqn{n} the total sample size
#' \item \eqn{n_j} the number of scores in category j
#' \item \eqn{k} the number of categories
#' \item \eqn{\bar{x}_j} the mean of the scores in category j
#' \item \eqn{SS_i} the sum of squares of i (sum of squared deviation of the mean)
#' \item \eqn{df_i} the degrees of freedom of i
#' \item \eqn{b} is between = factor = treatment = model
#' \item \eqn{w} is within = error (the variability within the groups)
#' }
#' 
#' Note that the Fisher-Pitman test (Pitman, 1937a, 1937b, 1938) uses a different approach but will lead to the same result. 
#' 
#' I'm not fully sure what the original source is for the Fisher test, but likely either of his
#' sources from 1918, 1921, 1925 or 1935.
#' 
#' **Choosing a test**
#' 
#' The classic/Fisher one-way ANOVA assumes the data is normally distributed and that the variances in each group are the same in the population (homoscedasticity). Many have tried to cover the situations when one or both of these conditions are not met.
#' 
#' Delacre et al. (2019) recommend to use the Welch ANOVA instead of the classic and Brown-Forsythe versions. How2stats (2018) give a slightly different recommendation based on Tomarken and Serlin (1986). They agree that usually the Welch ANOVA is preferred of the classic version, but if the average sample size is below six to still use the Brown-Forsythe.
#' 
#' The researchers in the previous paragraph did not take into consideration other approaches. A few comments found on those other methods.
#' 
#' According to Hartung et al. (2002, p. 225) the Cochran test is the standard test in meta-analysis, but should not be used, since it is always too liberal.
#' 
#' Schneider and Penfield (1997) looked at the Welch, Alexander-Govern and the James test (they ignored the Brown-Forsythe since they found it to perform worse than Welch or James), and concluded: “Under variance heterogeneity, Alexander-Govern’s approximation was not only comparable to the Welch test and the James second-order test but was superior, in certain instances, when coupled with the power results for those tests” (p. 285).
#' 
#' Cavus and Yazici (2020) compared many different tests. They showed that the Brown-Forsythe, Box correction, Cochran, Hartung-Agac-Makabi adjusted Welch, and Scott-Smith test, all do not perform well, compared to the Asiribo-Gurland correction, Alexander-Govern test, Özdemir-Kurt B2, Mehrotra modified Brown-Forsythe, and Welch.
#' 
#' I only came across the Johansen test in Algina et. al. (1991) and it appears to give the same results as the Welch test.
#' 
#' In my experience the one-way ANOVA is widely known and often discussed in textbooks. The Welch anova is gaining popularity. The Brown-Forsythe is already more obscure and some confuse it with the Brown-Forsythe test for variances. The James test and the Alexander-Govern are perhaps the least known and the Johansen even less than that (at least they were for me). So, although the Alexander-Govern test might be preferred over the Welch test, some researchers prefer to use a more commonly used test than a more obscure version. In the end it is up to you to decide on what might be the best test, and also depending on the importance of your research you might want to investigate which test fits your situation best, rather than taking my word for it.
#' 
#' Besides these, there are more methods, some using simulation (bootstrapping) (see Cavus and Yazici (2020) for a few of them), others using different techniques (see Yiğit and Gökpinar (2010) for a few more methods not in here).
#' 
#' @references 
#' Algina, J., Oshima, T. C., & Tang, K. L. (1991). Robustness of Yao’s, James’, and Johansen’s Tests under variance-covariance heteroscedasticity and nonnormality. *Journal of Educational Statistics, 16*(2), 125–139. doi:10.2307/1165116
#' 
#' Cavus, M., & Yazıcı, B. (2020). Testing the equality of normal distributed and independent groups’ means under unequal variances by doex package. *The R Journal, 12*(2), 134. doi:10.32614/RJ-2021-008
#' 
#' Delacre, M., Leys, C., Mora, Y. L., & Lakens, D. (2019). Taking parametric assumptions seriously: Arguments for the use of Welch’s F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology, 32*(1), 1–12. doi:10.5334/irsp.198
#' 
#' Fisher, R. A. (1921). On the “probable error” of a coefficient of correlation deduced from a small sample. *Metron, 1*, 3–32.
#' 
#' Hartung, J., Argaç, D., & Makambi, K. H. (2002). Small sample properties of tests on homogeneity in one-way anova and meta-analysis. *Statistical Papers, 43*(2), 197–235. doi:10.1007/s00362-002-0097-8
#' 
#' how2stats (Director). (2018, June 11). Welch’s F-test vs Brown-Forsythe F-test: Which Should You Use and When? https://youtu.be/jteKmatBgF8
#' 
#' Schneider, P. J., & Penfield, D. A. (1997). Alexander and Govern’s approximation: Providing an alternative to ANOVA under variance heterogeneity. *The Journal of Experimental Education, 65*(3), 271–286. doi:10.1080/00220973.1997.9943459
#' 
#' Tomarken, A. J., & Serlin, R. C. (1986). Comparison of ANOVA alternatives under variance heterogeneity and specific noncentrality structures. *Psychological Bulletin, 99*(1), 90–99. doi:10.1037/0033-2909.99.1.90
#' 
#' Yiğit, E., & Gökpinar, F. (2010). A simulation study on tests for one-way ANOVA under the unequal variance assumption. Communications, Faculty Of Science, University of Ankara, 15–34. doi:10.1501/Commua1_0000000660
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_fisher_owa <- function(nomField, scaleField, categories=NULL){
  dfr = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$nomField %in% categories),]}
  
  #overall mean and count
  xBar = mean(dfr$scaleField)
  n = nrow(dfr)
  
  #means and counts per category
  xBars = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=mean)$x
  ns = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=length)$x
  
  #number of categories
  k = length(ns)
  
  SSb = sum(ns*(xBars - xBar)**2)
  SSt = var(dfr$scaleField)*(n - 1)
  SSw = SSt - SSb
  
  df1 = k - 1
  df2 = n - k
  
  msb = SSb/df1
  msw = SSw/df2
  
  Fvalue = msb/msw
  
  
  pValue = 1 - pf(Fvalue, df1, df2)
  
  res = data.frame(matrix(nrow = 3, ncol = 6))  
  res[1,1] = "between"
  res[2,1] = "within"
  res[3,1] = "total"
  
  res[1,2] = SSb
  res[2,2] = SSw
  res[3,2] = SSt
  
  res[1,3] = df1
  res[2,3] = df2
  res[3,3] = n-1
  
  res[1,4] = msb
  res[2,4] = msw
  res[3,4] = NA
  
  res[1,5] = Fvalue
  res[2,5] = NA
  res[3,5] = NA
  
  res[1,6] = pValue
  res[2,6] = NA
  res[3,6] = NA
  
  colnames(res) = c("variance", "SS", "df", "MS", "F", "p-value")  
  
  return(res)
  
}