#' Post-Hoc Conover-Iman Test
#' @description 
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. A simple Bonferroni adjustment is also made for the multiple testing.
#' 
#' Other post-hoc tests that could be considered are Dunn, Nemenyi, Steel-Dwass, a pairwise Mann-Whitney U, or pairwise Mood-Median.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' 
#' @returns
#' A dataframe with:
#' \item{cat. 1}{one of the two categories being compared}
#' \item{cat. 2}{second of the two categories being compared}
#' \item{n1}{number of cat. 1. cases in comparison}
#' \item{n2}{number of cat. 2 cases in comparison}
#' \item{mean rank 1}{mean rank of cases in cat. 1, based on all cases (incl. categories not in comparison)}
#' \item{mean rank 2}{mean rank of cases in cat. 2, based on all cases (incl. categories not in comparison)}
#' \item{statistic}{the t-value of the test}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the p-value (significance)}
#' 
#' @details
#' The formula used is (Conover & Iman, 1979, p. 11):
#' \deqn{t_{1,2} = \frac{\bar{r}_1 - \bar{r}_2}{\sqrt{S^2\times\frac{n-1-T}{n-k}\times\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}}}
#' \deqn{df = n - k}
#' \deqn{sig. = 1 - T\left(\left|t_{1,2}\right|, df\right)}
#' 
#' With:
#' \deqn{S^2=\frac{\sum_{j=1}^k \sum_{i=1}^{n_j} r_{i,j}^2 - \frac{n\times\left(n+1\right)^2}{4}}{n-1}}
#' \deqn{T = \frac{\sum_{i=1}^k \frac{R_i^2}{n_i} - \frac{n\times\left(n+1\right)^2}{4}}{S^2}}
#' \deqn{R_i = \sum_{j=1}^{n_i} r_{i,j}}
#' 
#' Note that \eqn{S^2, T, k, n} are all based on all scores, including those not in the selected pair.
#' 
#' The formula can also be found in Conover (1980, pp. 230-231).
#' 
#' *Symbols used*
#' 
#' \itemize{
#'  \item \eqn{k}, the number of categories
#'  \item \eqn{n_i}, the number of scores in category i
#'  \item \eqn{r_{i,j}}, the rank of the j-th score in category i using all original scores (incl. those not in the comparison).
#'  \item \eqn{R_i}, the sum of the ranks in category i
#'  \item \eqn{\bar{r}_i}, the average of the ranks in category i, using all original scores (incl. those not in the comparison).
#'  \item \eqn{T\left(\dots\right)}, the cumulative distribution function of the Student t distribution.
#'  }
#' 
#' @references
#' Conover, W. J. (1980). *Practical nonparametric statistics* (2nd ed.). Wiley. 
#' 
#' Conover, W. J., & Iman, R. L. (1979). *On multiple-comparisons procedures* (LA-7677-MS; pp. 1-14). Los Alamos Scientific Laboratory.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_conover_iman <- function(catField, ordField, 
                            categories=NULL, levels=NULL){
  
  #set preset values
  tiesPresent = FALSE
  nSame = TRUE
  
  if (!is.null(levels)){
    myFieldOrd = factor(ordField, ordered = TRUE, levels = levels)
    ordField = as.numeric(myFieldOrd)
  }
  
  dfr = na.omit(data.frame(ordField, catField))
  
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$catField %in% categories),]}
  
  #add ranks to dataframe
  dfr["r"] = rank(dfr$ordField)
  
  #sample size
  n = nrow(dfr)
  
  #Sum of ranks and count for each category
  Rc = aggregate(dfr$r, by=list(group=dfr$catField), FUN=sum)$x
  nc = aggregate(dfr$r, by=list(group=dfr$catField), FUN=length)$x
  grNames = aggregate(dfr$r, by=list(group=dfr$catField), FUN=sum)$group
  
  k = length(nc)
  
  T = sum(Rc**2/nc)
  srj2 = sum(dfr$r**2)
  
  ff = n * (n + 1)**2 / 4
  s2 = (srj2 - ff) / (n - 1)
  T = (T - ff) / s2
  ff = s2 * (n - 1 - T) / (n - k)
  df = n - k
  
  #The Tests
  ncomp = k * (k - 1) / 2
  res = data.frame(matrix(nrow = ncomp, ncol = 9))
  resRow = 1
  
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      n1 = nc[i]
      n2 = nc[j]
      m1 = Rc[i] / n1
      m2 = Rc[j] / n2
      d = m1 - m2
      
      Var = ff * (1 / n1 + 1 / n2)      
      se = (Var)**0.5
      tVal = d/se
      p = 2*(1 - pt(abs(tVal), df))
      
      res[resRow, 1] = grNames[i]
      res[resRow, 2] = grNames[j]
      res[resRow, 3] = n1
      res[resRow, 4] = n2
      res[resRow, 5] = m1
      res[resRow, 6] = m2
      res[resRow, 7] = tVal
      res[resRow, 8] = df
      res[resRow, 9] = p
      if (res[resRow, 9] > 1){res[resRow, 9] = 1}
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("cat. 1", "cat. 2", "n1", "n2", "mean rank 1", "mean rank 2", "statistic", "df", "p-value")
  return (res)
  
}



