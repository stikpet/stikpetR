#' Post-Hoc Dunn Test
#' @description 
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. A simple Bonferroni adjustment is also made for the multiple testing.
#' 
#' Dunn (1964) describes two procedures. The first is his own, the second is one from Steel (1960) for comparison. The difference is that in Dunn's procedure the mean rank of each category is based on the scores of all categories, including those that are not being compared, while Steel's procedure re-calculates the mean rank for each category using only the scores from the two categories being compared. This later one wouls make it very similar to a pairwise Mann-Whitney U test (see ph_mann_whitney()).
#' 
#' Other post-hoc tests that could be considered are Nemenyi, Steel-Dwass, Conover, a pairwise Mann-Whitney U, or pairwise Mood-Median.
#' 
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
#' \item{statistic}{the z-value of the test}
#' \item{p-value}{the p-value (significance)}
#' \item{adj. p-value}{the Bonferroni adjusted p-value}
#' 
#' @details
#' The formula used (Dunn, 1964, p. 249):
#' \deqn{z_{1,2} = \frac{\bar{r}_1 - \bar{r}_2}{\sqrt{\sigma_m^2}}}
#' \deqn{sig. = 2\times\left(1 -\Phi\left(z\right)\right)}
#' 
#' With:
#' \deqn{\sigma_m^2 = \left(\frac{n\times\left(n+1\right)}{12} - \frac{T}{12\times\left(n - 1\right)}\right)\times\left(\frac{1}{n_1} + \frac{1}{n_2}\right)}
#' \deqn{T = \sum_{j=1}^k t_j^3 - t_j}
#' \deqn{\bar{r}_i = \frac{R_i}{n_i}}
#' \deqn{R_i = \sum_{j=1}^{n_i} r_{i,j}}
#' 
#' *Symbols used*
#' 
#' \itemize{
#' \item \eqn{k}, the number of categories
#' \item \eqn{t_j}, the frequency of the j-th unique rank.
#' \item \eqn{n_i}, the number of scores in category i
#' \item \eqn{r_{i,j}}, the rank of the j-th score in category i using all original scores (incl. those not in the comparison).
#' \item \eqn{R_i}, the sum of the ranks in category i
#' \item \eqn{\bar{r}_i}, the average of the ranks in category i, using all original scores (incl. those not in the comparison).
#' \item \eqn{\Phi\left(\dots\right)}, the cumulative distribution function of the standard normal distribution.
#' }
#' 
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics, 6*(3), 241–252. doi:10.1080/00401706.1964.10490181
#' 
#' Steel, R. G. D. (1960). A rank sum test for comparing all pairs of treatments. *Technometrics, 2*(2), 197–207. doi:10.1080/00401706.1960.10489894
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
ph_dunn <- function(catField, ordField, 
                    categories=NULL, levels=NULL){
  
  if (!is.null(levels)){
    myFieldOrd = factor(ordField, ordered = TRUE, levels = levels)
    ordField = as.numeric(myFieldOrd)
  }
  
  if (!is.null(categories)){
    if(!is.null(names(categories))){
      for (i in 1:length(categories)){catField[catField == unname(categories[i])] = names(categories)[i]}
      categories = names(categories)
    }
  }
  
  dfr = na.omit(data.frame(ordField, catField))
  
  #add ranks to dataframe
  dfr["r"] = rank(dfr$ordField)
  
  #sample size
  n = nrow(dfr)
  
  #Sum of ranks and count for each category
  Rc = aggregate(dfr$r, by=list(group=dfr$catField), FUN=sum)$x
  nc = aggregate(dfr$r, by=list(group=dfr$catField), FUN=length)$x
  grNames = aggregate(dfr$r, by=list(group=dfr$catField), FUN=sum)$group
  
  k = length(nc)
  
  #Ties correction
  #Frequencies of each rank
  tTab = table(dfr["r"])    
  #Ties adjustment
  t = sum(tTab**3 - tTab)   
  
  a = n * (n + 1) / 12 - t / (12 * (n - 1))
  
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
      Var = a * (1 / n1 + 1 / n2)
      se = (Var)**0.5
      Z = d / se
      p = 2*(1-pnorm(abs(Z)))
      
      res[resRow, 1] = grNames[i]
      res[resRow, 2] = grNames[j]
      res[resRow, 3] = n1
      res[resRow, 4] = n2
      res[resRow, 5] = m1
      res[resRow, 6] = m2
      res[resRow, 7] = Z
      res[resRow, 8] = p
      res[resRow, 9] = p * ncomp
      if (res[resRow, 9] > 1){res[resRow, 9] = 1}
      
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("cat. 1", "cat. 2", "n1", "n2", "mean rank 1", "mean rank 2", "statistic", "p-value", "adj. p-value")
  return (res)
  
}