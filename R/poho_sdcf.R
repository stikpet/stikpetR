#' Post-Hoc Steel-Dwass-Critchlow-Fligner Test
#' @description 
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. 
#' 
#' Other post-hoc tests that could be considered are Dunn, Nemenyi, Conover, a pairwise Mann-Whitney U, or pairwise Mood-Median.
#' 
#' Unlike the Dunn, Nemenyi and Conover-Iman test, this test re-calculates the mean ranks for each pair, using only the scores from the two categories. 
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
#' \item{statistic}{the q-value of the test}
#' \item{std. statistic}{the standardized q value}
#' \item{p-value}{the p-value (significance)}
#' 
#' 
#' @details
#' The formula used (Hollander & Wolfe, 1999, p. 241):
#' \deqn{q_{1,2} = \frac{\left| R_1 - E_1\right|}{\sqrt{\sigma^2}}}
#' 
#' With:
#' \deqn{R_1 = \sum_{i=1}^{n_1} r_{i,1}}
#' \deqn{n_{1,2} = n_1 + n_2}
#' \deqn{E_1 = \frac{n_1\times\left(n_{1,2} + 1\right)}{2}}
#' \deqn{\sigma^2 = \frac{n_1\times n_2}{12}\times\left(n_{1,2} + 1 - \frac{T}{n_{1,2}-1}\right)}
#' \deqn{T = \sum t_j^3 - t_j}
#' 
#' The p-value is then determined using (Critchlow & Fligner, 1991, p. 131):
#' \deqn{sig. = 1 - Q\left(q_{1,2}, k, df=\infty\right)}
#' 
#' Note that while looking at the R-code for this, posted by Shigenobu (n.d.), who references Nagata and Yoshida (1997), an alternative but same result equation for the variance can be used:
#' \deqn{\sigma^2 = \frac{n_1\times n_2}{n_{1,2}\times\left(n_{1,2} - 1\right)}\times\left(\sum_{i=1}^{n_1} r_{i,1}^2 + \sum_{i=1}^{n_2} r_{i,2}^2 - \frac{n_{1,2}\times\left(n_{1,2}+1\right)^2}{4}\right)}.
#' 
#' Steel (1960) and Dwass (1960) independently derived the basics for this test. Critchlow and Fligner (1991) added the case for larger samples using the Tukey Range Distribution, and in Hollander and Wolfe (1999) the version used here can be found, which includes a ties correction.
#' 
#' *Symbols used*
#' 
#' \itemize{
#'  \item \eqn{k}, the number of categories
#'  \item \eqn{t_j}, the frequency of the j-th unique rank.
#'  \item \eqn{n_i}, the number of scores in category i
#'  \item \eqn{r_{i,j}}, the rank of the j-th score in category i using only the scores from the two categories in the comparison.
#' \item \eqn{Q\left(\dots\right)}, the cumulative distribution function of the standardized range distribution.
#' }
#' 
#' @references
#' Critchlow, D. E., & Fligner, M. A. (1991). On distribution-free multiple comparisons in the one-way analysis of variance. *Communications in Statistics - Theory and Methods, 20*(1), 127–139. doi:10.1080/03610929108830487
#' 
#' Dwass, M. (1960). *Some k-sample rank-order tests*. In I. Olkin, S. G. Ghurye, W. Hoeffding, W. G. Madow, & H. B. Mann (Eds.), Contributions to probability and statistics; Essays in honor of Harold Hotelling. Stanford University Press.
#' 
#' Hollander, M., & Wolfe, D. A. (1999). *Nonparametric statistical methods* (2nd ed.). Wiley.
#' 
#' Nagata, Y., & Yoshida, M. (1997). *The Basics of Multiple Comparisons in Statistics*. Scientist Co.
#' 
#' Shigenobu. (2004, July 28). Multiple comparisons using the Steel-Dwass method. http://aoki2.si.gunma-u.ac.jp/R/Steel-Dwass.html
#' 
#' Steel, R. G. D. (1960). A rank sum test for comparing all pairs of treatments. Technometrics, 2(2), 197–207. doi:10.1080/00401706.1960.10489894
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_sdcf <- function(catField, ordField, 
                    categories=NULL, levels=NULL){
  
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
  
  grNames = aggregate(dfr$r, by=list(group=dfr$catField), FUN=sum)$group  
  k = length(grNames)
  
  #The Tests
  ncomp = k * (k - 1) / 2
  res = data.frame(matrix(nrow = ncomp, ncol = 9))
  resRow = 1
  
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      res[resRow, 1] = grNames[i]
      res[resRow, 2] = grNames[j]
      selCats = c(res[resRow, 1], res[resRow, 2])
      dfPair = dfr[dfr$catField %in% selCats, ]
      
      nPair = nrow(dfPair)
      dfPair["r"] = rank(dfPair$ordField)
      
      #Ties correction
      #Frequencies of each rank
      tTab = table(dfPair["r"])    
      
      #Ties adjustment
      t = sum(tTab**3 - tTab)
      
      #Sum of ranks and count for each category
      Rc = aggregate(dfPair$r, by=list(group=dfPair$catField), FUN=sum)$x
      nc = aggregate(dfPair$r, by=list(group=dfPair$catField), FUN=length)$x
      
      n1 = nc[1]
      n2 = nc[2]
      
      e = n1 * (nPair + 1) / 2
      s2 = n1 * n2 / 12 * (nPair + 1 - t / (nPair * (nPair - 1)))
      q = (Rc[1] - e) / (s2**0.5)
      
      res[resRow, 3] = n1
      res[resRow, 4] = n2
      res[resRow, 5] = Rc[1]/n1
      res[resRow, 6] = Rc[2]/n2
      res[resRow, 7] = q
      res[resRow, 8] = q*2**0.5
      res[resRow, 9] = 1 - ptukey(abs(q*2**0.5), nmeans=k, df=Inf)
      
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("cat. 1", "cat. 2", "n1", "n2", "mean rank 1", "mean rank 2", "statistic", "std. statistic", "p-value")
  return (res)
  
}