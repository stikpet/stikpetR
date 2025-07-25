#' Post-Hoc Nemenyi Test
#' @description 
#' This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).
#' 
#' The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal. 
#' 
#' Pohlert (2016) mentions the exact version should only be used if there are no ties, and suggest to use a chi-square alternative in case of ties. This is referred to by Zaiontz (n.d.-b) as the Schaich-Hamerle test (1984).
#' 
#' The ties correction is taken from Pohlert (2016).
#' 
#' Other post-hoc tests that could be considered are Dunn, Steel-Dwass, Conover, a pairwise Mann-Whitney U, or pairwise Mood-Median.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param version string, optional. version of the test to use. Either "auto" (default), "exact", "sh", "sh-ties".
#' 
#' @returns
#' A dataframe with:
#' \item{cat. 1}{one of the two categories being compared}
#' \item{cat. 2}{second of the two categories being compared}
#' \item{n1}{number of cat. 1. cases in comparison}
#' \item{n2}{number of cat. 2 cases in comparison}
#' \item{mean rank 1}{mean rank of cases in cat. 1, based on all cases (incl. categories not in comparison)}
#' \item{mean rank 2}{mean rank of cases in cat. 2, based on all cases (incl. categories not in comparison)}
#' \item{se}{the standard error used}
#' \item{statistic}{the z-value of the test}
#' \item{p-value}{the p-value (significance)}
#' 
#' 
#' @details
#' The formula used (Pohlert, 2016, p. 3):
#' \deqn{q_{1,2} = \frac{\bar{r}_1 - \bar{r}_2}{\sqrt{\frac{n\times\left(n+1\right)}{24}\times\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}}}
#' \deqn{sig. = 1 - Q\left(q_{1,2}, k, df=\infty\right)}
#' 
#' A chi-square distribution can also be used. Zaiontz (n.d.) and BrightStat (n.d.) refer to this as the Schaich-Hamerle test. 
#' 
#' The formula used then changes to (Sachs, 1982, p. 549):
#' \deqn{\chi_{1,2}^2 = \frac{\left(\bar{r}_1 - \bar{r}_2\right)^2}{\frac{n\times\left(n+1\right)}{12}\times\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{1,2}^2, df\right)}
#' 
#' A ties correction found in Pohlert (2016, p. 3) adjusts this to:
#' \deqn{\chi_{1,2}^2 = \frac{\left(\bar{r}_1 - \bar{r}_2\right)^2}{\left(1-T\right)\times\frac{n\times\left(n+1\right)}{12}\times\left(\frac{1}{n_1}+\frac{1}{n_2}\right)}}
#' \deqn{T = \frac{\sum t_i^3 - t_i}{n^3 - n}}
#' 
#' The original formula is most likely from Nemenyi (1963) and the Schaich and Hamerle (1984).
#' 
#' *Symbols used*
#' 
#' \itemize{
#'  \item \eqn{k}, the number of categories
#'  \item \eqn{t_j}, the frequency of the j-th unique rank.
#'  \item \eqn{n}, the total sample size, of all scores, incl. those not in the comparison
#'  \item \eqn{n_i}, the number of scores in category i
#'  \item \eqn{r_{i,j}}, the rank of the j-th score in category i, using all original scores (incl. those not in the comparison).
#'  \item \eqn{\bar{r}_i}, the average of the ranks in category i, using all original scores (incl. those not in the comparison).
#'  \item \eqn{Q\left(\dots\right)}, the cumulative distribution function of the standardized range distribution.
#'  \item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function of the chi-square distribution.
#' }
#' 
#' @references
#' BrightStat. (n.d.). Kruskal-Wallis test. BrightStat. Retrieved October 25, 2023, from https://secure.brightstat.com/index.php?p=c&d=1&c=2&i=7
#' 
#' Nemenyi, P. (1963). *Distribution-free Multiple Comparisons*. Princeton University.
#' 
#' Pohlert, T. (2016). The pairwise multiple comparison of mean ranks package (PMCMR). https://cran.r-hub.io/web/packages/PMCMR/vignettes/PMCMR.pdf
#' 
#' Sachs, L. (1982). *Applied statistics: A handbook of techniques*. Springer-Verlag.
#' 
#' Schaich, E., & Hamerle, A. (1984). *Verteilungsfreie statistische Prufverfahren*. Springer. doi:10.1007/978-3-642-70032-3
#' 
#' Zaiontz, C. (n.d.). Schaich-Hamerle Test after KW. Real Statistics Using Excel. Retrieved October 25, 2023, from https://real-statistics.com/one-way-analysis-of-variance-anova/kruskal-wallis-test/schaich-hamerle-test/
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_nemenyi <- function(catField, ordField, 
                       categories=NULL, levels=NULL, version="auto"){
  
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
  
  #Ties correction
  #Frequencies of each rank
  tTab = table(dfr["r"])    
  
  #Ties adjustment
  t = sum(tTab**3 - tTab)    
  t = t/(n**3 - n)
  
  if (version == "auto"){
    if (!tiesPresent && nSame){
      version = "exact"}
    else if (!tiesPresent && nSame){
      version = "sh"}
    else{
      version = "sh-ties"}
  }
  
  if (version == "exact"){
    ff = n * (n + 1) / 24
    df = n - k}
  else{
    ff = n * (n + 1) / 12
    df = k - 1}
  
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
      
      if (version == "exact"){
        Var = ff * (1 / n1 + 1 / n2)}
      else if (version == "sh-ties"){
        Var = (1 - t) * ff * (1 / n1 + 1 / n2)}
      else{
        Var = ff * (1 / n1 + 1 / n2)}
      
      se = (Var)**0.5
      
      if (version == "exact"){
        stat = d / se
        p = 1 - ptukey(abs(stat), k, 100001)}
      else{
        stat = (d / se)**2
        p = pchisq(stat, df, lower.tail=FALSE)}
      
      res[resRow, 1] = grNames[i]
      res[resRow, 2] = grNames[j]
      res[resRow, 3] = n1
      res[resRow, 4] = n2
      res[resRow, 5] = m1
      res[resRow, 6] = m2
      res[resRow, 7] = se
      res[resRow, 8] = stat
      res[resRow, 9] = p
      if (res[resRow, 9] > 1){res[resRow, 9] = 1}
      resRow = resRow + 1
    }
  }
  
  colnames(res) = c("cat. 1", "cat. 2", "n1", "n2", "mean rank 1", "mean rank 2", "se", "statistic", "p-value")
  return (res)
  
}



