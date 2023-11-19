#' Post-Hoc Pairwise Student T
#' @description 
#' This function performs pairwise independent samples Student t tests, for use after a one-way ANOVA, to determine which categories significantly differ from each other.
#' 
#' It differs slightly in the calculation of the standard error, than the version used by using ph_pairwise_is(nomField, scaleField, isTest = "student"). This version appears to be producing the same results as SPSS shows, when using a Bonferroni correction. SPSS refers to Winer (1962) for their procedures.
#' 
#' A simple Bonferroni correction is also applied.
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' 
#' @returns 
#' A dataframe with:
#' \item{category 1}{the first category in the pair}
#' \item{category 2}{the second category in the pair}
#' \item{n1}{sample size of first category}
#' \item{n2}{sample size of second category}
#' \item{mean 1}{arithmetic mean of scores in first category}
#' \item{mean 2}{arithmetic mean of scores in second category}
#' \item{sample diff.}{difference between the two arithmetic means}
#' \item{hyp diff.}{the hypothesized difference}
#' \item{statistic}{the test-statistic}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the unadjusted p-value (significance)}
#' \item{adj. p-value}{the Bonferroni adjusted p-values}
#' \item{test}{description of test used}
#' 
#' @details
#' The formula used:
#' \deqn{t_{1,2} = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{MS_w \times\left(\frac{1}{n_1}+ \frac{1}{n_2}\right)}}}
#' \deqn{df_w = n - k}
#' \deqn{sig. = 2\times\left(1 - T\left(\left|t_{1,2}\right|, df_w\right)\right)}
#' 
#' With:
#' \deqn{MS_w = \frac{SS_w}{df_w}}
#' \deqn{SS_w = \sum_{j=1}^k \sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}
#' \deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{x_{i,j}}, the i-th score in category j
#' \item \eqn{n}, the total sample size
#' \item \eqn{n_j}, the number of scores in category j
#' \item \eqn{k}, the number of categories
#' \item \eqn{\bar{x}_j}, the mean of the scores in category j
#' \item \eqn{MS_w}, the mean square within
#' \item \eqn{SS_w}, the sum of squares of within (sum of squared deviation of the mean)
#' \item \eqn{df_w}, the degrees of freedom of within
#' }
#' 
#' The Bonferroni adjustment is simply:
#' \deqn{p_{adj} = \min \left(p \times n_{comp}, 1\right)}
#' \deqn{n_{comp} = \frac{k\times\left(k-1\right)}{2}}
#' 
#' *Symbols used:*
#' 
#' \itemize{
#'  \item \eqn{n_{comp}}, number of comparisons (pairs)
#'  \item \eqn{k}, number of categories
#' }
#' 
#' 
#' @references
#' Winer, B. J. (1962). *Statistical principles in experimental design*. McGraw Hill.
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_pairwise_t <- function(nomField, scaleField, categories=NULL){
  
  dfr = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$nomField %in% categories),]}
  
  colNames = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=sum)$group
  nj = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=length)$x
  sj = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=sum)$x
  mj = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=mean)$x 
  
  k <- length(colNames)
  n = sum(nj)
  m = mean(dfr$scaleField)  
  
  ssb = sum(nj*(mj-m)**2)
  sst = var(dfr$scaleField)*(n-1)
  ssw = sst - ssb
  
  dfw = n - k
  msw = ssw/dfw  
  
  ncomp = k * (k - 1) / 2
  
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 13))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      selCats = c(colNames[i], colNames[j])
      res[resRow,1] = colNames[i]
      res[resRow,2] = colNames[j]
      res[resRow, 3] = nj[i]
      res[resRow, 4] = nj[j]
      res[resRow, 5] = mj[i]
      res[resRow, 6] = mj[j]
      res[resRow, 7] = res[resRow, 5] - res[resRow, 6]
      res[resRow, 8] = 0
      
      sej = (msw * (1 / res[resRow, 3] + 1 / res[resRow, 4]))**0.5
      tVal = res[resRow, 7]/sej
      res[resRow, 9] = tVal
      
      res[resRow, 10] = dfw
      
      pValue = 2*(1-pt(abs(tVal), dfw))
      res[resRow, 11] = pValue
      res[resRow, 12] = pValue*ncomp
      if (res[resRow,12] > 1){
        res[resRow,12] = 1}
      
      res[resRow, 13] = "Winer pairwise t"
      
      resRow = resRow + 1
    }
  }  
  
  colnames(res) = c("category 1", "category 2", "n1", "n2", "mean 1", "mean 2", "sample diff.", "hyp diff.", "statistic", "df", "p-value", "adj. p-value", "test")  
  return (res)
  
}