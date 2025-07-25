#' Post-Hoc Pairwise Independent Samples Test
#' @description 
#' This function can perform various pairwise independent samples tests, for use after a one-way ANOVA, to determine which categories significantly differ from each other.
#' 
#' A simple Bonferroni correction is also applied.
#' 
#' The independent samples tests that can be used are:
#' \itemize{
#' \item Student t, see ts_student_t_is() for details. An alternative version for this is available by using the ph_pairwise_t() function.
#' \item Welch t, see ts_welch_t_is() for details
#' \item Trimmed Mean / Yuen, see ts_trimmed_mean_is() for details
#' \item Z, see ts_z_is() for details
#' }
#' 
#' @param nomField the groups variable
#' @param scaleField the numeric scores variable
#' @param categories vector, optional. the categories to use from catField
#' @param isTest string, optional. The independent samples test to use. Either "student" (default), "welch", "trimmed", "yuen", "z"
#' @param trimProp float, optional. The trim proportion to use, if applicable. Default is 0.1.
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
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_pairwise_is <- function(nomField, scaleField, categories=NULL, isTest = "student", trimProp = 0.1){
  
  dfr = na.omit(data.frame(scaleField, nomField))
  #replace categories if provided
  if (!is.null(categories)){
    dfr = dfr[(dfr$nomField %in% categories),]}
  
  colNames = aggregate(dfr$scaleField, by=list(group=dfr$nomField), FUN=sum)$group
  k <- length(colNames)
  
  ncomp = k * (k - 1) / 2
  
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 13))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      selCats = c(colNames[i], colNames[j])
      res[resRow,1] = colNames[i]
      res[resRow,2] = colNames[j]
      
      sel2cat = c(colNames[i], colNames[j])
      if (isTest == "student"){
        isRes = ts_student_t_is(nomField, scaleField, sel2cat)}
      else if (isTest == "welch"){
        isRes = ts_welch_t_is(nomField, scaleField, sel2cat)}
      else if (isTest == "trimmed"){
        isRes = ts_trimmed_mean_is(nomField, scaleField, sel2cat, trimProp=trimProp, se="yuen-dixon")}
      else if (isTest == "yuen"){
        isRes = ts_trimmed_mean_is(nomField, scaleField, sel2cat, trimProp=trimProp, se="yuen")}
      else if (isTest == "z"){
        isRes = ts_z_is(nomField, scaleField, sel2cat)}
      
      res[resRow, 3] = isRes[1,1]
      res[resRow, 4] = isRes[1,2]
      res[resRow, 5] = isRes[1,3]
      res[resRow, 6] = isRes[1,4]
      res[resRow, 7] = isRes[1,5]
      res[resRow, 8] = isRes[1,6]
      res[resRow, 9] = isRes[1,7]
      
      if (isTest == "z"){
        res[resRow, 10] = NA
        res[resRow, 11] = isRes[1,8]}
      else{
        res[resRow, 10] = isRes[1,8]
        res[resRow, 11] = isRes[1,9]}
      
      res[resRow, 12] = res[resRow,11] * ncomp
      if (res[resRow,12] > 1){
        res.iloc[resRow,12] = 1}
      
      if (isTest == "z"){
        res[resRow, 13] = isRes[1,9]}
      else{
        res[resRow, 13] = isRes[1,10]}
      
      resRow = resRow + 1
    }
  }  
  
  colnames(res) = c("category 1", "category 2", "n1", "n2", "mean 1", "mean 2", "sample diff.", "hyp diff.", "statistic", "df", "p-value", "adj. p-value", "test")  
  return (res)
  
}



