#' Mood Median Test
#' @description 
#' This test looks if the median from different categories would be the same in the population. If not, at least one is different then at least one other category. A Kruskal-Wallis test (see ts_kruksal_wallis()) is very similar but checks the average ranks instead of median.
#' 
#' The test only looks at the number of scores above the overall median and those that are equal or below. A cross table is made with each category and the numbers below and above the overall median. From this table a test of independence can be used.
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param test string, optional. the test of independence to use. Default is "pearson". Other options are "pearson", "fisher", "freeman-tukey", "g", "mod-log", "neyman", "power"
#' @param cc : string, optional. method for continuity correction. Either NULL (default), "yates", "pearson", "williams"
#' @param lambd float or string, optional. either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3. Only applies to Power Divergence test. Other options include float, "cressie-read", "likelihood-ratio", "mod-log", "pearson", "freeman-tukey", "neyman"
#' 
#' @returns 
#' A dataframe with the results of the specified test.
#' 
#' @details 
#' The Mood Median test creates a 2xk cross table, with k being the number of categories. The two rows are one for the number of scores in that category that are above the overall median, and the second row the number of scores in that category that are equal or below the overall median.
#' 
#' A chi-square test of independence on this cross table can then be performed. There are quite some different options for this:
#' 
#' \itemize{
#' \item "pearson", will perform a Pearson chi-square test of independence using the ts_pearson_ind() function.
#' \item "fisher", will perform a Fisher exact test using the ts_fisher() function, but only if there are 2 categories, if there are more the test will be set to "pearson"
#' \item "freeman-tukey", will perform a Freeman-Tukey test of independence using the ts_freeman_tukey_ind() function
#' \item "g", will perform a G test of independence using the ts_g_ind() function
#' \item "mod-log", will perform a Mod-Log Likelihood test of independence using the ts_mod_log_likelihood_ind() function
#' \item "neyman", will perform a Neyman test of independence using the ts_neyman_ind() function
#' \item "power", will perform a Power Divergence test of independence using the ts_powerdivergence_ind() function.
#' }
#' 
#' The formula using the default Pearson test is:
#' \deqn{\chi_{M}^2 = \sum_{i=1}^2 \sum_{j=1}^k \frac{\left(F_{i,j}-E_{i,j}\right)^2}{E_{i,j}}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_{M}^2, df\right)}
#' 
#' With:
#' \deqn{E_{i,j} = \frac{R_i \times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^k F_{i,j}}
#' \deqn{C_j = \sum_{i=1}^2 F_{i,j}}
#' \deqn{n = \sum_{i=1}^2 \sum_{j=1}^k F_{i,j} = \sum_{i=1}^2 R_i = \sum_{j=1}^k C_j}
#' 
#' The original source for the formula is most likely Mood (1950), but the ones shown are based on Brown and Mood (1951).
#' 
#' *Symbols used:*
#' 
#' \itemize{
#' \item \eqn{k}, the number of categories (columns)
#' \item \eqn{F_{1,j}}, the number of scores is category j that are above the overall median
#' \item \eqn{F_{2,j}}, the number of scores is category j that are equal or below the overall median
#' \item \eqn{E_{i,j}}, the expected count in row i and column j.
#' \item \eqn{R_i}, the row total of row i 
#' \item \eqn{C_j}, the column total of column j
#' \item \eqn{n}, the overall total.
#' \item \eqn{df}, the degrees of freedom
#' \item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function of the chi-square distribution.
#' }
#' 
#' @references
#' 
#' Brown, G. W., & Mood, A. M. (1951). On median tests for linear hypotheses. Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability, 2, 159–167.
#' 
#' Mood, A. M. (1950). *Introduction to the theory of statistics*. McGraw-Hill.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'
#' @export
ts_mood_median <- function(catField, ordField, 
                            categories=NULL, levels=NULL, 
                            test="pearson", 
                            cc=c(NULL,"yates", "pearson", "williams"), 
                            lambd=2/3){
  #default for cc to none
  if (length(cc)>1) {cc = NULL}
  
  if (!is.null(levels)){
    myFieldOrd = factor(ordField, ordered = TRUE, levels = levels)
    ordField = as.numeric(myFieldOrd)
  }
  
  datFrame = na.omit(data.frame(catField, ordField))
  
  #replace categories if provided
  if (!is.null(categories)){
    datFrame = datFrame[(datFrame$catField %in% categories),]}
  
  med = median(datFrame$ordField)
  
  datTable = table(datFrame)
  k = nrow(datTable) #the number of groups
  
  obs = matrix(1, nrow=2, ncol=k)
  for (i in 1:k) {
    obs[1,i] = sum(datFrame$ordField[datFrame$catField == row.names(datTable)[i]] > med)
    obs[2,i] = sum(datFrame$ordField[datFrame$catField == row.names(datTable)[i]] <= med)
  }
  
  # New dataframe with simply above or below or equal.
  catArr = c()
  ordArr = c()
  arrRow = 1
  for (j in 1:k){
    for (i in 1:2){
      for (sc in 1:obs[i, j]){
        if (obs[i,j]>0){
          catArr[arrRow] = colnames(datTable)[j]
          ordArr[arrRow] = i+1
          arrRow = arrRow + 1
        }
      }
    }
  }
  
  #now for the test
  if (test=="fisher"){
    if (k>2){test = "pearson"}
    else{
      res = ts_fisher(catArr, ordArr)}
  }
  
  if (test=="freeman-tukey"){
    res = ts_freeman_tukey_ind(catArr, ordArr, cc=cc)}
  else if (test=="g"){
    res = ts_g_ind(catArr, ordArr, cc=cc)}
  else if (test=="mod-log"){
    res = ts_mod_log_likelihood_ind(catArr, ordArr, cc=cc)}
  else if (test=="neyman"){
    res = ts_neyman_ind(catArr, ordArr, cc=cc)}
  else if (test=="pearson"){
    res = ts_pearson_ind(catArr, ordArr, cc=cc)}
  else if (test=="power"){
    res = ts_powerdivergence_ind(catArr, ordArr, cc=cc, lambd=lambd)}
  
  
  return(res)
  
}